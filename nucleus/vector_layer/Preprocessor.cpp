/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2024 Lucas Dworschak
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "Preprocessor.h"

#include <QImage>
#include <QString>

#include <CDT.h>

#include "nucleus/Raster.h"
#include <nucleus/utils/bit_coding.h>

#include "nucleus/utils/rasterizer.h"

#include <mapbox/vector_tile.hpp>

#include "constants.h"

namespace nucleus::vector_layer {

// TODO here:
// - integrate style buffer into gl engine / shader
// - write style index to buffer per triangle
// - open the floodgates and visualize other polygons
// - probably improve shader so that everything is visualized correctly

// - ok if we allow all polygons we are over our size for each individual tile for some tiles...
// -> possible solution split each tile into smaller portions (to prevent too much empty data), but allow overflow -> tile x,y,z data is available in texture index 400, 401, and 408
// -> here we sucessfully split it up (at least it should be split up) -> next gl engine upload the split up thing than interpret it on the shader correctly

GpuVectorLayerTile preprocess(tile::Id id, const QByteArray& vector_tile_data, const Style& style)
// GpuVectorLayerTile preprocess(tile::Id, const QByteArray& vector_tile_data, const Style&)
{
    if (vector_tile_data.isEmpty())
        return {};

    // std::cout << "entry:\t" << id.zoom_level << "\t";

    // DEBUG polygons
    // const std::vector<std::vector<glm::vec2>> triangle_points = { { glm::vec2(10.5 / 64.0 * constants::grid_size, 30.5 / 64.0 * constants::grid_size),
    //     glm::vec2(30.5 / 64.0 * constants::grid_size, 10.5 / 64.0 * constants::grid_size),
    //     glm::vec2(50.5 / 64.0 * constants::grid_size, 50.5 / 64.0 * constants::grid_size) } };
    const std::vector<unsigned int> style_indices = { 1 };

    auto triangle_points = details::parse_tile(id, vector_tile_data, style);

    // // TODO somehow parse the data to lines and triangles

    auto processed_triangles = details::preprocess_triangles(triangle_points, style_indices);

    auto tile = create_gpu_tile(processed_triangles, processed_triangles);

    return tile;
}

} // namespace nucleus::vector_layer
namespace nucleus::vector_layer::details {

// PointCollectionVec2 parse_tile(tile::Id id, const QByteArray& vector_tile_data, const Style& style)
std::vector<PolygonData> parse_tile(tile::Id, const QByteArray& vector_tile_data, const Style&)
{
    const auto d = vector_tile_data.toStdString();
    const mapbox::vector_tile::buffer tile(d);

    bool first_layer = true;
    float grid_scale = 1.0; // scale between [0-grid_size]

    std::vector<PolygonData> polygons;
    std::vector<size_t> polygon_styles;
    std::vector<PolygonData> lines;

    for (const auto& layer_name : tile.layerNames()) {
        // std::cout << layer_name << std::endl;

        // if (layer_name != "NUTZUNG_L15_12")
        if (layer_name != "GEWAESSER_F_GEWF")
            continue; // DEBUG

        // TODO enable again
        // size_t style_index = style.layer_style_index(layer_name, id.zoom_level);
        // if (style_index == -1ul) // no style found -> we do not visualize it
        //     continue;

        const mapbox::vector_tile::layer layer = tile.getLayer(layer_name);
        std::size_t feature_count = layer.featureCount();

        if (first_layer) // TODO extent should be the same along all layers -> verify this
        {
            first_layer = false;
            grid_scale = float(constants::grid_size) / layer.getExtent();
        }

        for (std::size_t i = 0; i < feature_count; ++i) {
            const auto feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);
            auto props = feature.getProperties();

            // TODO style_index decision should also include props._symbols (and other variables defined there)

            if (feature.getType() == mapbox::vector_tile::GeomType::POLYGON) {
                PointCollectionVec2 geom = feature.getGeometries<PointCollectionVec2>(grid_scale);

                // for (size_t j = 0; j < geom.size(); ++j) {
                //     if (geom[j][0] == geom[j][geom[j].size() - 1]) {
                //         // duplicate vertex detected -> has to be removed
                //         geom[j].pop_back();
                //     }
                // }

                // construct edges from the current polygon and move them to the collection
                PolygonData p;
                for (size_t j = 0; j < geom.size(); ++j) {
                    auto edges = nucleus::utils::rasterizer::generate_neighbour_edges(geom[j], p.vertices.size());
                    p.edges.insert(p.edges.end(), edges.begin(), edges.end());
                    p.vertices.insert(p.vertices.end(), geom[j].begin(), geom[j].end());
                }
                polygons.push_back(p);

                // concat
                // polygons.reserve(polygons.size() + geom.size()); // calling reserve here significantly impacts performance -> see if it is possible to call it outside of the loop
                // polygons.insert(polygons.end(), geom.begin(), geom.end());

            } else if (feature.getType() == mapbox::vector_tile::GeomType::LINESTRING) {
                PointCollectionVec2 geom = feature.getGeometries<PointCollectionVec2>(grid_scale);

                // concat
                // lines.reserve(lines.size() + geom.size());
                // lines.insert(lines.end(), geom.begin(), geom.end());
            }
        }
    }

    return polygons;
}

// polygon describe the outer edge of a closed shape
// -> neighbouring vertices form an edge
// last vertex connects to first vertex
VectorLayerCollection preprocess_triangles(const std::vector<PolygonData> polygons, const std::vector<unsigned int>)
{
    VectorLayerCollection triangle_collection;
    triangle_collection.cell_to_temp = std::vector<std::unordered_set<uint32_t>>(constants::grid_size * constants::grid_size, std::unordered_set<uint32_t>());

    float thickness = 0.0f;

    // TODO static assert to make sure that float are indeed the same size as uint32_t -> so as to not remove some important bits

    size_t data_offset = 0;

    // create the triangles from polygons
    for (size_t i = 0; i < polygons.size(); ++i) {

        // TODO i think triangulize repeats points
        // -> we may be able to reduce the data array, if we increase the index array (if neccessary)

        // polygon to ordered triangles
        std::vector<glm::vec2> triangle_points = nucleus::utils::rasterizer::triangulize(polygons[i].vertices, polygons[i].edges, true);

        // add ordered triangles to collection
        for (size_t j = 0; j < triangle_points.size() / 3; ++j) {
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 0].x));
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 0].y));
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 1].x));
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 1].y));
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 2].x));
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 2].y));
            // triangle_collection.data.push_back(style_indices[i]);
            triangle_collection.data.push_back(1); // hardcoded for now
        }

        const auto cell_writer = [&triangle_collection, data_offset](glm::vec2 pos, int data_index) {
            // if in grid_size bounds than add data_index to vector
            if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(constants::grid_size), pos)))
                triangle_collection.cell_to_temp[int(pos.x) + constants::grid_size * int(pos.y)].insert(data_index + data_offset);
        };

        nucleus::utils::rasterizer::rasterize_triangle(cell_writer, triangle_points, thickness);

        // add to the data offset for the next polygon
        data_offset += triangle_points.size() / 3;
    }

    // std::cout << "tris: " << data_offset << " "; // DEBUG how many triangles were processed

    return triangle_collection;
}

// TODO
VectorLayerCollection preprocess_lines(const std::vector<PolygonData>, const std::vector<unsigned int>)
// VectorLayerCollection preprocess_lines(const std::vector<PolygonData> lines, const std::vector<unsigned int> style_indices)
{
    VectorLayerCollection line_collection;
    line_collection.cell_to_temp = std::vector<std::unordered_set<uint32_t>>(constants::grid_size * constants::grid_size, std::unordered_set<uint32_t>());

    // float thickness = 5.0f;

    // size_t data_offset = 0;

    // // create the triangles from polygons
    // for (size_t i = 0; i < lines.size(); ++i) {
    //     // create_lines(line_collection, lines[i], style_indices[i]);

    //     // TODO currently line points are duplicated, in the final shader -> ideally we only want to store one large point list for one line.
    //     // the problem here is how to also store meta data like style index efficiently (at start or end of line would be ideal but we would need another offset pointer to this location)
    //     // idea -> maybe store styleindex for lines and triangles separately (-> this would also probably be better for triangles)
    //     line_collection.data.reserve(lines[i].size() * 3);
    //     for (size_t j = 0; j < lines[i].size() - 1; j++) {
    //         glm::vec2 p0 = lines[i][j];
    //         glm::vec2 p1 = lines[i][j + 1];
    //         line_collection.data.push_back(*reinterpret_cast<uint32_t*>(&p0.x));
    //         line_collection.data.push_back(*reinterpret_cast<uint32_t*>(&p0.y));
    //         line_collection.data.push_back(*reinterpret_cast<uint32_t*>(&p1.x));
    //         line_collection.data.push_back(*reinterpret_cast<uint32_t*>(&p1.y));
    //         line_collection.data.push_back(style_indices[i]);
    //     }

    //     const auto cell_writer = [&line_collection, data_offset](glm::vec2 pos, int data_index) {
    //         if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(constants::grid_size), pos)))
    //             line_collection.cell_to_temp[int(pos.x) + constants::grid_size * int(pos.y)].insert(data_index + data_offset);
    //     };

    //     nucleus::utils::rasterizer::rasterize_line(cell_writer, lines[i], thickness);

    //     // add to the data offset for the next polyline
    //     // TODO test if this is the correct offset (we might be off by one here but it should suffice for now)
    //     data_offset += lines[i].size();
    // }

    return line_collection;
}

std::vector<std::shared_ptr<const nucleus::Raster<uint32_t>>> split_to_multi_raster(const std::vector<uint32_t>& data)
{
    std::vector<std::shared_ptr<const nucleus::Raster<uint32_t>>> multi_raster;

    constexpr auto squared_size = constants::data_size * constants::data_size;

    for (size_t i = 0; i < floor(double(data.size()) / double(squared_size)); i++) {
        auto layer_data = std::vector<uint32_t>(data.begin() + i * squared_size, data.begin() + (i + 1) * squared_size);
        const auto r = nucleus::Raster<uint32_t>(constants::data_size, std::move(layer_data));
        multi_raster.push_back(std::make_shared<const nucleus::Raster<uint32_t>>(r));
    }

    // add the last of the triangle data (will most likely not be of data_size)
    auto layer_data = std::vector<uint32_t>(data.begin() + multi_raster.size() * squared_size, data.end());
    layer_data.resize(squared_size, -1u); // resize and fill the last layer to the full size
    const auto r = nucleus::Raster<uint32_t>(constants::data_size, std::move(layer_data));
    multi_raster.push_back(std::make_shared<const nucleus::Raster<uint32_t>>(r));

    return multi_raster;
}

/*
 * Function condenses data and fills the GpuVectorLayerTile.
 * condensing:
 *      go over every cell and gather distinct indices
 *      example input(triangle indice for 6 cells): [1], [1], [1,2], [2,3], [2,3], [2,3]
 *      expected output: [[1], [1,2], [2,3]]
 *      Note: it is theoretically possible to further condense the above example to [[1,2,3]] where we can both point to 1, 12 and 23
 *      nevertheless for more complex entries this might be overkill and take more time to compute than it is worth -> only necessary if we need more buffer space
 * simultaneously we also generate the final grid for the tile that stores the offset and size for lookups into the index_bridge
 */
GpuVectorLayerTile create_gpu_tile(const VectorLayerCollection& triangle_collection, const VectorLayerCollection&)
{
    GpuVectorLayerTile tile;

    std::unordered_map<std::unordered_set<uint32_t>, uint32_t, Hasher> unique_entries;
    std::vector<uint32_t> grid;
    std::vector<uint32_t> index_bridge;

    uint32_t start_offset = 0;

    // TODO performance: early exit if not a single thing was written -> currently grid is all 0

    for (size_t i = 0; i < triangle_collection.cell_to_temp.size(); ++i) {
        // go through every cell
        if (triangle_collection.cell_to_temp[i].size() == 0) {
            grid.push_back(0); // no data -> only add an emtpy cell
        } else {
            if (unique_entries.contains(triangle_collection.cell_to_temp[i])) {
                // we already have such an entry -> get the offset_size from there
                grid.push_back(unique_entries[triangle_collection.cell_to_temp[i]]);
            } else {
                // we have to add a new element
                const auto offset_size = nucleus::utils::bit_coding::u24_u8_to_u32(start_offset, uint8_t(triangle_collection.cell_to_temp[i].size()));
                unique_entries[triangle_collection.cell_to_temp[i]] = offset_size;

                grid.push_back(offset_size);

                // add the indices to the bridge
                index_bridge.insert(index_bridge.end(), triangle_collection.cell_to_temp[i].begin(), triangle_collection.cell_to_temp[i].end());

                start_offset += triangle_collection.cell_to_temp[i].size();
            }
        }
    }
    // std::cout << "\t" << index_bridge.size() << ""; // << std::endl;
    // std::cout << "\t" << triangle_collection.data.size() << std::endl;

    constexpr auto max_layer_amount = 8; // 8 * constants::data_size should be sufficient, if not there should be 8 bit available (TODO look if this is still correct)

    assert(index_bridge.size() <= constants::data_size * constants::data_size * max_layer_amount);
    assert(triangle_collection.data.size() <= constants::data_size * constants::data_size * max_layer_amount);

    tile.grid_triangle = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(constants::grid_size, std::move(grid)));
    tile.grid_to_data = split_to_multi_raster(index_bridge);
    tile.data_triangle = split_to_multi_raster(triangle_collection.data);

    // TODO do the same for lines

    return tile;
}

} // namespace nucleus::vector_layer::details
