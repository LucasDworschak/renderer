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

GpuVectorLayerTile preprocess(tile::Id id, const QByteArray& vector_tile_data, const Style& style)
{
    if (vector_tile_data.isEmpty())
        return {};

    // DEBUG polygons
    // const std::vector<std::vector<glm::vec2>> triangle_points = { { glm::vec2(10.5 / 64.0 * constants::grid_size, 30.5 / 64.0 * constants::grid_size),
    //     glm::vec2(30.5 / 64.0 * constants::grid_size, 10.5 / 64.0 * constants::grid_size),
    //     glm::vec2(50.5 / 64.0 * constants::grid_size, 50.5 / 64.0 * constants::grid_size) } };
    const std::vector<unsigned int> style_indices = { 1 };

    auto tile_data = details::parse_tile(id, vector_tile_data, style);

    // qDebug() << id.coords.x << ", " << id.coords.y << " z: " << id.zoom_level;

    auto processed_triangles = details::preprocess_triangles(tile_data);

    auto tile = create_gpu_tile(processed_triangles, processed_triangles);

    return tile;
}

} // namespace nucleus::vector_layer
namespace nucleus::vector_layer::details {

TempDataHolder parse_tile(tile::Id id, const QByteArray& vector_tile_data, const Style& style)
{
    const auto d = vector_tile_data.toStdString();
    const mapbox::vector_tile::buffer tile(d);

    bool first_layer = true;
    float scale = 1.0;

    TempDataHolder data;

    for (const auto& layer_name : tile.layerNames()) {
        // qDebug() << layer_name << id.zoom_level;

        const mapbox::vector_tile::layer layer = tile.getLayer(layer_name);
        std::size_t feature_count = layer.featureCount();

        if (first_layer) {
            first_layer = false;
            data.extent = constants::tile_extent;

            scale = float(constants::tile_extent) / layer.getExtent();
        }

        // if (layer_name != "NUTZUNG_L15_12")
        // if (layer_name != "GEWAESSER_F_GEWF")
        // if (layer_name != "water")
        //     continue; // DEBUG

        for (std::size_t i = 0; i < feature_count; ++i) {
            const auto feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);
            auto props = feature.getProperties();

            const auto type = (feature.getType() == mapbox::vector_tile::GeomType::LINESTRING) ? "line" : "fill";
            const auto style_index = style.layer_style_index(layer_name, type, id.zoom_level, feature);
            if (style_index == -1u) // no style found -> we do not visualize it
                continue;

            if (feature.getType() == mapbox::vector_tile::GeomType::POLYGON) {
                PointCollectionVec2 geom = feature.getGeometries<PointCollectionVec2>(scale);

                // construct edges from the current polygon and move them to the collection
                PolygonData p;
                for (size_t j = 0; j < geom.size(); ++j) {
                    auto edges = nucleus::utils::rasterizer::generate_neighbour_edges(geom[j].size(), p.vertices.size());
                    p.edges.insert(p.edges.end(), edges.begin(), edges.end());
                    p.vertices.insert(p.vertices.end(), geom[j].begin(), geom[j].end());
                }
                data.polygons.push_back(p);
                data.polygon_styles.push_back(style_index);

                // concat
                // polygons.reserve(polygons.size() + geom.size()); // calling reserve here significantly impacts performance -> see if it is possible to call it outside of the loop
                // polygons.insert(polygons.end(), geom.begin(), geom.end());

            } else if (feature.getType() == mapbox::vector_tile::GeomType::LINESTRING) {
                PointCollectionVec2 geom = feature.getGeometries<PointCollectionVec2>(1.0);

                // concat
                // lines.reserve(lines.size() + geom.size());
                // lines.insert(lines.end(), geom.begin(), geom.end());
            }
        }
    }

    qDebug();

    return data;
}

glm::uvec3 pack_triangle_data(glm::vec2 a, glm::vec2 b, glm::vec2 c, uint32_t style_index)
{
    glm::uvec3 out;
    // coordinates 13 bits
    // style 14 bits

    constexpr uint32_t sign_mask = (1u << 31);
    constexpr uint32_t bitmask_12 = (1u << 12) - 1u;
    constexpr uint32_t bitmask_6 = (1u << 6) - 1u;
    constexpr uint32_t bitmask_2 = (1u << 2) - 1u;

    // the extent from the tile gives us e.g. 4096
    // in reality the coordinates can be a little over and a little under the extent -> something like [-128, 4096+128]
    // solution: move all coordinates by half the extent for storing, and move them back when retrieving
    // this way we are using only necessary the necessary bits
    a -= constants::tile_extent / 2.0;
    b -= constants::tile_extent / 2.0;
    c -= constants::tile_extent / 2.0;

    out.x = uint32_t(a.x) << (32u - 13u);
    out.x = out.x | (uint32_t(a.x) & sign_mask);
    out.x = out.x | ((uint32_t(a.y) & bitmask_12) << (32u - 26u));
    out.x = out.x | ((uint32_t(a.y) & sign_mask) >> (13u));

    out.y = uint32_t(b.x) << (32u - 13u);
    out.y = out.y | (uint32_t(b.x) & sign_mask);
    out.y = out.y | ((uint32_t(b.y) & bitmask_12) << (32u - 26u));
    out.y = out.y | ((uint32_t(b.y) & sign_mask) >> (13u));

    out.z = uint32_t(c.x) << (32u - 13u);
    out.z = out.z | (uint32_t(c.x) & sign_mask);
    out.z = out.z | ((uint32_t(c.y) & bitmask_12) << (32u - 26u));
    out.z = out.z | ((uint32_t(c.y) & sign_mask) >> (13u));

    out.x = out.x | ((style_index >> (14u - 6u)) & bitmask_6);
    out.y = out.y | ((style_index >> (14u - 12u)) & bitmask_6);
    out.z = out.z | ((style_index & bitmask_2) << 4u);

    // last 4 bits of out.z are empty

    return out;
}

std::tuple<glm::vec2, glm::vec2, glm::vec2, uint32_t> unpack_triangle_data(glm::uvec3 packed)
{
    constexpr uint32_t sign_mask_first = (1u << 31);
    constexpr uint32_t sign_mask_middle = (1u << (31u - 13u));
    constexpr uint32_t bitmask_12_first = ((1u << 12) - 1u) << (32u - 13u);
    constexpr uint32_t bitmask_12_middle = ((1u << 12) - 1u) << (32u - 26u);
    constexpr uint32_t bitmask_6 = (1u << 6) - 1u;
    constexpr uint32_t bitmask_2 = (1u << 2) - 1u;

    constexpr int32_t negative_bits = ((1u << 31) - 1u) << 12u;

    glm::vec2 a;
    glm::vec2 b;
    glm::vec2 c;
    uint32_t style_index;

    // a.x = ((packed.x & bitmask_12_first) >> (32u - 13u)) ^ -(int32_t((packed.x & sign_mask_first) >> 31));
    // a.x = int32_t((packed.x & bitmask_12_first) >> (32u - 13u)) | (negative_bits & -int32_t((packed.x & sign_mask_first) >> 31));
    // a.y = (packed.x & sign_mask_first);
    // b.x = (packed.x & sign_mask_first) >> 31u;
    // b.y = -int32_t((packed.x & sign_mask_first) >> 31);
    // c.x = negative_bits;
    // c.y = 0;

    a.x = int32_t((packed.x & bitmask_12_first) >> (32u - 13u)) | (negative_bits & -int32_t((packed.x & sign_mask_first) >> 31));
    a.y = int32_t((packed.x & bitmask_12_middle) >> (32u - 26u)) | (negative_bits & -int32_t((packed.x & sign_mask_middle) >> (31 - 13)));
    b.x = int32_t((packed.y & bitmask_12_first) >> (32u - 13u)) | (negative_bits & -int32_t((packed.y & sign_mask_first) >> 31));
    b.y = int32_t((packed.y & bitmask_12_middle) >> (32u - 26u)) | (negative_bits & -int32_t((packed.y & sign_mask_middle) >> (31 - 13)));
    c.x = int32_t((packed.z & bitmask_12_first) >> (32u - 13u)) | (negative_bits & -int32_t((packed.z & sign_mask_first) >> 31));
    c.y = int32_t((packed.z & bitmask_12_middle) >> (32u - 26u)) | (negative_bits & -int32_t((packed.z & sign_mask_middle) >> (31 - 13)));

    style_index = (packed.x & bitmask_6) << (14u - 6u);
    style_index = style_index | (packed.y & bitmask_6) << (14u - 12u);
    style_index = style_index | ((packed.z & (bitmask_2 << 4u)) >> 4u);

    // move the values back to the correct coordinates
    a += constants::tile_extent / 2.0;
    b += constants::tile_extent / 2.0;
    c += constants::tile_extent / 2.0;

    return { a, b, c, style_index };
}

// polygon describe the outer edge of a closed shape
// -> neighbouring vertices form an edge
// last vertex connects to first vertex
VectorLayerCollection preprocess_triangles(const TempDataHolder& data)
{
    VectorLayerCollection triangle_collection;
    triangle_collection.acceleration_grid = std::vector<std::unordered_set<uint32_t>>(constants::grid_size * constants::grid_size, std::unordered_set<uint32_t>());

    float thickness = 0.0f;

    size_t data_offset = 0;

    // create the triangles from polygons
    for (size_t i = 0; i < data.polygons.size(); ++i) {

        // TODO performance i think triangulize repeats points
        // -> we may be able to reduce the data array, if we increase the index array (if neccessary)

        // polygon to ordered triangles
        std::vector<glm::vec2> triangle_points = nucleus::utils::rasterizer::triangulize(data.polygons[i].vertices, data.polygons[i].edges, true);

        // add ordered triangles to collection
        for (size_t j = 0; j < triangle_points.size() / 3; ++j) {
            auto packed = pack_triangle_data(triangle_points[j * 3 + 0], triangle_points[j * 3 + 1], triangle_points[j * 3 + 2], data.polygon_styles[i]);

            triangle_collection.vertex_buffer.push_back(packed);
        }

        const auto cell_writer = [&triangle_collection, data_offset](glm::vec2 pos, int data_index) {
            // if in grid_size bounds than add data_index to vector
            if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(constants::grid_size), pos)))
                triangle_collection.acceleration_grid[int(pos.x) + constants::grid_size * int(pos.y)].insert(data_index + data_offset);
        };

        nucleus::utils::rasterizer::rasterize_triangle(cell_writer, triangle_points, thickness, float(constants::grid_size) / float(data.extent));

        // add to the data offset for the next polygon
        data_offset += triangle_points.size() / 3;
    }

    // std::cout << "tris: " << data_offset << " "; // DEBUG how many triangles were processed

    return triangle_collection;
}

// TODO
VectorLayerCollection preprocess_lines(const TempDataHolder&)
// VectorLayerCollection preprocess_lines(const std::vector<PolygonData> lines, const std::vector<unsigned int> style_indices)
{
    VectorLayerCollection line_collection;
    line_collection.acceleration_grid = std::vector<std::unordered_set<uint32_t>>(constants::grid_size * constants::grid_size, std::unordered_set<uint32_t>());

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
    //             line_collection.acceleration_grid[int(pos.x) + constants::grid_size * int(pos.y)].insert(data_index + data_offset);
    //     };

    //     nucleus::utils::rasterizer::rasterize_line(cell_writer, lines[i], thickness);

    //     // add to the data offset for the next polyline
    //     // TODO test if this is the correct offset (we might be off by one here but it should suffice for now)
    //     data_offset += lines[i].size();
    // }

    return line_collection;
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
    std::vector<uint32_t> acceleration_grid;
    std::vector<uint32_t> index_buffer;
    std::vector<glm::u32vec3> data_triangle = triangle_collection.vertex_buffer;

    uint32_t start_offset = 0;

    // size_t max = 0;

    // TODO performance: early exit if not a single thing was written -> currently grid is all 0

    for (size_t i = 0; i < triangle_collection.acceleration_grid.size(); ++i) {
        // go through every cell
        if (triangle_collection.acceleration_grid[i].size() == 0) {
            acceleration_grid.push_back(0); // no data -> only add an emtpy cell
        } else {
            if (unique_entries.contains(triangle_collection.acceleration_grid[i])) {
                // we already have such an entry -> get the offset_size from there
                acceleration_grid.push_back(unique_entries[triangle_collection.acceleration_grid[i]]);
            } else {
                // we have to add a new element
                // if (triangle_collection.acceleration_grid[i].size() > max)
                //     max = triangle_collection.acceleration_grid[i].size();

                // TODO enable assert again
                // assert(triangle_collection.acceleration_grid[i].size() < 256); // make sure that we are not removing indices we want to draw
                size_t index_buffer_size = triangle_collection.acceleration_grid[i].size();
                if (index_buffer_size > 255)
                    index_buffer_size = 255; // just cap it to 255 as we currently cant go any higher

                const auto offset_size = nucleus::utils::bit_coding::u24_u8_to_u32(start_offset, uint8_t(index_buffer_size));
                unique_entries[triangle_collection.acceleration_grid[i]] = offset_size;

                acceleration_grid.push_back(offset_size);

                // TODO possible performance
                // R32UI 32 bit / index
                // RG32UI 21bit / index -> 3 indices per pixel
                // offset_size of grid could also save 2 bits for where in the rg32ui it should start to pack them tightly

                // add the indices to the bridge
                index_buffer.insert(index_buffer.end(), triangle_collection.acceleration_grid[i].begin(), triangle_collection.acceleration_grid[i].end());

                start_offset += triangle_collection.acceleration_grid[i].size();
            }
        }
    }


    // find fitting cascades for the index and vertex data
    // both index and vertex buffer are set to the same size, as this makes things easier in the GpuMultiArrayHelper
    // TODO performance: check how much index and vertex buffer sizes differ from each other if they differ much you might want to use different cascade levels
    // -> but this causes another possible bottleneck by having to add more GpuArrayHelpers, so that we can accurately
    //      distinguish between index and vertex buffer in GpuMultiArrayHelper::generate_dictionary()
    uint8_t fitting_cascade_index = uint8_t(-1u);
    for (uint8_t i = 0; i < constants::data_size.size(); i++) {
        const auto data_size_sq = constants::data_size[i] * constants::data_size[i];
        if (index_buffer.size() <= data_size_sq && data_triangle.size() <= data_size_sq) {
            fitting_cascade_index = i;
            break;
        }
    }

    // qDebug() << "index: " << start_offset << fitting_cascade_index;

    // if assert is triggered -> consider adding a value to constants::data_size
    assert(fitting_cascade_index < constants::data_size.size());

    // resize data to buffer size
    index_buffer.resize(constants::data_size[fitting_cascade_index] * constants::data_size[fitting_cascade_index], -1u);
    data_triangle.resize(constants::data_size[fitting_cascade_index] * constants::data_size[fitting_cascade_index], glm::u32vec3(-1u));

    tile.triangle_acceleration_grid = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(constants::grid_size, std::move(acceleration_grid)));
    tile.triangle_index_buffer = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(constants::data_size[fitting_cascade_index], std::move(index_buffer)));
    tile.triangle_vertex_buffer = std::make_shared<const nucleus::Raster<glm::u32vec3>>(nucleus::Raster<glm::u32vec3>(constants::data_size[fitting_cascade_index], std::move(data_triangle)));

    tile.triangle_buffer_info = fitting_cascade_index;

    // TODO do the same for lines

    return tile;
}

} // namespace nucleus::vector_layer::details
