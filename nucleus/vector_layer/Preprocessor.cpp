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

// TODOs
// - line visualization

GpuVectorLayerTile preprocess(tile::Id id, const QByteArray& vector_tile_data, const Style& style)
{
    // qDebug() << id.coords.x << ", " << id.coords.y << " z: " << id.zoom_level;
    if (vector_tile_data.isEmpty())
        return {};
    // qDebug() << "b";

    // DEBUG polygons
    // const std::vector<std::vector<glm::vec2>> triangle_points = { { glm::vec2(10.5 / 64.0 * constants::grid_size, 30.5 / 64.0 * constants::grid_size),
    //     glm::vec2(30.5 / 64.0 * constants::grid_size, 10.5 / 64.0 * constants::grid_size),
    //     glm::vec2(50.5 / 64.0 * constants::grid_size, 50.5 / 64.0 * constants::grid_size) } };
    const std::vector<unsigned int> style_indices = { 1 };

    auto tile_data = details::parse_tile(id, vector_tile_data, style);

    auto layer_collection = details::preprocess_geometry(tile_data);

    auto tile = create_gpu_tile(layer_collection);

    return tile;
}

} // namespace nucleus::vector_layer
namespace nucleus::vector_layer::details {

std::vector<GeometryData> parse_tile(tile::Id id, const QByteArray& vector_tile_data, const Style& style)
{
    const auto d = vector_tile_data.toStdString();
    const mapbox::vector_tile::buffer tile(d);

    bool first_layer = true;
    float scale = 1.0;

    std::vector<GeometryData> data;

    uint32_t extent;

    const auto style_buffer = style.styles()->buffer();

    for (const auto& layer_name : tile.layerNames()) {
        // qDebug() << layer_name << id.zoom_level;

        const mapbox::vector_tile::layer layer = tile.getLayer(layer_name);
        std::size_t feature_count = layer.featureCount();

        if (first_layer) {
            first_layer = false;
            extent = constants::tile_extent;

            scale = float(constants::tile_extent) / layer.getExtent();
        }

        for (std::size_t i = 0; i < feature_count; ++i) {
            const auto feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);

            const auto type = (feature.getType() == mapbox::vector_tile::GeomType::LINESTRING) ? "line" : "fill";
            auto styles = style.indices(layer_name, type, id.zoom_level, feature);

            if (styles.size() == 0) // no styles found -> we do not visualize it
                continue;

            if (feature.getType() == mapbox::vector_tile::GeomType::POLYGON) {
                PointCollectionVec2 geom = feature.getGeometries<PointCollectionVec2>(scale);

                // construct edges from the current polygon and move them to the collection
                std::vector<glm::vec2> vertices;
                std::vector<glm::ivec2> all_edges;
                for (size_t j = 0; j < geom.size(); ++j) {
                    auto edges = nucleus::utils::rasterizer::generate_neighbour_edges(geom[j].size(), vertices.size());
                    all_edges.insert(all_edges.end(), edges.begin(), edges.end());
                    vertices.insert(vertices.end(), geom[j].begin(), geom[j].end());
                }

                // qDebug() << layer_name << styles.size();

                for (const auto& style : styles)
                    data.emplace_back(vertices, extent, style.first, style.second, true, all_edges, 0);
                // data.emplace_back(std::vector<glm::vec2>(vertices), extent, styles[0].first, styles[0].second, true, std::vector<glm::ivec2>(edges), 0);
                // data.push_back({ vertices, extent, styles[0].first, styles[0].second, true, all_edges, 0 });

            } else if (feature.getType() == mapbox::vector_tile::GeomType::LINESTRING) {
                PointCollectionVec2 geom = feature.getGeometries<PointCollectionVec2>(scale);

                for (size_t j = 0; j < geom.size(); ++j) {
                    // lines have to be split according to geom -> otherwise Bug #171
                    std::vector<glm::vec2> vertices;
                    vertices.insert(vertices.end(), geom[j].begin(), geom[j].end());

                    constexpr std::vector<glm::ivec2> no_edges; // necessary for emplace_back

                    // TODO performance -> instead of duplicating the vertice data we only want to duplicate the index data

                    for (const auto& style : styles) {
                        const auto line_width = float(style_buffer[style.first].z) / float(constants::style_precision);
                        data.emplace_back(vertices, extent, style.first, style.second, false, no_edges, line_width);
                    }
                }
            }
        }
    }

    // sort the geometry, so that the layer is in a descending order
    // layer in descending order is needed to have a top down drawing order
    std::sort(data.begin(), data.end(), [](GeometryData a, GeometryData b) { return a.layer > b.layer; });

    return data;
}

glm::uvec3 pack_line_data(glm::vec2 a, glm::vec2 b, uint32_t style_index)
{
    // TODO -> possibly need to store additional values here in the future -> refactor this and pack_triangle_data to pack it efficiently
    return pack_triangle_data(a, b, {}, style_index);
}

glm::uvec3 pack_triangle_data(glm::vec2 a, glm::vec2 b, glm::vec2 c, uint32_t style_index)
{
    glm::uvec3 out;

    // the extent from the tile gives us e.g. 4096
    // in reality the coordinates can be a little over and a little under the extent -> something like [-128, 4096+128]
    // solution: coordinate normalization (move from [-tile_extent - tile_extent] to [0 - 2*tile_extent])
    a += constants::tile_extent;
    b += constants::tile_extent;
    c += constants::tile_extent;

    // make sure that we do not remove bits from the coordinates
    assert((uint32_t(a.x) & ((1u << coordinate_bits) - 1u)) == uint32_t(a.x));
    assert((uint32_t(a.y) & ((1u << coordinate_bits) - 1u)) == uint32_t(a.y));
    assert((uint32_t(b.x) & ((1u << coordinate_bits) - 1u)) == uint32_t(b.x));
    assert((uint32_t(b.y) & ((1u << coordinate_bits) - 1u)) == uint32_t(b.y));
    assert((uint32_t(c.x) & ((1u << coordinate_bits) - 1u)) == uint32_t(c.x));
    assert((uint32_t(c.y) & ((1u << coordinate_bits) - 1u)) == uint32_t(c.y));

    // make sure that we do not remove bits from style_index
    assert((style_index & ((1u << all_style_bits) - 1u)) == style_index);

    out.x = uint32_t(a.x) << coordinate_shift1;
    out.x = out.x | ((uint32_t(a.y) & coordinate_bitmask) << coordinate_shift2);

    out.y = uint32_t(b.x) << coordinate_shift1;
    out.y = out.y | ((uint32_t(b.y) & coordinate_bitmask) << coordinate_shift2);

    out.z = uint32_t(c.x) << coordinate_shift1;
    out.z = out.z | ((uint32_t(c.y) & coordinate_bitmask) << coordinate_shift2);

    out.x = out.x | ((style_index >> style_shift1) & style_bitmask);
    out.y = out.y | ((style_index >> style_shift2) & style_bitmask);
    out.z = out.z | ((style_index & style_bitmask));

    return out;
}

std::tuple<glm::vec2, glm::vec2, glm::vec2, uint32_t> unpack_triangle_data(glm::uvec3 packed)
{
    glm::vec2 a;
    glm::vec2 b;
    glm::vec2 c;
    uint32_t style_index;

    a.x = int32_t((packed.x & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    a.y = int32_t((packed.x & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);
    b.x = int32_t((packed.y & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    b.y = int32_t((packed.y & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);
    c.x = int32_t((packed.z & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    c.y = int32_t((packed.z & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);

    style_index = (packed.x & style_bitmask) << style_shift1;
    style_index = style_index | (packed.y & style_bitmask) << style_shift2;
    style_index = style_index | (packed.z & style_bitmask);

    // move the values back to the correct coordinates
    a -= constants::tile_extent;
    b -= constants::tile_extent;
    c -= constants::tile_extent;

    return { a, b, c, style_index };
}

// polygon describe the outer edge of a closed shape
// -> neighbouring vertices form an edge
// last vertex connects to first vertex
VectorLayerCollection preprocess_geometry(const std::vector<GeometryData>& data)
{
    VectorLayerCollection layer_collection;
    layer_collection.acceleration_grid = std::vector<std::set<uint32_t>>(constants::grid_size * constants::grid_size, std::set<uint32_t>());

    uint32_t data_offset = 0;

    // create the triangles from polygons
    for (size_t i = 0; i < data.size(); ++i) {

        if (data[i].is_polygon) {

            // TODO performance i think triangulize repeats points
            // -> we may be able to reduce the data array, if we increase the index array (if neccessary)

            // polygon to ordered triangles
            std::vector<glm::vec2> triangle_points = nucleus::utils::rasterizer::triangulize(data[i].vertices, data[i].edges, true);

            // add ordered triangles to collection
            for (size_t j = 0; j < triangle_points.size() / 3; ++j) {
                auto packed = pack_triangle_data(triangle_points[j * 3 + 0], triangle_points[j * 3 + 1], triangle_points[j * 3 + 2], data[i].style);

                layer_collection.vertex_buffer.push_back(packed);
            }

            const auto cell_writer = [&layer_collection, data_offset](glm::vec2 pos, int data_index) {
                // if in grid_size bounds and not already present -> than add index to vector
                if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(constants::grid_size), pos))) {
                    // last bit of index indicates that this is a polygon
                    // !! IMPORTANT !! we have to use the lowest bit since set orders the input depending on key -> lines and polygons have to stay intermixed
                    const auto index = ((data_index + data_offset) << 1) | 1u;
                    layer_collection.acceleration_grid[int(pos.x) + constants::grid_size * int(pos.y)].insert(index);
                }
            };

            nucleus::utils::rasterizer::rasterize_triangle(cell_writer, triangle_points, 0.0f, float(constants::grid_size) / float(data[i].extent));

            // add to the data offset for the next polygon
            data_offset += triangle_points.size() / 3u;
        } else {
            // polylines
            for (size_t j = 0; j < data[i].vertices.size() - 1; ++j) {
                auto packed = pack_line_data(data[i].vertices[j], data[i].vertices[j + 1], data[i].style);
                layer_collection.vertex_buffer.push_back(packed);
            }

            const auto cell_writer = [&layer_collection, data_offset](glm::vec2 pos, int data_index) {
                // if in grid_size bounds and not already present -> than add index to vector
                if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(constants::grid_size), pos))) {
                    // last bit of index indicates that this is a line
                    // !! IMPORTANT !! we have to use the lowest bit since set orders the input depending on key -> lines and polygons have to stay intermixed
                    const auto index = ((data_index + data_offset) << 1); // | 0u;
                    layer_collection.acceleration_grid[int(pos.x) + constants::grid_size * int(pos.y)].insert(index);
                }
            };

            const auto scale = float(constants::grid_size) / float(data[i].extent);

            nucleus::utils::rasterizer::rasterize_line(cell_writer, data[i].vertices, data[i].line_width * scale, scale);

            data_offset += data[i].vertices.size() - 1u;
        }
    }

    // make sure that data_offset does not use more than 31 bits
    assert(data_offset < (1u << 31));

    // qDebug() << "num indices: " << data_offset; // DEBUG how many indices are used

    return layer_collection;
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
GpuVectorLayerTile create_gpu_tile(const VectorLayerCollection& layer_collection)
{
    GpuVectorLayerTile tile;

    std::unordered_map<std::set<uint32_t>, uint32_t, Hasher> unique_entries;
    std::vector<uint32_t> acceleration_grid;
    std::vector<uint32_t> index_buffer;
    std::vector<glm::u32vec3> data_triangle = layer_collection.vertex_buffer;

    uint32_t start_offset = 0;

    // size_t max = 0;

    // TODO performance: early exit if not a single thing was written -> currently grid is all 0

    for (size_t i = 0; i < layer_collection.acceleration_grid.size(); ++i) {
        // go through every cell
        if (layer_collection.acceleration_grid[i].size() == 0) {
            acceleration_grid.push_back(0); // no data -> only add an emtpy cell
        } else {
            if (unique_entries.contains(layer_collection.acceleration_grid[i])) {
                // we already have such an entry -> get the offset_size from there
                acceleration_grid.push_back(unique_entries[layer_collection.acceleration_grid[i]]);
            } else {
                // we have to add a new element
                // if (layer_collection.acceleration_grid[i].size() > max)
                //     max = layer_collection.acceleration_grid[i].size();

                // TODO enable assert again
                // assert(layer_collection.acceleration_grid[i].size() < 256); // make sure that we are not removing indices we want to draw
                size_t index_buffer_size = layer_collection.acceleration_grid[i].size();
                if (index_buffer_size > 255)
                    index_buffer_size = 255; // just cap it to 255 as we currently cant go any higher

                const auto offset_size = nucleus::utils::bit_coding::u24_u8_to_u32(start_offset, uint8_t(index_buffer_size));
                unique_entries[layer_collection.acceleration_grid[i]] = offset_size;

                acceleration_grid.push_back(offset_size);

                // TODO possible performance
                // R32UI 32 bit / index
                // RG32UI 21bit / index -> 3 indices per pixel
                // offset_size of grid could also save 2 bits for where in the rg32ui it should start to pack them tightly

                // add the indices to the bridge
                index_buffer.insert(index_buffer.end(), layer_collection.acceleration_grid[i].begin(), layer_collection.acceleration_grid[i].end());

                start_offset += layer_collection.acceleration_grid[i].size();
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
    if (fitting_cascade_index >= constants::data_size.size())
        qDebug() << index_buffer.size() << data_triangle.size() << "data does not fit";
    assert(fitting_cascade_index < constants::data_size.size());

    // resize data to buffer size
    index_buffer.resize(constants::data_size[fitting_cascade_index] * constants::data_size[fitting_cascade_index], -1u);
    data_triangle.resize(constants::data_size[fitting_cascade_index] * constants::data_size[fitting_cascade_index], glm::u32vec3(-1u));

    tile.acceleration_grid = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(constants::grid_size, std::move(acceleration_grid)));
    tile.index_buffer = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(constants::data_size[fitting_cascade_index], std::move(index_buffer)));
    tile.vertex_buffer = std::make_shared<const nucleus::Raster<glm::u32vec3>>(nucleus::Raster<glm::u32vec3>(constants::data_size[fitting_cascade_index], std::move(data_triangle)));

    tile.buffer_info = fitting_cascade_index;

    // TODO do the same for lines

    return tile;
}

} // namespace nucleus::vector_layer::details
