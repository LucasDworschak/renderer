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

#include "nucleus/Raster.h"
#include <nucleus/utils/bit_coding.h>

#include "nucleus/utils/rasterizer.h"

#include <mapbox/vector_tile.hpp>

#include "constants.h"

namespace nucleus::vector_layer {

GpuVectorLayerTile create_default_gpu_tile()
{
    GpuVectorLayerTile default_tile;
    default_tile.buffer_info = 0;
    default_tile.acceleration_grid = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(glm::uvec2(constants::grid_size), 0));
    default_tile.index_buffer = std::make_shared<const nucleus::Raster<uint32_t>>(
        nucleus::Raster<uint32_t>(glm::uvec2(constants::data_size[0] * constants::index_buffer_size_multiplier), -1u));
    default_tile.vertex_buffer
        = std::make_shared<const nucleus::Raster<glm::u32vec3>>(nucleus::Raster<glm::u32vec3>(glm::uvec2(constants::data_size[0]), glm::u32vec3(-1u)));

    return default_tile;
}

GpuVectorLayerTile preprocess(tile::Id id, const QByteArray& vector_tile_data, const Style& style)
{
    if (vector_tile_data.isEmpty())
        return {};
    // qDebug() << id.coords.x << ", " << id.coords.y << " z: " << id.zoom_level;

    // DEBUG polygons
    // const std::vector<std::vector<glm::vec2>> triangle_points = { { glm::vec2(10.5 / 64.0 * constants::grid_size, 30.5 / 64.0 * constants::grid_size),
    //     glm::vec2(30.5 / 64.0 * constants::grid_size, 10.5 / 64.0 * constants::grid_size),
    //     glm::vec2(50.5 / 64.0 * constants::grid_size, 50.5 / 64.0 * constants::grid_size) } };
    const std::vector<unsigned int> style_indices = { 1 };

    auto tile_data = details::parse_tile(id, vector_tile_data, style);

    const auto style_buffer = style.styles()->buffer();
    auto layer_collection = details::preprocess_geometry(tile_data, style_buffer);

    auto tile = create_gpu_tile(layer_collection);
    tile.id = id;

    return tile;
}

} // namespace nucleus::vector_layer
namespace nucleus::vector_layer::details {

std::vector<std::pair<uint32_t, uint32_t>> simplify_styles(std::vector<std::pair<uint32_t, uint32_t>> style_and_layer_indices, const std::vector<glm::u32vec4> style_buffer)
{
    // we get multiple styles that may have full opacity and the same width
    // creating render calls for both does not make sense -> we only want to draw the top layer
    // this function simplifys all the styles so that only the styles which actually have a change to be rendered will remain

    // order the styles so that we look at layer in descending order
    std::sort(style_and_layer_indices.begin(), style_and_layer_indices.end(), [](std::pair<uint32_t, uint32_t> a, std::pair<uint32_t, uint32_t> b) { return a.second > b.second; });

    std::vector<std::pair<uint32_t, uint32_t>> out_styles;
    int accummulative_opacity = 0;
    float width = 0.0;

    for (const auto& indices : style_and_layer_indices) {
        const auto style_data = style_buffer[indices.first >> 1];
        const float current_width = float(style_data.z) / float(constants::style_precision);
        const int current_opacity = style_data.x & 255;

        if (width < current_width) {
            // reset opacity
            accummulative_opacity = 0;
            width = current_width;
        }

        if (accummulative_opacity < 255) {
            accummulative_opacity += current_opacity;
            out_styles.push_back(indices);
        }
    }

    return out_styles;
}

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
            auto style_and_layer_indices = style.indices(layer_name, type, id.zoom_level, feature);
            style_and_layer_indices = simplify_styles(style_and_layer_indices, style_buffer);

            if (style_and_layer_indices.size() == 0) // no styles found -> we do not visualize it
                continue;

            const auto is_polygon = feature.getType() == mapbox::vector_tile::GeomType::POLYGON;
            PointCollectionVec2 geom = feature.getGeometries<PointCollectionVec2>(scale);
            data.emplace_back(std::vector<std::vector<glm::vec2>>(geom.begin(), geom.end()), extent, style_and_layer_indices, is_polygon);
        }
    }

    return data;
}

uint32_t pack_index_buffer(uint32_t geometry_index, uint32_t style_index, bool polygon) { return (geometry_index << (constants::style_bits + 1)) | (style_index << 1) | (polygon ? 1u : 0u); }

glm::uvec3 pack_line_data(glm::vec2 a, glm::vec2 b)
{
    // TODO -> possibly need to store additional values here in the future -> refactor this and pack_triangle_data to pack it efficiently
    return pack_triangle_data(a, b, { 0, 0 });
}

glm::uvec3 pack_triangle_data(glm::vec2 a, glm::vec2 b, glm::vec2 c)
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

    out.x = uint32_t(a.x) << coordinate_shift1;
    out.x = out.x | ((uint32_t(a.y) & coordinate_bitmask) << coordinate_shift2);

    out.y = uint32_t(b.x) << coordinate_shift1;
    out.y = out.y | ((uint32_t(b.y) & coordinate_bitmask) << coordinate_shift2);

    out.z = uint32_t(c.x) << coordinate_shift1;
    out.z = out.z | ((uint32_t(c.y) & coordinate_bitmask) << coordinate_shift2);

    return out;
}

std::tuple<glm::vec2, glm::vec2, glm::vec2> unpack_triangle_data(glm::uvec3 packed)
{
    glm::vec2 a;
    glm::vec2 b;
    glm::vec2 c;

    a.x = int32_t((packed.x & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    a.y = int32_t((packed.x & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);
    b.x = int32_t((packed.y & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    b.y = int32_t((packed.y & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);
    c.x = int32_t((packed.z & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    c.y = int32_t((packed.z & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);

    // move the values back to the correct coordinates
    a -= constants::tile_extent;
    b -= constants::tile_extent;
    c -= constants::tile_extent;

    return { a, b, c };
}

// polygon describe the outer edge of a closed shape
// -> neighbouring vertices form an edge
// last vertex connects to first vertex
VectorLayerCollection preprocess_geometry(const std::vector<GeometryData>& data, const std::vector<glm::u32vec4> style_buffer)
{
    VectorLayerCollection layer_collection;
    // layer_collection.acceleration_grid = std::vector<std::set<uint32_t>>(constants::grid_size * constants::grid_size, std::set<uint32_t>());
    layer_collection.acceleration_grid = std::vector<std::map<uint32_t, uint32_t>>(constants::grid_size * constants::grid_size, std::map<uint32_t, uint32_t>());

    uint32_t data_offset = 0;

    // create the triangles from polygons
    for (size_t i = 0; i < data.size(); ++i) {

        if (data[i].is_polygon) {

            // TODO performance i think triangulize repeats points
            // -> we may be able to reduce the data array, if we increase the index array (if neccessary)

            // polygon to ordered triangles
            std::vector<glm::vec2> triangle_points = nucleus::utils::rasterizer::triangulize(data[i].vertices, true);

            // add ordered triangles to collection
            for (size_t j = 0; j < triangle_points.size() / 3; ++j) {
                auto packed = pack_triangle_data(triangle_points[j * 3 + 0], triangle_points[j * 3 + 1], triangle_points[j * 3 + 2]);

                layer_collection.vertex_buffer.push_back(packed);
            }

            for (const auto& style_layer : data[i].style_and_layer_indices) {
                const auto cell_writer = [&layer_collection, data_offset, style_layer](glm::vec2 pos, int data_index) {
                    // if in grid_size bounds and not already present -> than add index to vector
                    if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(constants::grid_size), pos))) {
                        // index:
                        // 18 bits for the index in geometry buffer
                        // 12 bits for style (see constants if still the same)
                        // 1 bit to encode if this is a polygon or not
                        const auto index = pack_index_buffer((data_index + data_offset), style_layer.first, true);
                        const auto layer = (style_layer.second << (32 - constants::style_bits)) | (data_index + data_offset);
                        layer_collection.acceleration_grid[int(pos.x) + constants::grid_size * int(pos.y)][layer] = index;
                    }
                };
                const auto scale = float(constants::grid_size) / float(data[i].extent);
                nucleus::utils::rasterizer::rasterize_triangle(cell_writer, triangle_points, 0.0f, scale);
            }
            // add to the data offset for the next polygon
            data_offset += triangle_points.size() / 3u;

        } else {
            // polylines
            for (size_t j = 0; j < data[i].vertices.size(); ++j) {
                for (size_t k = 0; k < data[i].vertices[j].size() - 1; ++k) {
                    auto packed = pack_line_data(data[i].vertices[j][k], data[i].vertices[j][k + 1]);
                    layer_collection.vertex_buffer.push_back(packed);
                }
            }

            for (const auto& style_layer : data[i].style_and_layer_indices) {
                const auto line_width = float(style_buffer[style_layer.first >> 1].z) / float(constants::style_precision);

                const auto cell_writer = [&layer_collection, data_offset, style_layer](glm::vec2 pos, int data_index) {
                    // if in grid_size bounds and not already present -> than add index to vector
                    if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(constants::grid_size), pos))) {
                        // index:
                        // 18 bits for the index in geometry buffer
                        // 13 bits for style (see constants if still the same)
                        // 1 bit to encode if this is a polygon or not
                        const auto index = pack_index_buffer((data_index + data_offset), style_layer.first, false);
                        const auto layer = (style_layer.second << (32 - constants::style_bits)) | (data_index + data_offset);
                        layer_collection.acceleration_grid[int(pos.x) + constants::grid_size * int(pos.y)][layer] = index;
                    }
                };

                const auto scale = float(constants::grid_size) / float(data[i].extent);

                for (size_t j = 0; j < data[i].vertices.size(); ++j) {
                    // TODO: according to Task #151 -> we doubled the line width that goes into the acceleration structure because we are looking at tiles that are bigger
                    // Nevertheless, we artificially worsened the performance by introducing more cells where a line could be (although it is only there on specific zoom levels)
                    // This performance issue will be solved with Task #198 (mipmaps)
                    nucleus::utils::rasterizer::rasterize_line(cell_writer, data[i].vertices[j], line_width * scale * 2.0, scale);
                }
            }
            for (size_t j = 0; j < data[i].vertices.size(); ++j) {
                data_offset += data[i].vertices[j].size() - 1u;
            }
        }
    }

    // make sure that data_offset is not truncated
    assert(data_offset < (1u << (32 - constants::style_bits - 1)));

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

    std::unordered_map<std::vector<uint32_t>, uint32_t, Hasher> unique_entries;
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
            auto cell_indices = std::vector<uint32_t>();
            std::transform(
                layer_collection.acceleration_grid[i].crbegin(), layer_collection.acceleration_grid[i].crend(), std::back_inserter(cell_indices), [](const std::pair<uint32_t, uint32_t> pair) {
                    return pair.second;
                });

            if (unique_entries.contains(cell_indices)) {
                // we already have such an entry -> get the offset_size from there
                acceleration_grid.push_back(unique_entries[cell_indices]);
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
                unique_entries[cell_indices] = offset_size;

                acceleration_grid.push_back(offset_size);

                // TODO possible performance
                // change R32UI to RGBA32UI to save 4 indices at the same time
                // -> we most likely want to test a few indices one after another and this might be better
                // also possible to use a 2 bits in the offset_size to indicate the start (to pack the indices tightly even if the previous indices belong to previous cell)

                // add the indices to the bridge
                index_buffer.insert(index_buffer.end(), cell_indices.begin(), cell_indices.end());

                start_offset += cell_indices.size();
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
        if (index_buffer.size() <= data_size_sq * constants::index_buffer_size_multiplier * constants::index_buffer_size_multiplier && data_triangle.size() <= data_size_sq) {
            fitting_cascade_index = i;
            break;
        }
    }

    // qDebug() << "cascade_data:" << index_buffer.size() << data_triangle.size();

    // qDebug() << "index: " << start_offset << fitting_cascade_index;

    // qDebug() << "indices: " << index_buffer.size() << data_triangle.size();

    // if assert is triggered -> consider adding a value to constants::data_size
    if (fitting_cascade_index >= constants::data_size.size())
        qDebug() << index_buffer.size() << data_triangle.size() << "data does not fit";
    assert(fitting_cascade_index < constants::data_size.size());

    // resize data to buffer size
    index_buffer.resize(
        constants::data_size[fitting_cascade_index] * constants::data_size[fitting_cascade_index] * constants::index_buffer_size_multiplier * constants::index_buffer_size_multiplier, -1u);
    data_triangle.resize(constants::data_size[fitting_cascade_index] * constants::data_size[fitting_cascade_index], glm::u32vec3(-1u));

    tile.acceleration_grid = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(constants::grid_size, std::move(acceleration_grid)));
    tile.index_buffer
        = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(constants::data_size[fitting_cascade_index] * constants::index_buffer_size_multiplier, std::move(index_buffer)));
    tile.vertex_buffer = std::make_shared<const nucleus::Raster<glm::u32vec3>>(nucleus::Raster<glm::u32vec3>(constants::data_size[fitting_cascade_index], std::move(data_triangle)));

    tile.buffer_info = fitting_cascade_index;

    // TODO do the same for lines

    return tile;
}

} // namespace nucleus::vector_layer::details
