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

#include <unordered_set>

#include "constants.h"
#include "nucleus/Raster.h"
#include "nucleus/utils/rasterizer.h"
#include <nucleus/utils/bit_coding.h>

#include <earcut.hpp>
#include <glm/gtx/closest_point.hpp>
#include <mapbox/vector_tile.hpp>

// TODO alternative to class approach is going back to function approach and passing through a struct/class that holds all temporary heap data
// this object is passed through all functions and cleared afterwards (when it has time)
// furthermore multiple of those data holder structs can be constructed and put into a pool so that multithreading could be utilized

// allow vec2 points for earcut input
namespace mapbox {
namespace util {

    template <>
    struct nth<0, glm::vec2> {
        inline static auto get(const glm::vec2& t) { return t.x; };
    };
    template <>
    struct nth<1, glm::vec2> {
        inline static auto get(const glm::vec2& t) { return t.y; };
    };

    template <>
    struct nth<0, nucleus::vector_layer::ClipperPoint> {
        inline static auto get(const nucleus::vector_layer::ClipperPoint& t) { return t.x; };
    };
    template <>
    struct nth<1, nucleus::vector_layer::ClipperPoint> {
        inline static auto get(const nucleus::vector_layer::ClipperPoint& t) { return t.y; };
    };

} // namespace util
} // namespace mapbox

namespace nucleus::vector_layer {

Preprocessor::Preprocessor(Style&& style)
    : m_style(std::move(style))
    , m_style_buffer(m_style.styles()->buffer())
{
    generate_preprocess_grid();
}

const Style& Preprocessor::style() { return m_style; }
size_t Preprocessor::processed_amount() { return m_processed_amount; }

GpuVectorLayerTile Preprocessor::create_default_gpu_tile()
{
    GpuVectorLayerTile default_tile;
    default_tile.buffer_info = 0;
    default_tile.acceleration_grid = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(glm::uvec2(constants::grid_size), 0));
    default_tile.geometry_buffer
        = std::make_shared<const nucleus::Raster<glm::u32vec2>>(nucleus::Raster<glm::u32vec2>(glm::uvec2(constants::data_size[0]), glm::u32vec2(-1u)));

    return default_tile;
}

GpuVectorLayerTile Preprocessor::preprocess(tile::Id id, const QByteArray& vector_tile_data)
{
    if (vector_tile_data.isEmpty())
        return {};
    // qDebug() << id.coords.x << ", " << id.coords.y << " z: " << id.zoom_level;

    // DEBUG polygons
    // const std::vector<std::vector<glm::vec2>> triangle_points = { { glm::vec2(10.5 / 64.0 * constants::grid_size, 30.5 / 64.0 * constants::grid_size),
    //     glm::vec2(30.5 / 64.0 * constants::grid_size, 10.5 / 64.0 * constants::grid_size),
    //     glm::vec2(50.5 / 64.0 * constants::grid_size, 50.5 / 64.0 * constants::grid_size) } };
    // const std::vector<unsigned int> style_indices = { 1 };

    auto tile_data = parse_tile(id, vector_tile_data);

    preprocess_geometry(tile_data);

    // return if not a single thing was written
    if (m_processed_amount == 0u)
        return {};

    auto tile = create_gpu_tile();
    tile.id = id;

    return tile;
}

std::vector<std::pair<uint32_t, uint32_t>> Preprocessor::simplify_styles(
    std::vector<std::pair<uint32_t, uint32_t>>* style_and_layer_indices, const std::vector<glm::u32vec4>& style_buffer)
{
    // we get multiple styles that may have full opacity and the same width
    // creating render calls for both does not make sense -> we only want to draw the top layer
    // this function simplifys all the styles so that only the styles which actually have a change to be rendered will remain

    // TODO this sort should happen at creation of the vector not here
    // order the styles so that we look at layer in descending order
    std::sort(style_and_layer_indices->begin(), style_and_layer_indices->end(), [](std::pair<uint32_t, uint32_t> a, std::pair<uint32_t, uint32_t> b) {
        return a.second > b.second;
    });

    std::vector<std::pair<uint32_t, uint32_t>> out_styles;
    int accummulative_opacity = 0;
    float width = 0.0;

    for (const auto& indices : *style_and_layer_indices) {
        const auto style_data = style_buffer[indices.first >> 1];
        const float current_width = float(style_data.z) / float(constants::style_precision);
        const int current_opacity = style_data.x & 255;

        if (current_opacity == 0)
            continue; // we dont care about 0 opacity geometry

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

VectorLayers Preprocessor::parse_tile(tile::Id id, const QByteArray& vector_tile_data)
{
    const auto d = vector_tile_data.toStdString();
    const mapbox::vector_tile::buffer tile(d);

    // bool first_layer = true;
    // uint32_t extent;
    constexpr float scale = 1.0;

    constexpr auto cell_scale = float(constants::grid_size) / float(constants::tile_extent);

    VectorLayers data;

    std::array<int, constants::max_style_expression_keys> temp_values;

    for (const auto& layer_name : tile.layerNames()) {
        // qDebug() << layer_name << id.zoom_level;

        const mapbox::vector_tile::layer layer = tile.getLayer(layer_name);
        std::size_t feature_count = layer.featureCount();

        assert(layer.getExtent() == constants::tile_extent);

        for (std::size_t i = 0; i < feature_count; ++i) {
            const auto feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);

            const auto type = (feature.getType() == mapbox::vector_tile::GeomType::LINESTRING) ? 0 : 1;
            // qDebug() << layer_name;
            auto style_and_layer_indices = m_style.indices(layer_name, type, id.zoom_level, feature, &temp_values);
            style_and_layer_indices = simplify_styles(&style_and_layer_indices, m_style_buffer);

            if (style_and_layer_indices.size() == 0) // no styles found -> we do not visualize it
                continue;

            const auto is_polygon = feature.getType() == mapbox::vector_tile::GeomType::POLYGON;
            PointCollectionVec2 geom = feature.getGeometries<PointCollectionVec2>(scale);

            radix::geometry::Aabb2i aabb({ constants::grid_size * 2, constants::grid_size * 2 }, { -constants::grid_size * 2, -constants::grid_size * 2 });
            // calculate bounds
            std::vector<ClipperRect> bounds;
            if (is_polygon) {
                // only needed for polygons
                bounds.reserve(geom.size());
                for (size_t j = 0; j < geom.size(); j++) {
                    const auto bound = Clipper2Lib::GetBounds(geom[j]);

                    aabb.expand_by({ std::floor(float(bound.left) * cell_scale), std::floor(float(bound.top) * cell_scale) });
                    aabb.expand_by({ std::ceil(float(bound.right) * cell_scale), std::ceil(float(bound.bottom) * cell_scale) });
                    bounds.push_back(bound);
                }

                // TODO if aabb does not intersect with grid at all -> continue;
            }

            for (const auto& style_layer : style_and_layer_indices) {
                data[style_layer.second].emplace_back(ClipperPaths(geom.begin(), geom.end()), bounds, aabb, style_layer, is_polygon);
            }
        }
    }

    return data;
}

glm::u32vec2 Preprocessor::pack_line_data(glm::i64vec2 a, glm::i64vec2 b, uint16_t style_index) { return pack_triangle_data({ a, b, b, style_index, false }); }

glm::u32vec2 Preprocessor::pack_triangle_data(VectorLayerData data)
{
    glm::u32vec2 packed_data;

    data.a += geometry_offset;
    data.b += geometry_offset;
    data.c += geometry_offset;

    // if (data.a.x < 0 || data.a.x > 255 || data.a.y < 0 || data.a.y > 255 || data.b.x < 0 || data.b.x > 255 || data.b.y < 0 || data.b.y > 255 || data.c.x < 0
    //     || data.c.x > 255 || data.c.y < 0 || data.c.y > 255)
    //     qDebug() << geometry_offset << data.a.x << data.a.y << data.b.x << data.b.y << data.c.x << data.c.y;

    // make sure that we do not remove bits from the coordinates
    assert((uint32_t(data.a.x) & coordinate_bitmask) == uint32_t(data.a.x));
    assert((uint32_t(data.a.y) & coordinate_bitmask) == uint32_t(data.a.y));
    assert((uint32_t(data.b.x) & coordinate_bitmask) == uint32_t(data.b.x));
    assert((uint32_t(data.b.y) & coordinate_bitmask) == uint32_t(data.b.y));
    assert((uint32_t(data.c.x) & coordinate_bitmask) == uint32_t(data.c.x));
    assert((uint32_t(data.c.y) & coordinate_bitmask) == uint32_t(data.c.y));
    // assert(style_layer < ((1u << style_bits) - 1u));
    assert(constants::style_bits + 1 <= available_style_bits); // make sure that the stylebits we need are available here

    packed_data.x = uint32_t(data.a.x) << coordinate_shift1;
    packed_data.x = packed_data.x | ((uint32_t(data.a.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.x = packed_data.x | ((uint32_t(data.b.x) & coordinate_bitmask) << coordinate_shift3);
    packed_data.x = packed_data.x | ((uint32_t(data.b.y) & coordinate_bitmask) << coordinate_shift4);

    packed_data.y = uint32_t(data.c.x) << coordinate_shift1;
    packed_data.y = packed_data.y | ((uint32_t(data.c.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.y = packed_data.y | ((data.style_index << 1) | ((data.is_polygon) ? 1u : 0u));
    // alternative only for neceesary for shader testing
    // packed_data.y = packed_data.y | ((data.style_index << 2) | (((data.should_blend) ? 1u : 0u) << 1) | ((data.is_polygon) ? 1u : 0u));

    return packed_data;
}

VectorLayerData Preprocessor::unpack_data(glm::uvec2 packed_data)
{
    VectorLayerData unpacked_data;

    unpacked_data.a.x = int((packed_data.x & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    unpacked_data.a.y = int((packed_data.x & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);
    unpacked_data.b.x = int((packed_data.x & (coordinate_bitmask << coordinate_shift3)) >> coordinate_shift3);
    unpacked_data.b.y = int((packed_data.x & (coordinate_bitmask << coordinate_shift4)) >> coordinate_shift4);
    unpacked_data.c.x = int((packed_data.y & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    unpacked_data.c.y = int((packed_data.y & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);

    const uint32_t style_and_blend = (packed_data.y & ((1u << available_style_bits) - 1u)) >> 1;
    // NOTE: the should_blend bit is only necessary for the shader -> on cpu side we do not need to separate them
    unpacked_data.style_index = style_and_blend; // >> 1;
    // unpacked_data.should_blend = (style_and_blend & 1u) == 1u;

    unpacked_data.is_polygon = (packed_data.y & 1u) == 1u;

    unpacked_data.a -= geometry_offset;
    unpacked_data.b -= geometry_offset;
    unpacked_data.c -= geometry_offset;

    return unpacked_data;
}

std::pair<uint32_t, uint32_t> Preprocessor::get_split_index(uint32_t index, const std::vector<uint32_t>& polygon_sizes)
{
    // first index test different since we use the previous size in the for loop
    if (index < polygon_sizes[0]) {
        return { 0, index };
    }

    for (uint32_t i = 1; i < polygon_sizes.size(); i++) {
        if (index < polygon_sizes[i]) {
            return { i, index - polygon_sizes[i - 1] };
        }
    }

    // should not happen -> the index does not match a valid polygon point
    assert(false);
    return { 0, 0 };
}

// returns how many triangles have been generated;
size_t Preprocessor::triangulize_earcut(const ClipperPaths& polygon_points, VectorLayerCell* temp_cell, const std::pair<uint32_t, uint32_t>& style_layer)
{
    const auto& indices = mapbox::earcut<uint32_t>(polygon_points);

    // move all polygons sizes to single accumulated array -> so that we can match the index
    std::vector<uint32_t> polygon_sizes;
    polygon_sizes.reserve(polygon_points.size());
    uint32_t previous_size = 0;
    for (size_t i = 0; i < polygon_points.size(); ++i) {
        polygon_sizes.push_back(previous_size + polygon_points[i].size());
        previous_size += polygon_points[i].size();
    }

    // geometry_buffer.reserve(geometry_buffer.size() + indices.size());

    for (size_t i = 0; i < indices.size() / 3; ++i) {
        const auto& ind0 = get_split_index(indices[i * 3 + 0], polygon_sizes);
        const auto& ind1 = get_split_index(indices[i * 3 + 1], polygon_sizes);
        const auto& ind2 = get_split_index(indices[i * 3 + 2], polygon_sizes);

        const auto& p0 = polygon_points[ind0.first][ind0.second];
        const auto& p1 = polygon_points[ind1.first][ind1.second];
        const auto& p2 = polygon_points[ind2.first][ind2.second];

        const auto& data = nucleus::vector_layer::Preprocessor::pack_triangle_data({ { p0.x, p0.y }, { p1.x, p1.y }, { p2.x, p2.y }, style_layer.first, true });

        (*temp_cell).emplace_back(data);
    }

    // returns how many triangles have been generated;
    return indices.size() / 3;
}

void Preprocessor::generate_preprocess_grid()
{
    // make sure that bits we have declared for the cell width are exactly the same amount as the bits we need
    assert(constants::tile_extent / constants::grid_size <= cell_width);

    assert(cell_width < max_cell_width - 50); // make sure that we have plenty of space for big lines

    constexpr auto clipper_margin = 1; // since clipper sometimes returns shapes slightly outside rect -> we need a small margin

    std::vector<PreprocessCell> grid;
    for (int y = 0; y < constants::grid_size; y++) {
        for (int x = 0; x < constants::grid_size; x++) {
            const auto rect = ClipperRect(x * cell_width - constants::aa_border,
                y * cell_width - constants::aa_border,
                (x + 1) * cell_width + constants::aa_border,
                (y + 1) * cell_width + constants::aa_border);
            // for clip lines we are using +0 since we are adding the max_cell_width
            const auto rect_lines = ClipperRect(x * cell_width - geometry_offset + clipper_margin + constants::aa_border,
                y * cell_width - geometry_offset + clipper_margin + constants::aa_border,
                x * cell_width - geometry_offset + max_cell_width - clipper_margin - constants::aa_border,
                y * cell_width - geometry_offset + max_cell_width - clipper_margin - constants::aa_border);

            grid.emplace_back(PreprocessCell { RectClip(rect), RectClipLines(rect_lines), rect, VectorLayerCell(), false });

            // qDebug() << rect.left << rect.top << rect.right << rect.bottom;
        }
    }

    m_preprocess_grid = nucleus::Raster<PreprocessCell>(constants::grid_size, std::move(grid));
}

bool Preprocessor::fully_covers(const ClipperPaths& solution, const ClipperRect& rect)
{
    if (solution.size() != 1 || solution[0].size() != 4)
        return false; // the solution either has holes or does not fully cover the rect since we do not have exactly 4 points

    bool top_left = false;
    bool top_right = false;
    bool bottom_left = false;
    bool bottom_right = false;
    for (int i = 0; i < 4; i++) {
        if (solution[0][i].x <= rect.left && solution[0][i].y <= rect.top)
            top_left = true;
        if (solution[0][i].x >= rect.right && solution[0][i].y <= rect.top)
            top_right = true;
        else if (solution[0][i].x <= rect.left && solution[0][i].y >= rect.bottom)
            bottom_left = true;
        else if (solution[0][i].x >= rect.right && solution[0][i].y >= rect.bottom)
            bottom_right = true;
    }

    // qDebug() << top_left << top_right << bottom_left << bottom_right;

    return top_left && top_right && bottom_left && bottom_right;
}

// checks if a line fully covers the cell with the given line width
// if it does the result is the index of the line index where this is the case, else -1u is returned
size_t Preprocessor::line_fully_covers(const ClipperPaths& solution, float line_width, const ClipperRect& rect)
{

    const auto rect_center = glm::vec2(rect.left + (rect.right - rect.left) / 2, rect.top + (rect.bottom - rect.top) / 2); // TODO could be calculated once
    const auto diagonal_dist_from_center = glm::distance(glm::vec2(rect.left, rect.top), rect_center); // TODO can be calculated outside once

    // qDebug() << "dist_from_center" << diagonal_dist_from_center;

    // no line will ever fill the cell
    if (line_width < diagonal_dist_from_center) {
        return -1u;
    }

    for (size_t i = 0; i < solution.size(); i++) {
        const auto line = solution[i];

        // const auto y_diff = line[1].y - line[0].y;
        // const auto x_diff = line[1].x - line[0].x;

        // // dist to center
        // const auto line_dist_to_furthest_point
        //     = abs(y_diff * rect_center.x - x_diff * rect_center.y + float(line[1].x) * float(line[0].y) - float(line[1].y) * float(line[0].x))
        //         / sqrt(y_diff * y_diff + x_diff * x_diff)
        //     + diagonal_dist_from_center;

        const auto line_dist_to_furthest_point
            = glm::distance(glm::closestPointOnLine(rect_center, glm::vec2(line[0].x, line[0].y), glm::vec2(line[1].x, line[1].y)), rect_center)
            + diagonal_dist_from_center;

        if (line_width >= line_dist_to_furthest_point)
            return size_t(i);
    }

    return -1u;
}

// polygon describe the outer edge of a closed shape
// -> neighbouring vertices form an edge
// last vertex connects to first vertex
void Preprocessor::preprocess_geometry(const VectorLayers& layers)
{

    // TODOs
    // - accelleration grid vector<map> to array<map> since size is given from start
    // -- maybe array<array<vector> -> since we also know the max amount of layer_indices and we can then also create the std::vectors with a default size
    // beforehand

    // reset grid
    for (auto& cell : m_preprocess_grid) {
        cell.is_done = false;

        cell.cell_data.clear();
    }

    m_processed_amount = 0;
    constexpr auto scale = float(constants::grid_size) / float(constants::tile_extent);

    for (auto it = layers.crbegin(); it != layers.crend(); ++it) {

        auto& data = it->second;

        for (size_t i = 0; i < data.size(); ++i) {
            if (data[i].is_polygon) {

                const auto& style_layer = data[i].style_layer;
                const auto& vertices = data[i].vertices;
                const auto& bounds = data[i].bounds;

                m_preprocess_grid.visit(data[i].aabb, [this, &vertices, &bounds, &style_layer](glm::uvec2, PreprocessCell& cell) {
                    if (cell.is_done) {
                        // qDebug() << "cell_done";
                        return;
                    }

                    cell.clipper.ExecuteRepeated(vertices, bounds, &m_clipper_result);

                    if (m_clipper_result.empty())
                        return;

                    if (fully_covers(m_clipper_result, cell.rect)) {
                        cell.is_done = true;

                        // we only need one triangle that covers the whole cell
                        // this triangle is a bit larger to cover the whole cell even with antialiasing
                        const auto& data = nucleus::vector_layer::Preprocessor::pack_triangle_data(
                            { { -cell_width, -cell_width }, { -cell_width, cell_width * 4 }, { cell_width * 4, -cell_width }, style_layer.first, true });

                        cell.cell_data.push_back(data);
                        m_processed_amount++;

                    } else {
                        // anchor clipped paths to cell origin
                        Clipper2Lib::TranslatePathsInPlace(&m_clipper_result, -cell.rect.left, -cell.rect.top);

                        m_processed_amount += triangulize_earcut(m_clipper_result, &cell.cell_data, style_layer);
                    }
                });

            } else {

                const auto line_width = float(m_style_buffer[data[i].style_layer.first >> 1].z) / float(constants::style_precision);

                std::unordered_map<glm::uvec2, std::unordered_set<glm::uvec2, Hasher>, Hasher> cell_list;

                const auto cell_writer = [&cell_list](glm::vec2 pos, const glm::uvec2& data_index) {
                    // if in grid_size bounds and not already present -> than add index to vector
                    if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(constants::grid_size), pos))) {
                        cell_list[glm::uvec2(pos)].insert(data_index);
                    }
                };

                // TODO: according to Task #151 -> we doubled the line width that goes into the acceleration structure because we are looking at tiles
                // that are bigger
                // Nevertheless, we artificially worsened the performance by introducing more cells where a line could be (although it is only there on
                // specific zoom levels)
                // This performance issue will be solved with Task #198 (mipmaps)
                nucleus::utils::rasterizer::rasterize_lines(cell_writer, data[i].vertices, line_width * scale * 2.0, scale);

                const auto& style_layer = data[i].style_layer;
                const auto& vertices = data[i].vertices;

                for (const auto& [cell_pos, indices] : cell_list) {
                    auto& cell = m_preprocess_grid.pixel(cell_pos);

                    if (cell.is_done) {
                        // qDebug() << "cell_done";
                        continue;
                    }

                    auto shapes = ClipperPaths {};
                    for (const auto& index : indices) {

                        // assert(data[i].vertices[0].size() > size_t(index + 1));
                        shapes.emplace_back(ClipperPath { { long(vertices[index.x][index.y].x), long(vertices[index.x][index.y].y) },
                            { long(vertices[index.x][index.y + 1].x), long(vertices[index.x][index.y + 1].y) } });
                    }

                    // ClipperPaths solution;

                    cell.clipper_lines.ExecuteRepeated(shapes, &m_clipper_result);

                    if (m_clipper_result.empty())
                        continue;

                    // TODO this would set cells to is_done if a line fully covers it. but:
                    //  1) this does not fully work for some zoom levels -> investigate why we have some white cells in certain circumstances
                    //      -> currently I think it has something to do with the * 2.0 multiplier when selecting cells, and/or the blending between 2 zooms
                    //  2) this does not really work for the most troublesome zoom levels (like 13/14) since the lines are too small to be much of concern
                    // there
                    //      -> maybe more useful when we reduce the grid cell size and or for optimizing close cells

                    // size_t full_cover_line_index = line_fully_covers(m_clipper_result, line_width, cell.rect);
                    // if (full_cover_line_index != -1u) {
                    //     cell.is_done = true;
                    // }

                    // TODO move translate before clipping and only clip with one single cell
                    // anchor clipped paths to cell origin
                    Clipper2Lib::TranslatePathsInPlace(&m_clipper_result, -cell.rect.left, -cell.rect.top);

                    for (const auto& line : m_clipper_result) {

                        const auto packed_data
                            = nucleus::vector_layer::Preprocessor::pack_line_data({ line[0].x, line[0].y }, { line[1].x, line[1].y }, style_layer.first);
                        cell.cell_data.push_back(packed_data);
                    }
                    m_processed_amount += m_clipper_result.size();
                }
            }
        }

        // qDebug() << geometry_amount;
    }
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
GpuVectorLayerTile Preprocessor::create_gpu_tile()
{
    GpuVectorLayerTile tile;

    uint fitting_cascade_index = uint(-1u);
    for (uint i = 0; i < constants::data_size.size(); i++) {
        if (m_processed_amount <= constants::data_size[i] * constants::data_size[i]) {
            fitting_cascade_index = i;
            break;
        }
    }
    // if assert is triggered -> consider adding a value to constants::data_size
    if (fitting_cascade_index >= constants::data_size.size())
        qDebug() << m_processed_amount << "data does not fit";
    assert(fitting_cascade_index < constants::data_size.size());

    std::vector<uint32_t> acceleration_grid;
    std::vector<glm::u32vec2> geometry_buffer;
    geometry_buffer.reserve(constants::data_size[fitting_cascade_index] * constants::data_size[fitting_cascade_index]);

    // size_t max = 0;

    for (const auto& cell : m_preprocess_grid) {
        // go through every cell
        if (cell.cell_data.size() == 0) {
            acceleration_grid.push_back(0); // no data -> only add an emtpy cell
        } else {

            auto start_offset = geometry_buffer.size();

            geometry_buffer.insert(geometry_buffer.end(), cell.cell_data.cbegin(), cell.cell_data.cend());

            auto cell_size = geometry_buffer.size() - start_offset;

            // we have to add a new element
            // if (cell_size > max)
            //     max = cell_size;

            // TODO enable assert again
            // assert(cell_size< 256); // make sure that we are not removing indices we want to draw
            if (cell_size > 255)
                cell_size = 255; // just cap it to 255 as we currently cant go any higher

            if (cell_size == 0)
                acceleration_grid.push_back(0); // no data -> only add an emtpy cell
            else {
                const auto offset_size = nucleus::utils::bit_coding::u24_u8_to_u32(start_offset, uint8_t(cell_size));
                acceleration_grid.push_back(offset_size);
            }

            // TODO possible performance
            // change geometry buffer to RGBA32UI to save 2 geometries at the same time
            // -> we most likely want to test a few geometries one after another and this might be better
            // also possible to use a 1 bit in the offset_size to indicate the start (to pack the data tightly even if the previous geometry belongs to
            // previous cell)
        }
    }

    // qDebug() << "max: " << max;

    // make sure that the buffer size is still like we expected and resize the data to actual buffer size
    assert(geometry_buffer.size() <= constants::data_size[fitting_cascade_index] * constants::data_size[fitting_cascade_index]);
    geometry_buffer.resize(constants::data_size[fitting_cascade_index] * constants::data_size[fitting_cascade_index], glm::u32vec2(-1u));

    tile.acceleration_grid = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(constants::grid_size, std::move(acceleration_grid)));
    tile.geometry_buffer = std::make_shared<const nucleus::Raster<glm::u32vec2>>(
        nucleus::Raster<glm::u32vec2>(constants::data_size[fitting_cascade_index], std::move(geometry_buffer)));

    tile.buffer_info = fitting_cascade_index;

    return tile;
}

} // namespace nucleus::vector_layer::details
