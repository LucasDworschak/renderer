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

const std::shared_ptr<const nucleus::Raster<glm::u32vec2>> Preprocessor::style() { return m_style.visible_styles(); }
bool Preprocessor::update_visible_styles() { return m_style.update_visible_styles(); }

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
    m_processed_amount = 0;
    if (vector_tile_data.isEmpty())
        return {};
    // qDebug() << " z: " << id.zoom_level << "coords: " << id.coords.x << ", " << id.coords.y << (id.scheme == Scheme::Tms);

    // DEBUG polygons
    // const std::vector<std::vector<glm::vec2>> triangle_points = { { glm::vec2(10.5 / 64.0 * constants::grid_size, 30.5 / 64.0 * constants::grid_size),
    //     glm::vec2(30.5 / 64.0 * constants::grid_size, 10.5 / 64.0 * constants::grid_size),
    //     glm::vec2(50.5 / 64.0 * constants::grid_size, 50.5 / 64.0 * constants::grid_size) } };
    // const std::vector<unsigned int> style_indices = { 1 };

    auto tile_data = parse_tile(id, vector_tile_data);

    preprocess_geometry(tile_data, id.zoom_level);

    // return if not a single thing was written
    if (m_processed_amount == 0u)
        return {};

    auto tile = create_gpu_tile();
    tile.id = id;

    return tile;
}

// calculates the area using the shoelace formula (surveyors formula) (https://en.wikipedia.org/wiki/Shoelace_formula)
// NOTE: as we are using this formula to determine the winding order of the polygon -> the area might be negative
// positive area -> exterior polygon ring
// negative area -> interior polygon ring (a hole in the polygon)
float Preprocessor::polygon_area(const ClipperPath& vertices)
{
    float area = 0.0f;

    for (size_t i = 0; i < vertices.size() - 1; ++i) {
        const ClipperPoint& current = vertices[i];
        const ClipperPoint& next = vertices[i + 1];

        area += current.x * next.y - next.x * current.y;
    }

    return area * 0.5f;
}

// goes over all the vertices and detects if a polygon is an exterior or interior(hole) polygon
// if it is an exterior polygon and is not the first polygon checked, this polygon and the following interior polygons are separated.
// the exterior/interior check is done using the polygon_area
// example: e=exterior polygon; i=interior(hole) polygon; paranthesis ()=ClipperPaths
// input: (e i i e e i e i)
// output: (e i i), (e), (e i), (e i)
std::vector<ClipperPaths> Preprocessor::separate_vertex_groups(const ClipperPaths& vertices)
{
    std::vector<ClipperPaths> output;
    ClipperPaths* current_path;

    for (const auto& polygon : vertices) {

        if (output.size() == 0 || polygon_area(polygon) > 0) {
            // we got an exterior polygon
            // create a new entry
            // note the first polygon should always be an exterior
            current_path = &output.emplace_back();
        }
        current_path->push_back(polygon);
    }

    return output;
}

std::vector<StyleLayerIndex> Preprocessor::simplify_styles(
    std::vector<StyleLayerIndex>* style_and_layer_indices, const uint zoom_level, const std::vector<glm::u32vec2>& style_buffer)
{
    // we get multiple styles that may have full opacity and the same width
    // creating render calls for both does not make sense -> we only want to draw the top layer
    // this function simplifys all the styles so that only the styles which actually have a change to be rendered will remain

    // TODO this sort should happen at creation of the vector not here
    // order the styles so that we look at layer in descending order
    std::sort(
        style_and_layer_indices->begin(), style_and_layer_indices->end(), [](StyleLayerIndex a, StyleLayerIndex b) { return a.layer_index > b.layer_index; });
    std::vector<StyleLayerIndex> out_styles;
    int accummulative_opacity = 0;
    float width = 0.0;

    for (const auto& indices : *style_and_layer_indices) {
        const auto style_data_lower = style_buffer[Style::get_style_index(indices.style_index, zoom_level) - 1];
        const auto style_data_higher = style_buffer[Style::get_style_index(indices.style_index, zoom_level)];
        const float lower_width = Style::get_style_width(style_data_lower);
        const float lower_opacity = style_data_lower.x & 255;
        const float higher_width = Style::get_style_width(style_data_higher);
        const bool uses_dashes = Style::uses_dashes(style_data_higher);
        const float higher_opacity = style_data_higher.x & 255;

        // by mixing the lower and higher style -> we get a value that better represents a real world example
        // this is neccessary for e.g. 1 landcover style that stops at z12 and another that starts at z13
        // -> we need to render both, because rendering only one at z13 would falsely represent a fade to alpha 0 between z12 and z13
        // by using z12.5 for the current opacity and width, we can make sure that any can be countered by the fading in the opposite direction
        const float current_width = (lower_width + higher_width) / 2.0;
        const int current_opacity = (lower_opacity + higher_opacity) / 2.0;

        if (current_opacity == 0)
            continue; // we dont care about 0 opacity geometry

        if (width < current_width) {
            // reset opacity
            accummulative_opacity = 0;
            width = current_width;
        }

        if (accummulative_opacity < 255) {
            if (!uses_dashes) // dashes do not count for accummulative opacity
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
            style_and_layer_indices = simplify_styles(&style_and_layer_indices, id.zoom_level, m_style_buffer);

            if (style_and_layer_indices.size() == 0) // no styles found -> we do not visualize it
                continue;

            const auto is_polygon = feature.getType() == mapbox::vector_tile::GeomType::POLYGON;

            const auto scale = (is_polygon) ? constants::scale_polygons : constants::scale_lines;
            const auto cell_scale = float(constants::grid_size) / (float(constants::tile_extent) * scale);
            PointCollectionVec2 geom = feature.getGeometries<PointCollectionVec2>(scale);

            std::vector<GeometryData> all_geometry_data;
            GeometryData* current_geom_data = nullptr;
            if (is_polygon) {
                for (size_t j = 0; j < geom.size(); j++) {

                    bool exterior_ring = polygon_area(geom[j]) > 0.0;

                    if (exterior_ring || j == 0) {
                        // reset the current_geom_data for every exterior ring we find
                        current_geom_data = &all_geometry_data.emplace_back();
                        current_geom_data->aabb = radix::geometry::Aabb2i({ constants::grid_size * 2 * scale, constants::grid_size * 2 * scale },
                            { -constants::grid_size * 2 * scale, -constants::grid_size * 2 * scale });
                    }

                    const auto bound = Clipper2Lib::GetBounds(geom[j]);

                    current_geom_data->aabb.expand_by({ std::floor(float(bound.left) * cell_scale), std::floor(float(bound.top) * cell_scale) });
                    current_geom_data->aabb.expand_by({ std::ceil(float(bound.right) * cell_scale), std::ceil(float(bound.bottom) * cell_scale) });
                    current_geom_data->bounds.push_back(bound);

                    current_geom_data->vertices.push_back(ClipperPath(geom[j].begin(), geom[j].end()));
                }
            } else {
                current_geom_data = &all_geometry_data.emplace_back();
                current_geom_data->vertices = ClipperPaths(geom.begin(), geom.end());
            }

            for (const auto& style_layer : style_and_layer_indices) {
                const auto opacity_lower = m_style_buffer[Style::get_style_index(style_layer.style_index, id.zoom_level) - 1].x & 255;
                const auto opacity_higher = m_style_buffer[Style::get_style_index(style_layer.style_index, id.zoom_level)].x & 255;
                const auto full_opaque = opacity_lower + opacity_higher == (255 + 255);

                for (const auto& geom_data : all_geometry_data) {
                    data[style_layer.layer_index].emplace_back(geom_data.vertices, geom_data.bounds, geom_data.aabb, style_layer, is_polygon, full_opaque);
                }
            }
        }
    }

    return data;
}

glm::u32vec2 Preprocessor::pack_line_data(glm::i64vec2 a, glm::i64vec2 b, uint16_t style_index, bool line_cap0, bool line_cap1)
{

    glm::u32vec2 data = pack_triangle_data({ a, b, b, style_index, false });

    if (line_cap0)
        data.y |= line_cap0_mask;
    if (line_cap1)
        data.y |= line_cap1_mask;

    return data;
}

/*
 * in order to minimize the divergence between lines and polygons, we pack the data as follows:
 *
 * triangle packing:
 * data:
 * x0 | y0 | x1 | y1
 * x2 | y2 | free | is_polygon | style_index
 * bits:
 * 8 | 8 | 8 | 8
 * 8 | 8 | 3 | 1 | 12
 *
 * line packing:
 * data:
 * x0_0 | y0_0 | x1_0 | y1_0
 * x0_1 | y0_1 | x1_1 | y1_1 | free | line_cap0 | line_cap1 | is_polygon | style_index
 * bits:
 * 8 | 8 | 8 | 8
 * 4 | 4 | 4 | 4 | 1 | 1 | 1 | 1 | 12
 * NOTE: a line stores the data in two separate locations x0_0 x0_1
 *       coordinate "_0" are the 8 least significant bits, and "_1" are the more significant bits of the coordinate
 */
glm::u32vec2 Preprocessor::pack_triangle_data(VectorLayerData data)
{

    glm::u32vec2 packed_data;

    if (data.is_polygon) {
        data.a += geometry_offset_polygons;
        data.b += geometry_offset_polygons;
        data.c += geometry_offset_polygons;

        // make sure that we do not remove bits from the coordinates
        assert((uint(data.a.x) & coordinate_bitmask) == uint(data.a.x));
        assert((uint(data.a.y) & coordinate_bitmask) == uint(data.a.y));
        assert((uint(data.b.x) & coordinate_bitmask) == uint(data.b.x));
        assert((uint(data.b.y) & coordinate_bitmask) == uint(data.b.y));
        assert((uint(data.c.x) & coordinate_bitmask) == uint(data.c.x));
        assert((uint(data.c.y) & coordinate_bitmask) == uint(data.c.y));

    } else {
        data.a += geometry_offset_line;
        data.b += geometry_offset_line;

        // make sure that we do not remove bits from the coordinates
        assert((uint(data.a.x) & coordinate_bitmask_lines) == uint(data.a.x));
        assert((uint(data.a.y) & coordinate_bitmask_lines) == uint(data.a.y));
        assert((uint(data.b.x) & coordinate_bitmask_lines) == uint(data.b.x));
        assert((uint(data.b.y) & coordinate_bitmask_lines) == uint(data.b.y));
    }

    // if (data.a.x < 0 || data.a.x > max_cell_width_polygons || data.a.y < 0 || data.a.y > max_cell_width_polygons || data.b.x < 0
    //     || data.b.x > max_cell_width_polygons || data.b.y < 0 || data.b.y > max_cell_width_polygons || data.c.x < 0 || data.c.x > max_cell_width_polygons
    //     || data.c.y < 0 || data.c.y > max_cell_width_polygons)
    //     qDebug() << geometry_offset_polygons << data.a.x << data.a.y << data.b.x << data.b.y << data.c.x << data.c.y;

    // style_bits + 1 since one bit is used to determine if it is a line or a polygon
    assert(constants::style_bits + 1 <= available_style_bits); // make sure that the stylebits we need are available here

    packed_data.x = uint(data.a.x) << coordinate_shift1;
    packed_data.x = packed_data.x | ((uint(data.a.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.x = packed_data.x | ((uint(data.b.x) & coordinate_bitmask) << coordinate_shift3);
    packed_data.x = packed_data.x | ((uint(data.b.y) & coordinate_bitmask) << coordinate_shift4);

    if (data.is_polygon) {
        packed_data.y = uint(data.c.x) << coordinate_shift1;
        packed_data.y = packed_data.y | ((uint(data.c.y) & coordinate_bitmask) << coordinate_shift2);
    } else {
        packed_data.y = ((uint(data.a.x) >> constants::coordinate_bits_polygons) << coordinate_shift1_lines);
        packed_data.y = packed_data.y | ((uint(data.a.y) >> constants::coordinate_bits_polygons) << coordinate_shift2_lines);
        packed_data.y = packed_data.y | ((uint(data.b.x) >> constants::coordinate_bits_polygons) << coordinate_shift3_lines);
        packed_data.y = packed_data.y | ((uint(data.b.y) >> constants::coordinate_bits_polygons) << coordinate_shift4_lines);
    }

    const uint is_polygon = (data.is_polygon ? 1u : 0u) << constants::style_bits;

    packed_data.y = packed_data.y | is_polygon | data.style_index;

    return packed_data;
}

VectorLayerData Preprocessor::unpack_data(glm::uvec2 packed_data)
{
    VectorLayerData unpacked_data;

    unpacked_data.a.x = int((packed_data.x & (coordinate_bitmask_shift1)) >> coordinate_shift1);
    unpacked_data.a.y = int((packed_data.x & (coordinate_bitmask_shift2)) >> coordinate_shift2);
    unpacked_data.b.x = int((packed_data.x & (coordinate_bitmask_shift3)) >> coordinate_shift3);
    unpacked_data.b.y = int((packed_data.x & (coordinate_bitmask_shift4)) >> coordinate_shift4);

    glm::uvec2 c;
    c.x = (packed_data.y & (coordinate_bitmask_shift1)) >> coordinate_shift1;
    c.y = (packed_data.y & (coordinate_bitmask_shift2)) >> coordinate_shift2;

    unpacked_data.style_index = packed_data.y & ((1u << constants::style_bits) - 1u);

    unpacked_data.is_polygon = (packed_data.y & is_polygon_bitmask) != 0u;

    if (unpacked_data.is_polygon) {
        unpacked_data.a -= geometry_offset_polygons;
        unpacked_data.b -= geometry_offset_polygons;
        unpacked_data.c = glm::ivec2(c) - geometry_offset_polygons;
    } else {
        // unpack most significant coordinates of the line and add them to the unpacked lines
        unpacked_data.a.x = unpacked_data.a.x | int(((c.x & (remaining_coordinates_bitmask_shift)) << remaining_coordinate_bits_lines));
        unpacked_data.a.y = unpacked_data.a.y | int(((c.x & remaining_coordinate_bitmask_lines) << constants::coordinate_bits_polygons));
        unpacked_data.b.x = unpacked_data.b.x | int(((c.y & (remaining_coordinates_bitmask_shift)) << remaining_coordinate_bits_lines));
        unpacked_data.b.y = unpacked_data.b.y | int(((c.y & remaining_coordinate_bitmask_lines) << constants::coordinate_bits_polygons));

        unpacked_data.a -= geometry_offset_line;
        unpacked_data.b -= geometry_offset_line;
    }

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
size_t Preprocessor::triangulize_earcut(const ClipperPaths& polygon_points, VectorLayerCell* temp_cell, const StyleLayerIndex& style_layer)
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

        const auto& data
            = nucleus::vector_layer::Preprocessor::pack_triangle_data({ { p0.x, p0.y }, { p1.x, p1.y }, { p2.x, p2.y }, style_layer.style_index, true });

        (*temp_cell).emplace_back(data);
    }

    // returns how many triangles have been generated;
    return indices.size() / 3;
}

void Preprocessor::generate_preprocess_grid()
{
    assert(cell_width_polygons < (1 << (constants::coordinate_bits_polygons)));
    assert(cell_width_lines < (1 << (constants::coordinate_bits_lines)));
    assert(cell_width_lines < max_cell_width_line);

    std::vector<PreprocessCell> grid;
    for (int y = 0; y < constants::grid_size; y++) {
        for (int x = 0; x < constants::grid_size; x++) {
            const auto rect = ClipperRect(x * cell_width_polygons - cell_width_polygons * constants::aa_border,
                y * cell_width_polygons - cell_width_polygons * constants::aa_border,
                (x + 1) * cell_width_polygons + cell_width_polygons * constants::aa_border,
                (y + 1) * cell_width_polygons + cell_width_polygons * constants::aa_border);
            const auto rect_lines = ClipperRect(x * cell_width_lines - cell_width_lines * constants::aa_border,
                y * cell_width_lines - cell_width_lines * constants::aa_border,
                (x + 1) * cell_width_lines + cell_width_lines * constants::aa_border,
                (y + 1) * cell_width_lines + cell_width_lines * constants::aa_border);

            grid.emplace_back(PreprocessCell { RectClip(rect), rect, rect_lines, VectorLayerCell(), false });

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
void Preprocessor::preprocess_geometry(const VectorLayers& layers, const uint zoom_level)
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

    for (auto it = layers.crbegin(); it != layers.crend(); ++it) {

        auto& data = it->second;

        for (size_t i = 0; i < data.size(); ++i) {
            if (data[i].is_polygon) {

                const auto& style_layer = data[i].style_layer;
                const auto& vertices = data[i].vertices;
                const auto& bounds = data[i].bounds;
                const auto& check_fully_covers = data[i].full_opaque;

                m_preprocess_grid.visit(data[i].aabb, [this, &vertices, &bounds, &style_layer, &check_fully_covers](glm::uvec2, PreprocessCell& cell) {
                    if (cell.is_done) {
                        // qDebug() << "cell_done";
                        return;
                    }

                    cell.clipper.ExecuteRepeated(vertices, bounds, &m_clipper_result);

                    if (m_clipper_result.empty())
                        return;

                    if (check_fully_covers && fully_covers(m_clipper_result, cell.rect_polygons)) {
                        cell.is_done = true;

                        // we are packing two triangles -> packing only one triangle would result in problems with multisample antialiasing
                        // we tried packing only one triangle but artificially scaling it up in the shader if a "is_full" flag was transmitted, but this causes
                        // more performance problems

                        const glm::ivec2 a = { -geometry_offset_polygons, -geometry_offset_polygons };
                        const glm::ivec2 b = { max_cell_width_polygons - geometry_offset_polygons - 1, -geometry_offset_polygons };
                        const glm::ivec2 c = { -geometry_offset_polygons, max_cell_width_polygons - geometry_offset_polygons - 1 };
                        const glm::ivec2 d = { max_cell_width_polygons - geometry_offset_polygons - 1, max_cell_width_polygons - geometry_offset_polygons - 1 };

                        const auto& data1 = nucleus::vector_layer::Preprocessor::pack_triangle_data({ a, b, c, style_layer.style_index, true });
                        const auto& data2 = nucleus::vector_layer::Preprocessor::pack_triangle_data({ d, b, c, style_layer.style_index, true });

                        cell.cell_data.push_back(data1);
                        cell.cell_data.push_back(data2);
                        m_processed_amount += 2;

                    } else {
                        // anchor clipped paths to cell origin
                        Clipper2Lib::TranslatePathsInPlace(&m_clipper_result,
                            -cell.rect_polygons.left - cell_width_polygons * constants::aa_border,
                            -cell.rect_polygons.top - cell_width_polygons * constants::aa_border);

                        const auto vertex_groups = separate_vertex_groups(m_clipper_result);

                        for (const auto& vertices : vertex_groups)
                            m_processed_amount += triangulize_earcut(vertices, &cell.cell_data, style_layer);
                    }
                });

            } else {
                constexpr auto scale = float(constants::grid_size) / (float(constants::tile_extent) * constants::scale_lines);

                float line_width = 0;
                // use the line width of the previous style

                line_width
                    = Style::get_style_width(m_style_buffer[Style::get_style_index(data[i].style_layer.style_index, zoom_level) - 1]) + constants::aa_lines;

                std::unordered_map<glm::uvec2, std::unordered_set<glm::uvec2, Hasher>, Hasher> cell_list;

                const auto cell_writer = [&cell_list](glm::vec2 pos, const glm::uvec2& data_index) {
                    // if in grid_size bounds and not already present -> than add index to vector
                    if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(constants::grid_size), pos))) {
                        cell_list[glm::uvec2(pos)].insert(data_index);
                    }
                };

                nucleus::utils::rasterizer::rasterize_lines(cell_writer, data[i].vertices, line_width * scale * 1.0 * constants::scale_lines, scale);

                const auto& style_layer = data[i].style_layer;
                const auto& vertices = data[i].vertices;

                for (const auto& [cell_pos, indices] : cell_list) {
                    auto& cell = m_preprocess_grid.pixel(cell_pos);

                    if (cell.is_done) {
                        // qDebug() << "cell_done";
                        continue;
                    }

                    for (const auto& index : indices) {
                        const auto packed_data = nucleus::vector_layer::Preprocessor::pack_line_data(
                            { long(vertices[index.x][index.y].x) - cell.rect_lines.left - cell_width_lines * constants::aa_border,
                                long(vertices[index.x][index.y].y) - cell.rect_lines.top - cell_width_lines * constants::aa_border },
                            { long(vertices[index.x][index.y + 1].x) - cell.rect_lines.left - cell_width_lines * constants::aa_border,
                                long(vertices[index.x][index.y + 1].y) - cell.rect_lines.top - cell_width_lines * constants::aa_border },
                            style_layer.style_index,
                            index.y == 0,
                            index.y == vertices[index.x].size() - 2); // -2 because the index we get does not use the last element of a line segment
                        cell.cell_data.push_back(packed_data);
                        m_processed_amount++;
                    }
                }
            }
        }

        // qDebug() << m_processed_amount;
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

            // TODO if we stay with 512x512 geometry buffer size -> we can add two bits to size and remove those from offset
            // -> this way there is no need to artificially cap the geometries to 255
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
