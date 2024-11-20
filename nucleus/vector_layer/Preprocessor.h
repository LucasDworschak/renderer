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

#pragma once

#include <vector>

#include <glm/glm.hpp>

#include "nucleus/tile/types.h"

namespace nucleus::vector_layer {

using namespace nucleus::tile;

union Line {
    // sizeof(glm::vec2) == 64 bits
    // sizeof(uint32_t) == 32 bits
    // -> 64 + 64 + 32 = 5 uint32_t
    struct {
        glm::vec2 line_start_vertex;
        glm::vec2 line_end_vertex;

        uint32_t style_index;

    } data;

    uint32_t packed[5];

    Line() = default;
    Line(glm::vec2 line_start_vertex, glm::vec2 line_end_vertex, uint32_t style_index)
    {
        data.line_start_vertex = line_start_vertex;
        data.line_end_vertex = line_end_vertex;
        data.style_index = style_index;
    }
};

union Triangle {
    // sizeof(glm::vec2) == 64 bits
    // sizeof(uint32_t) == 32 bits
    // -> 64 + 64 + 64 + 32 = 7 uint32_t
    struct {
        glm::vec2 top_vertex;
        glm::vec2 middle_vertex;
        glm::vec2 bottom_vertex;

        uint32_t style_index;

    } data;

    uint32_t packed[7];

    Triangle() = default;
    Triangle(glm::vec2 top_vertex, glm::vec2 middle_vertex, glm::vec2 bottom_vertex, uint32_t style_index)
    {
        data.top_vertex = top_vertex;
        data.middle_vertex = middle_vertex;
        data.bottom_vertex = bottom_vertex;
        data.style_index = style_index;
    }
};

using VectorLayerGrid = std::vector<std::unordered_set<uint32_t>>;

struct VectorLayerLineCollection {
public:
    std::vector<Line> lines;
    VectorLayerGrid cell_to_data;
};

struct VectorLayerTriangleCollection {
public:
    std::vector<Triangle> triangles;
    VectorLayerGrid cell_to_data;
};

// helpers for catch2
inline std::ostream& operator<<(std::ostream& os, const glm::vec2& v) { return os << "{ " << v.x << ", " << v.y << " }"; }

inline std::ostream& operator<<(std::ostream& os, const Triangle& t)
{
    return os << "{ top: " << t.data.top_vertex << ", middle: " << t.data.middle_vertex << ", bottom: " << t.data.bottom_vertex << ", style: " << std::to_string(t.data.style_index) << " }";
}

inline std::ostream& operator<<(std::ostream& os, const Line& t)
{
    return os << "{ start: " << t.data.line_start_vertex << ", end: " << t.data.line_end_vertex << ", style: " << std::to_string(t.data.style_index) << " }";
}

inline bool operator==(const Triangle& t1, const Triangle& t2)
{
    return t1.data.top_vertex == t2.data.top_vertex && t1.data.middle_vertex == t2.data.middle_vertex && t1.data.bottom_vertex == t2.data.bottom_vertex && t1.data.style_index == t2.data.style_index;
}

inline bool operator==(const Line& l1, const Line& l2)
{
    return l1.data.line_start_vertex == l2.data.line_start_vertex && l1.data.line_end_vertex == l2.data.line_end_vertex && l1.data.style_index == l2.data.style_index;
}

class Preprocessor {

public:
    Preprocessor(nucleus::tile::Id id);
    GpuVectorLayerTile preprocess(const tile::Data data);

    VectorLayerTriangleCollection preprocess_triangles(const std::vector<std::vector<glm::vec2>> polygons, const std::vector<unsigned int> style_indices);
    VectorLayerLineCollection preprocess_lines(const std::vector<std::vector<glm::vec2>> lines, const std::vector<unsigned int> style_indices);
    nucleus::Raster<uint8_t> visualize_grid(const VectorLayerGrid& grid);

private:
    // VectorLayer m_processed_tile;

    const glm::uvec2 m_grid_size = { 64, 64 };
    std::vector<int> m_x_values_per_y_step;
    int m_tile_up_direction;
    tile::SrsBounds m_tile_bounds;

    void create_triangles(VectorLayerTriangleCollection& triangle_collection, const std::vector<glm::vec2> polygon_points, unsigned int style_index);
    void create_lines(VectorLayerLineCollection& line_collection, const std::vector<glm::vec2> line_points, unsigned int style_index);

    Triangle create_ordered_triangle(glm::vec2 triangle_vertex_a, glm::vec2 triangle_vertex_b, glm::vec2 triangle_vertex_c, unsigned int style_index);

    inline std::pair<glm::vec2, int> calculate_dda_steps(const glm::vec2 line);

    void dda_line(VectorLayerGrid& grid, const glm::vec2 origin, const glm::vec2 line, const glm::vec2 thickness_normal, unsigned int data_index, int fill_direction, bool is_triangle);
    void dda_triangle(VectorLayerTriangleCollection& triangle_collection, unsigned int triangle_index, float thickness);
    void add_end_cap(VectorLayerGrid& grid, const glm::vec2 position, unsigned int data_index, float thickness);
    void write_to_cell(VectorLayerGrid& grid, const glm::vec2 current_position, unsigned int data_index);
};
} // namespace nucleus::vector_layer
