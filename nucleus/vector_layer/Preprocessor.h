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

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <glm/glm.hpp>

#include "nucleus/vector_tile/types.h"

namespace nucleus::vector_layer {

using namespace nucleus::vector_tile;

// struct DataHolder {
//     std::vector<glm::vec2> m_vertices;
//     std::vector<Triangle> m_triangles;

//     std::vector<std::unordered_set<uint8_t>> m_cell_to_triangle;
// };

// TODO refactor to use vectorlayer instead. Vector layer is passed in in the preprocess function. it holds all the vertice and triangle data
// VectorLayer

class Preprocessor {

public:
    Preprocessor();
    void preprocess(const nucleus::tile::Id tile_id, const std::vector<std::vector<glm::vec2>> polygons, const std::vector<unsigned int> style_indices);
    void visualize_grid();

private:
    VectorLayer m_processed_tile;
    // std::unordered_map<tile::Id, std::unordered_set<unsigned int>, tile::Id::Hasher> m_triangle_neighbours;

    const glm::uvec2 m_grid_size = { 64, 64 };
    std::vector<int> m_x_values_per_y_step;
    int m_tile_up_direction;
    tile::SrsBounds m_tile_bounds;

    void create_triangles(const std::vector<glm::vec2> polygon_points, unsigned int style_index);

    Triangle create_ordered_triangle(unsigned int triangle_index_a, unsigned int triangle_index_b, unsigned int triangle_index_c, unsigned int style_index);
    void dda_line(const glm::vec2 top_point, const glm::vec2 bottom_point, unsigned int data_index, int fill_direction);
    void dda_triangle(unsigned int triangle_index);

    void write_to_cell(glm::vec2 current_position, std::array<glm::vec2, 3> points, unsigned int data_index, float distance);
    // void write_to_cell(glm::vec2 position, unsigned int data_index);
};
} // namespace nucleus::vector_layer
