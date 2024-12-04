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

// union Line {
//     // sizeof(glm::vec2) == 64 bits
//     // sizeof(uint32_t) == 32 bits
//     // -> 64 + 64 + 32 = 5 uint32_t
//     struct {
//         glm::vec2 line_start_vertex;
//         glm::vec2 line_end_vertex;

//         uint32_t style_index;

//     } data;

//     uint32_t packed[5];

//     Line() = default;
//     Line(glm::vec2 line_start_vertex, glm::vec2 line_end_vertex, uint32_t style_index)
//     {
//         data.line_start_vertex = line_start_vertex;
//         data.line_end_vertex = line_end_vertex;
//         data.style_index = style_index;
//     }
// };

// union Triangle {
//     // sizeof(glm::vec2) == 64 bits
//     // sizeof(uint32_t) == 32 bits
//     // -> 64 + 64 + 64 + 32 = 7 uint32_t
//     struct {
//         glm::vec2 top_vertex;
//         glm::vec2 middle_vertex;
//         glm::vec2 bottom_vertex;

//         uint32_t style_index;

//     } data;

//     uint32_t packed[7];

//     Triangle() = default;
//     Triangle(glm::vec2 top_vertex, glm::vec2 middle_vertex, glm::vec2 bottom_vertex, uint32_t style_index)
//     {
//         data.top_vertex = top_vertex;
//         data.middle_vertex = middle_vertex;
//         data.bottom_vertex = bottom_vertex;
//         data.style_index = style_index;
//     }
// };

template <class T> inline void hash_combine(std::size_t& seed, T const& v) { seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2); }

struct Hasher {
    size_t operator()(const std::unordered_set<uint32_t>& seq) const
    {
        size_t seed = 0;
        for (const uint32_t& i : seq) {
            hash_combine<uint32_t>(seed, i);
        }

        return seed;
    }
};

using VectorLayerGrid = std::vector<std::unordered_set<uint32_t>>;

struct VectorLayerCollection {
public:
    std::vector<uint32_t> data;
    std::vector<uint32_t> index_bridge;
    std::vector<uint32_t> cell_to_index_bridge;
    VectorLayerGrid cell_to_temp;
};

class Preprocessor {

public:
    Preprocessor(nucleus::tile::Id id);
    GpuVectorLayerTile preprocess(const tile::Data data);

    VectorLayerCollection preprocess_triangles(const std::vector<std::vector<glm::vec2>> polygons, const std::vector<unsigned int> style_indices);
    VectorLayerCollection preprocess_lines(const std::vector<std::vector<glm::vec2>> lines, const std::vector<unsigned int> style_indices);

    void condense_data(VectorLayerCollection& triangle_collection, VectorLayerCollection& line_collection, GpuVectorLayerTile& tile);

private:
    const glm::uvec2 m_grid_size = { 64, 64 };
    std::vector<int> m_x_values_per_y_step;
    int m_tile_up_direction;
    tile::SrsBounds m_tile_bounds;


};
} // namespace nucleus::vector_layer
