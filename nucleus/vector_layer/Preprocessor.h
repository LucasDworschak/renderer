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

#include "Style.h"
#include <radix/hasher.h>

namespace nucleus::vector_layer {

using namespace nucleus::tile;

namespace details {

    /////////////////////////////////////////////
    // constants for data packing/unpacking
    // coordinates 14 bits
    //
    constexpr int all_bits = 32;
    constexpr int coordinate_bits = 14;

    constexpr int style_bits_per_coord = all_bits - (2 * coordinate_bits);

    constexpr int coordinate_shift1 = all_bits - coordinate_bits;
    constexpr int coordinate_shift2 = all_bits - (2 * coordinate_bits);

    constexpr uint32_t coordinate_bitmask = (1u << coordinate_bits) - 1u;
    //
    // end constants for data packing/unpacking
    /////////////////////////////////////////////

    class PointCollectionVec2 : public std::vector<std::vector<glm::vec2>> {
    public:
        using coordinate_type = float;
        static inline bool check_limits = false;
        static inline bool round = false;
        template <class... Args>
        PointCollectionVec2(Args&&... args)
            : std::vector<std::vector<glm::vec2>>(std::forward<Args>(args)...)
        {
        }
    };

    struct Hasher {
        size_t operator()(const std::vector<uint32_t>& seq) const
        {
            size_t seed = 0;
            for (const uint32_t& i : seq) {
                radix::hasher::hash_combine<uint32_t>(seed, i);
            }

            return seed;
        }
    };

    using VectorLayerGrid = std::vector<std::map<uint32_t, uint32_t>>;

    struct VectorLayerCollection {
    public:
        std::vector<glm::u32vec3> vertex_buffer;
        std::vector<uint32_t> index_buffer;
        VectorLayerGrid acceleration_grid;
    };

    struct GeometryData {
        std::vector<std::vector<glm::vec2>> vertices;
        uint32_t extent;
        std::vector<std::pair<uint32_t, uint32_t>> style_and_layer_indices;
        bool is_polygon;
    };

    VectorLayerCollection preprocess_geometry(const std::vector<GeometryData>& data, const std::vector<glm::u32vec4> style_buffer);

    GpuVectorLayerTile create_gpu_tile(const VectorLayerCollection& layer_collection);

    std::vector<GeometryData> parse_tile(tile::Id id, const QByteArray& vector_tile_data, const Style& style);

    std::vector<std::pair<uint32_t, uint32_t>> simplify_styles(std::vector<std::pair<uint32_t, uint32_t>> styles, const std::vector<glm::u32vec4> style_buffer);

    uint32_t pack_index_buffer(uint32_t geometry_index, uint32_t style_index, bool polygon);

    glm::uvec3 pack_triangle_data(glm::vec2 a, glm::vec2 b, glm::vec2 c);
    glm::uvec3 pack_line_data(glm::vec2 a, glm::vec2 b);
    std::tuple<glm::vec2, glm::vec2, glm::vec2> unpack_triangle_data(glm::uvec3 packed);
} // namespace details

GpuVectorLayerTile preprocess(tile::Id id, const QByteArray& vector_tile_data, const Style& style);

GpuVectorLayerTile create_default_gpu_tile();

} // namespace nucleus::vector_layer
