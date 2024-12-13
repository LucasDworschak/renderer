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

namespace nucleus::vector_layer {

using namespace nucleus::tile;

namespace details {

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

    constexpr glm::uvec2 grid_size = { 64, 64 };

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

    VectorLayerCollection preprocess_triangles(const PointCollectionVec2 polygons, const std::vector<unsigned int> style_indices);
    VectorLayerCollection preprocess_lines(const PointCollectionVec2 lines, const std::vector<unsigned int> style_indices);

    GpuVectorLayerTile create_gpu_tile(const VectorLayerCollection& triangle_collection, const VectorLayerCollection& line_collection);

    PointCollectionVec2 parse_tile(tile::Id id, const QByteArray& vector_tile_data, const Style& style);
} // namespace details

GpuVectorLayerTile preprocess(tile::Id id, const QByteArray& vector_tile_data, const Style& style);

} // namespace nucleus::vector_layer
