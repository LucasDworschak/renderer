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

    struct PolygonData {
        std::vector<glm::vec2> vertices;
        std::vector<glm::ivec2> edges;
    };

    using VectorLayerGrid = std::vector<std::unordered_set<uint32_t>>;

    struct VectorLayerCollection {
    public:
        std::vector<uint32_t> vertex_buffer;
        std::vector<uint32_t> index_buffer;
        VectorLayerGrid acceleration_grid;
    };

    struct TempDataHolder {
        std::vector<PolygonData> polygons;
        uint32_t extent;
        std::vector<size_t> polygon_styles;
        std::vector<PolygonData> lines;
        std::vector<size_t> line_styles;
    };

    VectorLayerCollection preprocess_triangles(const TempDataHolder data);
    VectorLayerCollection preprocess_lines(const TempDataHolder data);

    GpuVectorLayerTile create_gpu_tile(const VectorLayerCollection& triangle_collection, const VectorLayerCollection& line_collection);

    TempDataHolder parse_tile(tile::Id id, const QByteArray& vector_tile_data, const Style& style);

    /**
     * 96 bits -> rgb32UI
     * 2*3*13=78 bits for all coordinate values
     * 14 bits for style
     * 78+14 bits = 92 -> 4 bits remain
     */
    glm::uvec3 pack_triangle_data(glm::vec2 a, glm::vec2 b, glm::vec2 c, uint32_t style_index);
    std::tuple<glm::vec2, glm::vec2, glm::vec2, uint32_t> unpack_triangle_data(glm::uvec3 packed);
} // namespace details

GpuVectorLayerTile preprocess(tile::Id id, const QByteArray& vector_tile_data, const Style& style);

} // namespace nucleus::vector_layer
