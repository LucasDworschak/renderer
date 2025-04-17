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

#include <clipper2/clipper.h>

namespace nucleus::vector_layer {

using namespace nucleus::tile;

namespace details {

    /////////////////////////////////////////////
    // constants for data packing/unpacking
    // coordinates 14 bits
    //
    constexpr int all_bits = 32;
    constexpr int coordinate_bits = 8;
    constexpr int coordinate_bits_buffer = 2; // those bits can be used to address coordinates outside cell

    constexpr int available_style_bits = all_bits - (2 * coordinate_bits);

    constexpr int coordinate_shift1 = all_bits - coordinate_bits;
    constexpr int coordinate_shift2 = all_bits - (2 * coordinate_bits);
    constexpr int coordinate_shift3 = all_bits - (3 * coordinate_bits);
    constexpr int coordinate_shift4 = all_bits - (4 * coordinate_bits);

    constexpr uint32_t coordinate_bitmask = (1u << coordinate_bits) - 1u;
    constexpr int32_t geometry_offset = coordinate_bitmask >> coordinate_bits_buffer;
    //
    // end constants for data packing/unpacking
    /////////////////////////////////////////////

    struct VectorLayerData {
        glm::ivec2 a;
        glm::ivec2 b;
        glm::ivec2 c;

        uint32_t style_index;
        // bool should_blend;
        bool is_polygon;
    };

    struct ClipperRect {
        Clipper2Lib::RectClip64 clipper;
        Clipper2Lib::Rect64 rect;
        bool is_done;
    };

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

    using VectorLayerCell = std::map<uint32_t, std::vector<glm::u32vec2>>;
    using VectorLayerGrid = std::vector<VectorLayerCell>;

    struct VectorLayerMeta {
        VectorLayerGrid temp_grid;
        size_t geometry_amount;
    };

    struct GeometryData {
        std::vector<std::vector<glm::vec2>> vertices;
        uint32_t extent;
        std::vector<std::pair<uint32_t, uint32_t>> style_and_layer_indices;
        bool is_polygon;
    };

    VectorLayerMeta preprocess_geometry(const std::vector<GeometryData>& data, const std::vector<glm::u32vec4> style_buffer);

    nucleus::Raster<ClipperRect> generate_clipper2_grid(unsigned steps);
    std::pair<uint32_t, uint32_t> get_split_index(uint32_t index, const std::vector<uint32_t>& polygon_sizes);

    size_t triangulize_earcut(
        const Clipper2Lib::Paths64& polygon_points, VectorLayerCell* temp_cell, const std::vector<std::pair<uint32_t, uint32_t>>& style_and_layer_indices);

    GpuVectorLayerTile create_gpu_tile(const VectorLayerMeta& meta);

    std::vector<GeometryData> parse_tile(tile::Id id, const QByteArray& vector_tile_data, const Style& style);

    std::vector<std::pair<uint32_t, uint32_t>> simplify_styles(std::vector<std::pair<uint32_t, uint32_t>> styles, const std::vector<glm::u32vec4> style_buffer);

    glm::u32vec2 pack_triangle_data(VectorLayerData data);
    glm::u32vec2 pack_line_data(glm::i64vec2 a, glm::i64vec2 b, uint16_t style_layer);
    VectorLayerData unpack_triangle_data(glm::uvec2 packed_data);

} // namespace details

GpuVectorLayerTile preprocess(tile::Id id, const QByteArray& vector_tile_data, const Style& style);

GpuVectorLayerTile create_default_gpu_tile();

} // namespace nucleus::vector_layer
