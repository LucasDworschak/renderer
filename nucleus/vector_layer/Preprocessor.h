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

#include <clipper2/clipper.h>
#include <glm/glm.hpp>
#include <radix/hasher.h>

#include "nucleus/tile/types.h"
#include "nucleus/vector_layer/Style.h"
#include "nucleus/vector_layer/constants.h"

namespace nucleus::vector_layer {

using namespace nucleus::tile;

/////////////////////////////////////////////
// constants for data packing/unpacking

constexpr int available_style_bits = constants::all_bits - (2 * constants::coordinate_bits_polygons);

constexpr uint coordinate_bitmask = (1u << constants::coordinate_bits_polygons) - 1u;
constexpr uint coordinate_bitmask_lines = (1u << constants::coordinate_bits_lines) - 1u;
constexpr int remaining_coordinate_bits_lines = constants::coordinate_bits_lines - constants::coordinate_bits_polygons;
constexpr uint remaining_coordinate_bitmask_lines = (1u << remaining_coordinate_bits_lines) - 1u;
constexpr uint is_polygon_bitmask = (1u << constants::style_bits);

constexpr int coordinate_shift1 = constants::all_bits - constants::coordinate_bits_polygons;
constexpr int coordinate_shift2 = constants::all_bits - (2 * constants::coordinate_bits_polygons);
constexpr int coordinate_shift3 = constants::all_bits - (3 * constants::coordinate_bits_polygons);
constexpr int coordinate_shift4 = constants::all_bits - (4 * constants::coordinate_bits_polygons);

constexpr uint coordinate_bitmask_shift1 = coordinate_bitmask << coordinate_shift1;
constexpr uint coordinate_bitmask_shift2 = coordinate_bitmask << coordinate_shift2;
constexpr uint coordinate_bitmask_shift3 = coordinate_bitmask << coordinate_shift3;
constexpr uint coordinate_bitmask_shift4 = coordinate_bitmask << coordinate_shift4;
constexpr uint remaining_coordinates_bitmask_shift = remaining_coordinate_bitmask_lines << remaining_coordinate_bits_lines;

// for packed.y -> we need to divide coordinate_bits_polygons by 2
constexpr int coordinate_shift1_lines = constants::all_bits - int(0.5 * constants::coordinate_bits_polygons);
constexpr int coordinate_shift2_lines = constants::all_bits - int(1.0 * constants::coordinate_bits_polygons);
constexpr int coordinate_shift3_lines = constants::all_bits - int(1.5 * constants::coordinate_bits_polygons);
constexpr int coordinate_shift4_lines = constants::all_bits - int(2.0 * constants::coordinate_bits_polygons);

constexpr int cell_width_polygons = int((constants::tile_extent * constants::scale_polygons) / constants::grid_size);
constexpr int cell_width_lines = int((constants::tile_extent * constants::scale_lines) / constants::grid_size);

constexpr int max_cell_width_polygons = (1 << (constants::coordinate_bits_polygons));
constexpr int geometry_offset_polygons = (max_cell_width_polygons - cell_width_polygons) / 2;
constexpr int max_cell_width_line = (1 << (constants::coordinate_bits_lines));
constexpr int geometry_offset_line = (max_cell_width_line - cell_width_lines) / 2;

constexpr uint line_cap0_mask = 1u << (constants::style_bits + 1);
constexpr uint line_cap1_mask = 1u << (constants::style_bits + 2);

//
// end constants for data packing/unpacking
/////////////////////////////////////////////

struct VectorLayerData {
    glm::ivec2 a;
    glm::ivec2 b;
    glm::ivec2 c;

    uint32_t style_index;
    bool is_polygon;
};

using ClipperResolution = int16_t;
using ClipperPoint = Clipper2Lib::Point<ClipperResolution>;
using ClipperPath = Clipper2Lib::Path<ClipperResolution>;
using ClipperPaths = Clipper2Lib::Paths<ClipperResolution>;
using ClipperRect = Clipper2Lib::Rect<ClipperResolution>;
using RectClip = Clipper2Lib::RectClip<ClipperResolution>;
using RectClipLines = Clipper2Lib::RectClipLines<ClipperResolution>;

class PointCollectionVec2 : public ClipperPaths {
public:
    using coordinate_type = ClipperResolution;
    static inline bool check_limits = false;
    static inline bool round = false;
    template <class... Args>
    PointCollectionVec2(Args&&... args)
        : ClipperPaths(std::forward<Args>(args)...)
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

    size_t operator()(const glm::uvec2& coords) const
    {
        size_t seed = 0;
        radix::hasher::hash_combine<uint32_t>(seed, coords.x);
        radix::hasher::hash_combine<uint32_t>(seed, coords.y);

        return seed;
    }
};

struct GeometryData {
    ClipperPaths vertices;
    std::vector<ClipperRect> bounds;

    radix::geometry::Aabb2i aabb;
    StyleLayerIndex style_layer;
    bool is_polygon;
    bool full_opaque;
};

using VectorLayers = std::map<uint32_t, std::vector<GeometryData>>;
using VectorLayerCell = std::vector<glm::u32vec2>;

struct PreprocessCell {
    RectClip clipper;
    ClipperRect rect_polygons;
    ClipperRect rect_lines;
    VectorLayerCell cell_data;
    bool is_done;
};

class Preprocessor {
public:
    Preprocessor(Style&& style);

    GpuVectorLayerTile preprocess(tile::Id id, const QByteArray& vector_tile_data);

    VectorLayers parse_tile(tile::Id id, const QByteArray& vector_tile_data);
    void preprocess_geometry(const VectorLayers& layers, const uint zoom_level);
    GpuVectorLayerTile create_gpu_tile();

    static glm::u32vec2 pack_triangle_data(VectorLayerData data);
    static glm::u32vec2 pack_line_data(glm::i64vec2 a, glm::i64vec2 b, uint16_t style_layer, bool line_cap0, bool line_cap1);
    static VectorLayerData unpack_data(glm::uvec2 packed_data);

    static bool fully_covers(const ClipperPaths& solution, const ClipperRect& rect);
    static size_t line_fully_covers(const ClipperPaths& solution, float line_width, const ClipperRect& rect);

    static GpuVectorLayerTile create_default_gpu_tile();

    static std::vector<StyleLayerIndex> simplify_styles(
        std::vector<StyleLayerIndex>* styles, const uint zoom_level, const std::vector<glm::u32vec2>& style_buffer);

    const std::shared_ptr<const nucleus::Raster<glm::u32vec2>> style();
    bool update_visible_styles();
    size_t processed_amount();

private:
    std::pair<uint32_t, uint32_t> get_split_index(uint32_t index, const std::vector<uint32_t>& polygon_sizes);

    size_t triangulize_earcut(const ClipperPaths& polygon_points, VectorLayerCell* temp_cell, const StyleLayerIndex& style_layer);

    void generate_preprocess_grid();

    Style m_style;
    const std::vector<glm::u32vec2> m_style_buffer;

    nucleus::Raster<PreprocessCell> m_preprocess_grid;

    size_t m_processed_amount;

    ClipperPaths m_clipper_result; // temporary data. we hold it here to avoid allocations.
};

} // namespace nucleus::vector_layer
