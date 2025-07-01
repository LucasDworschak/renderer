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

#include <array>
#include <cmath>
#include <stdint.h>

#include <glm/glm.hpp>

namespace nucleus::vector_layer::constants {
// sizes are all only one side -> and have to be squared to get the actual amount of data stored in the buffer
// if values here change -> you also need to change them in the shader
constexpr auto grid_size = 64;
// constexpr auto grid_size = 128;
constexpr auto data_size = std::array<uint32_t, 4> { 64u, 128u, 256u, 512u }; // needs to be in ascending order
// constexpr auto data_size = std::array<uint32_t, 4> { 128u, 64u, 512u, 1024u }; // needs to be in ascending order
// how many array layers tiles per data size
// NOTE: -1u is used to say that we should use the upper limit determined by renderingcontext, if renderingcontext gives us a lower value than set, it is
// automatically lowered IMPORTANT: only set -1u for the first values since those values will be combined to only one array_helper
constexpr auto array_layer_tile_amount = std::array<uint32_t, 4> { -1u, -1u, 1524u, 256u };
constexpr auto tile_extent = 1024;
constexpr auto scale_polygons = 4.0;
constexpr auto scale_lines = 1.0;
constexpr auto mipmap_levels = 4;

// if you change one of the following setting you also need to change the other (that they match is asserted in style.cpp)
// one bit is used to signal if it should blend with next style or not
constexpr auto style_bits = 13;
constexpr auto style_buffer_size = 64;

// number is used to multiply float numbers in style to get int values for storage ( e.g. style_precision:100 means 1.253 -> 125 -> 1.250)
constexpr uint16_t style_precision = 100;
// ensure that we know how large a line can get
constexpr uint16_t max_line_width = 48;

constexpr float line_width_multiplier = 1.25;
constexpr glm::uvec2 style_zoom_range = glm::uvec2(8u, 18u);

constexpr float small_line_px = 0.5; // how many px on a 256x256 resolution
constexpr unsigned small_line_zoom_threshold = 4; // until what zoom level do we use the small line simplification

// there is only a limited amount of style expressions that are tested by the StyleExpression code
// -> we can store all the values in a map and have an even faster comparison between values
// openstreetmap 14
// qwant 10
// osm-bright 11
constexpr auto max_style_expression_keys = 14;

// data unpacking/packing constants
constexpr int all_bits = 32; // per output channel
constexpr int coordinate_bits_polygons = 8;
constexpr int coordinate_bits_lines = 12;
constexpr int aa_border = 2;

// array helper constants
constexpr int array_helper_all_bits = 16;
constexpr int array_helper_buffer_info_bits = 2;

} // namespace nucleus::vector_layer::constants
