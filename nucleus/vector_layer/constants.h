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
constexpr auto tile_extent = 4096;

// if you change one of the following setting you also need to change the other (that they match is asserted in style.cpp)
// one bit is used to signal if it should blend with next style or not
constexpr auto style_bits = 13;
constexpr auto style_buffer_size = 64;

// number is used to multiply float numbers in style to get int values for storage ( e.g. style_precision:100 means 1.253 -> 125 -> 1.250)
constexpr uint16_t style_precision = 100;

constexpr uint8_t line_width_multiplier = 5;
constexpr glm::uvec2 style_zoom_range = glm::uvec2(8u, 18u);
constexpr uint8_t style_zoom_blend_steps = 4u; // how many zoom steps below min_zoom should be created to allow blending

constexpr float small_line_px = 2.0; // how many px on a 256x256 resolution
constexpr unsigned small_line_zoom_threshold = 16; // until what zoom level do we use the small line simplification

} // namespace nucleus::vector_layer::constants
