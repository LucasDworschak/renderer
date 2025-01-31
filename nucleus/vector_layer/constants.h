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
#include <stdint.h>

namespace nucleus::vector_layer::constants {
// sizes are all only one side -> and have to be squared to get the actual amount of data stored in the buffer
// if values here change -> you also need to change them in the shader
constexpr auto grid_size = 64; // 64
constexpr auto data_size = std::array<uint32_t, 4> { 64u, 128u, 256u, 512u }; // needs to be in ascending order
// how many array layers quads per data size
// NOTE: -1u is used to say that we should use the upper limit determined by renderingcontext, if renderingcontext gives us a lower value than set, it is automatically lowered
// IMPORTANT: only set -1u for the first values since those values will be combined to only one array_helper
constexpr auto array_layer_quad_amount = std::array<uint32_t, 4> { -1u, -1u, 256u, 32u };
constexpr auto tile_extent = 4096;
constexpr auto style_buffer_size = 64;

} // namespace nucleus::vector_layer::constants
