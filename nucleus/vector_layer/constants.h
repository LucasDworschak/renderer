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

namespace nucleus::vector_layer::constants {

// sizes are all only one side -> and have to be squared to get the actual amount of data stored in them

constexpr auto grid_size = 16; // 64
// 280=7*40 -> if 7 is part of the size (regarding prime factorization), we dont have to worry for corner cases where part of data is on another tile (since the data currently is also 7 vars long)
// constexpr auto data_size = 280;
constexpr auto data_size = 512;
constexpr auto style_data_size = 64;

} // namespace nucleus::vector_layer::constants
