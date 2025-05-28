/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2025 Lucas Dworschak
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

#include <nucleus/Raster.h>
#include <nucleus/tile/GpuArrayHelper.h>
#include <nucleus/tile/types.h>
#include <nucleus/vector_layer/constants.h>

namespace nucleus::vector_layer {

// calculates the first index of array_layer_quad_amount that is not -1
constexpr static size_t custom_array_layer_index()
{
    for (size_t i = 0; i < constants::array_layer_tile_amount.size(); i++) {
        if (constants::array_layer_tile_amount[i] != -1u)
            return i;
    }
    return constants::array_layer_tile_amount.size(); // no custom tiles found -> return max
}

constexpr uint8_t helper_size = (constants::array_layer_tile_amount.size() - custom_array_layer_index() + 1);

class GpuMultiArrayHelper {
public:
    struct Dictionary {
        nucleus::Raster<glm::u32vec2> packed_ids;
        nucleus::Raster<glm::u16vec2> layers;
    };
    struct MultiLayerInfo {
        tile::Id id;
        uint16_t index1;
        uint16_t index2;
    };
    GpuMultiArrayHelper();

    /// returns index in texture array
    glm::u16vec2 add_tile(const tile::Id& tile_id, uint8_t buffer_info);
    void remove_tile(const tile::Id& tile_id);
    void set_tile_limit(unsigned new_limit);
    constexpr static uint8_t buffer_amount() { return helper_size; };
    uint16_t layer_amount(uint8_t buffer_index) const;
    unsigned int n_occupied() const;
    Dictionary generate_dictionary() const;
    MultiLayerInfo layer(tile::Id tile_id) const;

private:
    uint8_t buffer_to_helper(uint8_t buffer) const;

    // array size is how many non -1u values are in the array_layer_quad_amount array + 1 for the default size
    // -> meaning the first gpu array helper is used for grid + -1u buffer
    // -> while the buffer with custom max layer amounts have individual gpuarrayhelpers
    std::array<nucleus::tile::GpuArrayHelper, helper_size> helpers;

}; // namespace nucleus::vector_layer

} // namespace nucleus::vector_layer
