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

#include <QDebug>

#include "GpuMultiArrayHelper.h"
#include <nucleus/srs.h>

// TODO: with the latest change of how the dictionary was created, we might be able to refactor this class a bit
// - since we no longer call generate_dictionary -> we can remove this function from here and from the individual helpers
// - since generate_dictionary is not used anymore we do not care anymore about hash conflicts -> we only need two helpers, one for acceleration structure and
// one for the data buffer (i think)

namespace nucleus::vector_layer {
GpuMultiArrayHelper::GpuMultiArrayHelper()
{
    // if this assert fails, it means we have to increase bits_for_buffer_info
    // in theory we should have 2*5 bits available with quad limit being 512
    assert(constants::array_layer_tile_amount.size() <= 4);
}

glm::u16vec2 GpuMultiArrayHelper::add_tile(const tile::Id& id, uint8_t buffer_info)
{
    constexpr auto buffer_offset = constants::array_helper_all_bits - constants::array_helper_buffer_info_bits;

    std::array<uint16_t, 3> layers;

    // active_helpers checks which buffer per tile should be fully used (other GpuArrayHelpers are still called but with false flag which is only used to have the same hashing collisions)
    std::array<bool, 3> active_helpers = { true, false, false };
    active_helpers[buffer_to_helper(buffer_info)] = true;

    for (uint8_t i = 0; i < buffer_amount(); i++) {
        layers[i] = helpers[i].add_tile(id, active_helpers[i], buffer_info, buffer_offset);
    }

    return { layers[0], layers[buffer_to_helper(buffer_info)] };
}

uint8_t GpuMultiArrayHelper::buffer_to_helper(uint8_t buffer) const
{
    if (constants::array_layer_tile_amount[buffer] == -1u)
        return 0;

    assert(buffer - custom_array_layer_index() + 1 > 0);

    return buffer - custom_array_layer_index() + 1;
}

void GpuMultiArrayHelper::remove_tile(const tile::Id& tile_id)
{
    for (auto& helper : helpers) {
        if (helper.contains_tile(tile_id))
            helper.remove_tile(tile_id);
    }
}

void GpuMultiArrayHelper::set_tile_limit(unsigned int new_limit)
{
    // if we have a limit greater than this assert, we cannot use the first bits to encode the buffer_info
    // -> either find a better way to encode buffer_info, or lower buffer amount.
    // Ultimately you can also refactor everything to send buffer info separately, increase layers to uint32_t or use rgb textures to encode the extra information in the b channel
    assert(new_limit <= (1 << (constants::array_helper_all_bits - constants::array_helper_buffer_info_bits))); // for 2 bits quad limit is 4096

    for (uint8_t i = 0; i < buffer_amount(); i++) {
        uint8_t offset = 0;
        if (i > 0)
            offset = custom_array_layer_index() - 1; // offset for buffer that use custom helpers

        // use the lower limit -> either limit of array_layer_quad_amount or the limit passed by the argument
        helpers[i].set_tile_limit((constants::array_layer_tile_amount[i + offset] < new_limit) ? constants::array_layer_tile_amount[i + offset] : new_limit);
    }
}

// constexpr uint8_t GpuMultiArrayHelper::buffer_amount() { return buffer_amount(); }
uint16_t GpuMultiArrayHelper::layer_amount(uint8_t buffer_index) const { return helpers[buffer_to_helper(buffer_index)].size(); }

unsigned GpuMultiArrayHelper::n_occupied() const
{
    // index is always filled -> while theoretically we also write sometimes to a second arrayHelper
    // it makes more sense to only provide the max amount of a helper that fills/deletes everytime add/remove tile is called.
    return unsigned(helpers[0].n_occupied());
}

GpuMultiArrayHelper::MultiLayerInfo GpuMultiArrayHelper::layer(tile::Id tile_id) const
{
    const auto layer_info_common = helpers[0].layer(tile_id);

    if (layer_info_common.id == tile::Id {})
        return { {}, 0, 0 }; // may be empty during startup.

    uint16_t index2 = layer_info_common.index;

    for (uint8_t i = 1; i < buffer_amount(); i++) { // we do not look int 0th buffer (since we use it as default)
        // go through every helper and execute layer with the id we know exists.
        const auto layer_info = helpers[i].layer(layer_info_common.id);

        if (layer_info.index != uint16_t(-1)) {
            assert(layer_info.id == layer_info_common.id);
            // we found the correct helper where we put our data
            index2 = layer_info.index;
            break;
        }
    }

    return { layer_info_common.id, layer_info_common.index, index2 };
}

GpuMultiArrayHelper::Dictionary GpuMultiArrayHelper::generate_dictionary() const
{
    GpuMultiArrayHelper::Dictionary dict;
    std::array<std::vector<uint16_t>, buffer_amount()> all_layers;
    glm::uvec2 layer_size;

    for (uint8_t i = 0; i < buffer_amount(); i++) {
        auto [packed_ids, layers] = helpers[i].generate_dictionary();

        all_layers[i] = std::move(layers.buffer());

        if (i == 0) {
            layer_size = layers.size();
            dict.packed_ids = std::move(packed_ids);
        }
    }

    dict.layers = { layer_size, glm::u16vec2(uint16_t(-1)) };
    std::vector<glm::u16vec2>& out_data = dict.layers.buffer();

    for (size_t i = 0; i < all_layers[0].size(); ++i) {
        out_data[i].x = all_layers[0][i]; // x value is always for the acceleration grid, which uses the first helper
        // y value uses the first layer that is not -1
        for (uint8_t j = buffer_amount(); j-- > 0;) { // backwards loop through all helpers
            if (all_layers[j][i] != uint16_t(-1)) {
                out_data[i].y = all_layers[j][i];

                break;
            }
        }
    }

    return dict;
}

} // namespace nucleus::vector_layer
