/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2024 Adam Celarek
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

#include "GpuArrayHelper.h"
#include <nucleus/srs.h>

namespace nucleus::tile {
GpuArrayHelper::GpuArrayHelper() { }

/**
 * adds tile id to m_id_to_array_layer
 * add_to_array can be set to false to prevent the id actually be added to the array
 * -> this is useful to allow multiple array helpers to be used (with less than max layer amount) and have the hashing collisions be the same for all helpers
 * if additional info is set; the value will be orred with the id_to_array_layer map value
 */
uint16_t GpuArrayHelper::add_tile(const tile::Id& id, bool add_to_array, uint8_t additonal_info, uint8_t additional_info_bitshift)
{
    assert(!m_id_to_array_layer.contains(id));

    uint16_t array_layer = uint16_t(-1u);

    if (add_to_array) {
        const auto t = std::find(m_array.begin(), m_array.end(), tile::Id { unsigned(-1), {} });
        assert(t != m_array.end());
        *t = id;

        // returns index in texture array
        array_layer = uint16_t(t - m_array.begin());
    }

    assert(m_array.size() - 1 <= (1u << additional_info_bitshift)); // make sure that the bitshifted additional data does not overlap array_layer data

    m_id_to_array_layer.emplace(id, array_layer | (additonal_info << additional_info_bitshift));

    return array_layer;
}

void GpuArrayHelper::remove_tile(const tile::Id& tile_id)
{
    assert(m_id_to_array_layer.contains(tile_id));
    const auto layer = m_id_to_array_layer.at(tile_id);
    m_id_to_array_layer.erase(tile_id);

    // only try to remove it from the array if the layer is actually set
    if (layer != uint16_t(-1)) {
        const auto t = std::find(m_array.begin(), m_array.end(), tile_id);
        assert(t != m_array.end()); // removing a tile that's not here. likely there is a race.
        *t = tile::Id { unsigned(-1), {} };
    }
}

bool GpuArrayHelper::contains_tile(const tile::Id& tile_id) { return m_id_to_array_layer.contains(tile_id); }

void GpuArrayHelper::set_quad_limit(unsigned int new_limit)
{
    assert(m_array.empty());
    m_array.resize(new_limit * 4);
    std::fill(m_array.begin(), m_array.end(), tile::Id { unsigned(-1), {} });
}

unsigned GpuArrayHelper::size() const { return unsigned(m_array.size()); }

GpuArrayHelper::Dictionary GpuArrayHelper::generate_dictionary() const
{
    const auto hash_to_pixel = [](uint16_t hash) { return glm::uvec2(hash & 255, hash >> 8); };
    nucleus::Raster<glm::u32vec2> packed_ids({ 256, 256 }, glm::u32vec2(-1, -1));
    nucleus::Raster<uint16_t> layers({ 256, 256 }, 0);
    for (const auto& [id, layer] : m_id_to_array_layer) {
        auto hash = nucleus::srs::hash_uint16(id);
        while (packed_ids.pixel(hash_to_pixel(hash)) != glm::u32vec2(-1, -1))
            hash++;

        packed_ids.pixel(hash_to_pixel(hash)) = nucleus::srs::pack(id);
        layers.pixel(hash_to_pixel(hash)) = layer;
    }

    return { packed_ids, layers };
}
} // namespace nucleus::tile
