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

#pragma once

#include "types.h"
#include <nucleus/Raster.h>

namespace nucleus::tile {

class GpuArrayHelper {
public:
    struct Dictionary {
        nucleus::Raster<glm::u32vec2> packed_ids;
        nucleus::Raster<uint16_t> array_layers;
    };
    struct LayerInfo {
        tile::Id id;
        uint16_t index;
    };

    GpuArrayHelper();

    /// returns index in texture array
    uint16_t add_tile(const tile::Id& tile_id, bool add_to_array = true, uint8_t additonal_info = 0, uint8_t additional_info_bitshift = 16);
    void remove_tile(const tile::Id& tile_id);
    bool contains_tile(const tile::Id& tile_id) const;
    void set_tile_limit(unsigned new_limit);
    unsigned size() const;
    unsigned int n_occupied() const;
    Dictionary generate_dictionary() const;
    LayerInfo layer(Id tile_id) const;

private:
    std::vector<tile::Id> m_array;
    tile::IdMap<uint16_t> m_id_to_array_layer;
};

} // namespace nucleus::tile
