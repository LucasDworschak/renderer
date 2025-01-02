/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2024 Adam Celarek
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

#include "types.h"
#include <nucleus/Raster.h>

namespace nucleus::tile {

class GpuArrayHelper {
public:
    struct Dictionary {
        nucleus::Raster<glm::u32vec2> packed_ids;
        nucleus::Raster<uint16_t> layers;
    };
    struct TileIDWithSubLayer {
        tile::Id id;
        uint8_t sub_layer;

        operator std::tuple<unsigned, unsigned, unsigned, unsigned, unsigned>() const { return std::make_tuple(id.zoom_level, id.coords.x, id.coords.y, unsigned(id.scheme), sub_layer); }

        bool operator==(const TileIDWithSubLayer& other) const = default;
        using Hasher = typename radix::hasher::for_tuple<unsigned, unsigned, unsigned, unsigned, unsigned>;
    };

    template <typename T> using IdMapSub = std::unordered_map<TileIDWithSubLayer, T, TileIDWithSubLayer::Hasher>;

    GpuArrayHelper();

    /// returns index in texture array
    unsigned add_tile(const tile::Id& tile_id, const uint8_t sub_layer = 0);
    void remove_tile(const tile::Id& tile_id, const uint8_t sub_layer = 0);
    bool contains_tile(const tile::Id& tile_id, const uint8_t sub_layer = 0);
    void set_quad_limit(unsigned new_limit);
    unsigned size() const;
    Dictionary generate_dictionary() const;

private:
    std::vector<TileIDWithSubLayer> m_array;
    IdMapSub<unsigned> m_id_to_layer;
};

} // namespace nucleus::tile
