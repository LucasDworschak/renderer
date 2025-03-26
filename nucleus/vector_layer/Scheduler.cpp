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

#include "Scheduler.h"
#include <nucleus/vector_tile/parse.h>

#include "nucleus/vector_layer/Preprocessor.h"

namespace nucleus::vector_layer {

Scheduler::Scheduler(const Scheduler::Settings& settings, Style&& style)
    : nucleus::tile::Scheduler(settings)
    , m_style(std::move(style))
    , m_needs_to_update_style(true)
{
    m_default_tile = nucleus::vector_layer::create_default_gpu_tile();
}

Scheduler::~Scheduler() = default;

void Scheduler::transform_and_emit(const std::vector<tile::DataQuad>& new_quads, const std::vector<tile::Id>& deleted_quads)
{
    std::vector<nucleus::tile::GpuVectorLayerTile> new_gpu_tiles;
    new_gpu_tiles.reserve(new_quads.size() * 4);

    for (const auto& quad : new_quads) {
        for (unsigned i = 0; i < 4; ++i) {

            GpuVectorLayerTile gpu_tile = nucleus::vector_layer::preprocess(quad.tiles[i].id, *quad.tiles[i].data, m_style);
            if (gpu_tile.id != quad.tiles[i].id) {
                gpu_tile = m_default_tile;
                gpu_tile.id = quad.tiles[i].id;
            }
            new_gpu_tiles.push_back(gpu_tile);
        }
    }

    emit gpu_tiles_updated(deleted_quads, new_gpu_tiles);
}

void Scheduler::set_enabled(bool new_enabled)
{
    nucleus::tile::Scheduler::set_enabled(new_enabled);

    if (m_needs_to_update_style) {
        emit style_updated(m_style.styles());
        m_needs_to_update_style = false;
    }
}

} // namespace nucleus::vector_layer
