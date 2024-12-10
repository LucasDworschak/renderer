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

Scheduler::Scheduler(std::string name, QObject* parent)
    : nucleus::tile::Scheduler(std::move(name), 256, parent)
{
}
Scheduler::~Scheduler() = default;

void Scheduler::transform_and_emit(const std::vector<tile::DataQuad>& new_quads, const std::vector<tile::Id>& deleted_quads)
{
    std::vector<nucleus::tile::GpuVectorLayerQuad> new_gpu_quads;
    new_gpu_quads.reserve(new_quads.size());

    std::transform(new_quads.cbegin(), new_quads.cend(), std::back_inserter(new_gpu_quads), [](const auto& quad) {
        // create GpuQuad based on cpu quad
        GpuVectorLayerQuad gpu_quad;
        gpu_quad.id = quad.id;
        assert(quad.n_tiles == 4);
        for (unsigned i = 0; i < 4; ++i) {
            gpu_quad.tiles[i] = nucleus::vector_layer::preprocess(*quad.tiles[i].data);
            gpu_quad.tiles[i].id = quad.tiles[i].id;
        }
        return gpu_quad;
    });

    emit gpu_quads_updated(new_gpu_quads, deleted_quads);
}

bool Scheduler::is_ready_to_ship(const nucleus::tile::DataQuad& quad) const
{
    assert(m_geometry_ram_cache);
    return m_geometry_ram_cache->contains(quad.id);
}

void Scheduler::set_geometry_ram_cache(nucleus::tile::MemoryCache* new_geometry_ram_cache) { m_geometry_ram_cache = new_geometry_ram_cache; }

} // namespace nucleus::vector_layer
