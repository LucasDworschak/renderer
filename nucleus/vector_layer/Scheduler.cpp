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

Scheduler::Scheduler(QObject* parent)
    : nucleus::tile::Scheduler(256, parent)
    // , m_style(":/vectorlayerstyles/basemap.json")
    , m_style(":/vectorlayerstyles/openstreetmap.json")
// , m_style(":/vectorlayerstyles/qwant.json")
// , m_style(":/vectorlayerstyles/osm-bright.json")
{
    connect(&m_style, &Style::load_finished, this, &Scheduler::enable_scheduler);
}
Scheduler::~Scheduler() = default;

void Scheduler::transform_and_emit(const std::vector<tile::DataQuad>& new_quads, const std::vector<tile::Id>& deleted_quads)
{
    std::vector<nucleus::tile::GpuVectorLayerQuad> new_gpu_quads;
    new_gpu_quads.reserve(new_quads.size());

    // const auto style = m_style;

    std::transform(new_quads.cbegin(), new_quads.cend(), std::back_inserter(new_gpu_quads), [this](const auto& quad) {
        // create GpuQuad based on cpu quad
        GpuVectorLayerQuad gpu_quad;
        gpu_quad.id = quad.id;

        assert(quad.n_tiles == 4);
        for (unsigned i = 0; i < 4; ++i) {

            gpu_quad.tiles[i] = nucleus::vector_layer::preprocess(quad.tiles[i].id, *quad.tiles[i].data, m_style);
            gpu_quad.tiles[i].id = quad.tiles[i].id;
        }
        return gpu_quad;
    });

    emit gpu_quads_updated(new_gpu_quads, deleted_quads);
}

// TODO not quite sure if this is the correct way we want to do this yet..
// especially with changing vector layer source this would definitely cause problems
void Scheduler::load_style() { m_style.load(); }

void Scheduler::enable_scheduler(std::shared_ptr<const nucleus::Raster<glm::u32vec4>> styles)
{
    qDebug() << "vectorlayer style loaded";

    set_enabled(true);
    emit style_updated(styles);
}

bool Scheduler::is_ready_to_ship(const nucleus::tile::DataQuad& quad) const
{
    assert(m_geometry_ram_cache);
    return m_geometry_ram_cache->contains(quad.id);
}

void Scheduler::set_geometry_ram_cache(nucleus::tile::MemoryCache* new_geometry_ram_cache) { m_geometry_ram_cache = new_geometry_ram_cache; }

} // namespace nucleus::vector_layer
