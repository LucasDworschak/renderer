/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2024 Adam Celarek
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

#include <nucleus/tile/Scheduler.h>
#include <nucleus/vector_tile/types.h>

#include "Style.h"

namespace nucleus::vector_layer {

class Scheduler : public nucleus::tile::Scheduler {
    Q_OBJECT
public:
    explicit Scheduler(std::string name, QObject* parent = nullptr);
    ~Scheduler() override;

    void set_geometry_ram_cache(nucleus::tile::MemoryCache* new_geometry_ram_cache);

    void load_style();

public slots:
    void enable_scheduler(std::shared_ptr<const nucleus::Raster<glm::u32vec4>>, std::shared_ptr<const nucleus::Raster<glm::u32vec4>>);

signals:
    void gpu_quads_updated(const std::vector<nucleus::tile::GpuVectorLayerQuad>& new_quads, const std::vector<tile::Id>& deleted_quads);
    void style_updated(std::shared_ptr<const nucleus::Raster<glm::u32vec4>> fill_styles, std::shared_ptr<const nucleus::Raster<glm::u32vec4>> line_styles);

protected:
    void transform_and_emit(const std::vector<tile::DataQuad>& new_quads, const std::vector<tile::Id>& deleted_quads) override;
    bool is_ready_to_ship(const nucleus::tile::DataQuad& quad) const override;

private:
    nucleus::tile::MemoryCache* m_geometry_ram_cache = nullptr;

    Style m_style;
};

} // namespace nucleus::vector_layer
