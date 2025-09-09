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
#include "nucleus/vector_layer/Preprocessor.h"

namespace nucleus::vector_layer {

class Scheduler : public nucleus::tile::Scheduler {
    Q_OBJECT
public:
    explicit Scheduler(const Scheduler::Settings& settings, Style&& style);
    ~Scheduler() override;

    void set_enabled(bool new_enabled) override;

signals:
    void gpu_tiles_updated(const std::vector<tile::Id>& deleted_tiles, const std::vector<nucleus::tile::GpuVectorLayerTile>& new_tiles);

    void style_updated(std::shared_ptr<const nucleus::Raster<glm::u32vec2>> style_buffer);

protected:
    void transform_and_emit(const std::vector<tile::DataQuad>& new_quads, const std::vector<tile::Id>& deleted_quads) override;

private:
    bool m_needs_to_update_style;
    nucleus::tile::GpuVectorLayerTile m_default_tile;

    Preprocessor m_preprocessor;
};

} // namespace nucleus::vector_layer
