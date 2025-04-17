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

#include <QObject>
#include <nucleus/Raster.h>
#include <nucleus/tile/DrawListGenerator.h>
#include <nucleus/tile/types.h>
#include <nucleus/vector_layer/GpuMultiArrayHelper.h>

namespace camera {
class Definition;
}

class QOpenGLShaderProgram;
class QOpenGLBuffer;
class QOpenGLVertexArrayObject;

namespace gl_engine {
class ShaderRegistry;
class ShaderProgram;
class Texture;
class TileGeometry;

class VectorLayer : public QObject {
    Q_OBJECT
public:
    explicit VectorLayer(QObject* parent = nullptr);
    void init(ShaderRegistry* shader_registry); // needs OpenGL context
    void draw(const TileGeometry& tile_geometry, const nucleus::camera::Definition& camera, const std::vector<nucleus::tile::TileBounds>& draw_list) const;

    unsigned tile_count() const;

public slots:
    void update_gpu_tiles(const std::vector<nucleus::tile::Id>& deleted_tiles, const std::vector<nucleus::tile::GpuVectorLayerTile>& new_tiles);
    void set_tile_limit(unsigned new_limit);
    void update_style(std::shared_ptr<const nucleus::Raster<glm::u32vec4>> styles);

private:
    std::shared_ptr<ShaderProgram> m_shader;

    std::unique_ptr<Texture> m_acceleration_grid_texture;
    std::vector<std::unique_ptr<Texture>> m_geometry_buffer_texture;
    std::unique_ptr<Texture> m_styles_texture;

    std::unique_ptr<Texture> m_instanced_zoom;
    std::unique_ptr<Texture> m_instanced_array_index;

    nucleus::vector_layer::GpuMultiArrayHelper m_gpu_multi_array_helper;

    std::shared_ptr<const nucleus::Raster<glm::u32vec4>> m_styles;

    bool m_initialized;

    // nucleus::tile::IdMap<std::shared_ptr<const std::vector<uint32_t>>> m_id_to_data_bridge;
    // nucleus::tile::IdMap<std::shared_ptr<const std::vector<uint32_t>>> m_id_to_triangle_data;
};
} // namespace gl_engine
