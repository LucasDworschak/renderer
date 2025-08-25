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

#include "helpers.h"
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
class TextureLayer;

struct IdLayer {
    std::vector<nucleus::tile::TileBounds> bounds;
    std::vector<uint> layer;
};

class VectorLayer : public QObject {
    Q_OBJECT
public:
    explicit VectorLayer(unsigned fallback_resolution = 256, QObject* parent = nullptr);

    void init(ShaderRegistry* shader_registry); // needs OpenGL context
    void draw(const TileGeometry& tile_geometry,
        const TextureLayer& texture_layer,
        const nucleus::camera::Definition& camera,
        const std::vector<nucleus::tile::TileBounds>& draw_list) const;

    unsigned tile_count() const;

    void update_max_vector_geometry(unsigned int new_max_vector_geometry);
    void set_defines(const std::unordered_map<QString, QString>& defines);
    static std::unordered_map<QString, QString> default_defines();

    bool check_fallback_textures();
signals:
    void fallback_textures_rendered();
public slots:
    void update_gpu_tiles(const std::vector<nucleus::tile::Id>& deleted_tiles, const std::vector<nucleus::tile::GpuVectorLayerTile>& new_tiles);
    void set_tile_limit(unsigned new_limit);
    void update_style(std::shared_ptr<const nucleus::Raster<glm::u32vec4>> styles);
    void generate_fallback_mipmaps();

private:
    bool m_initialized;

    const unsigned m_fallback_resolution;
    int m_max_vector_geometry = 8;

    unsigned m_fallback_framebuffer = unsigned(-1);

    std::shared_ptr<ShaderProgram> m_shader;
    std::shared_ptr<ShaderProgram> m_fallback_shader;

    std::unique_ptr<Texture> m_acceleration_grid_texture;
    std::vector<std::unique_ptr<Texture>> m_geometry_buffer_texture;
    std::unique_ptr<Texture> m_styles_texture;

    std::unique_ptr<Texture> m_fallback_texture_array;

    std::unique_ptr<Texture> m_instanced_zoom;
    std::unique_ptr<Texture> m_instanced_array_index;

    nucleus::vector_layer::GpuMultiArrayHelper m_gpu_multi_array_helper;

    std::shared_ptr<const nucleus::Raster<glm::u32vec4>> m_styles;

    std::unordered_map<QString, QString> m_defines;

    helpers::ScreenQuadGeometry m_screen_quad_geometry;

    std::vector<nucleus::tile::Id> m_vector_on_gpu;
    std::vector<nucleus::tile::Id> m_fallback_on_gpu;

    static constexpr int max_fallback_renders_per_frame = 32;

    bool m_fallback_render_possible;

    void update_fallback_textures(const IdLayer& fallbacks_to_render);
    void setup_buffers(std::shared_ptr<ShaderProgram> shader, const std::vector<nucleus::tile::TileBounds>& draw_list) const;
};
} // namespace gl_engine
