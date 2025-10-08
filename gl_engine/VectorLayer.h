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
    nucleus::tile::TileBounds bounds;
    uint layer;
};

struct FallbackMeta {
    bool finished;
    uint64_t first_check;
    uint64_t last_check;
    unsigned last_ortho_zoom;
};

class VectorLayer : public QObject {
    Q_OBJECT
public:
    explicit VectorLayer(unsigned fallback_resolution = 256, QObject* parent = nullptr);

    void init(ShaderRegistry* shader_registry); // needs OpenGL context
    void draw(const TileGeometry& tile_geometry,
        const nucleus::camera::Definition& camera,
        const std::vector<nucleus::tile::TileBounds>& draw_list) const;

    unsigned tile_count() const;

    void update_max_vector_geometry(unsigned int new_max_vector_geometry);
    void set_defines(const std::unordered_map<QString, QString>& defines);
    static std::unordered_map<QString, QString> default_defines();

    bool check_fallback_textures();

    void set_texture_layer(TextureLayer* texture_layer);

public slots:
    void update_gpu_tiles(const std::vector<nucleus::tile::Id>& deleted_tiles, const std::vector<nucleus::tile::GpuVectorLayerTile>& new_tiles);
    void set_tile_limit(unsigned new_limit);
    void update_style(std::shared_ptr<const nucleus::Raster<glm::u32vec2>> styles);

private:
    bool m_initialized;

    const unsigned m_fallback_resolution;
    int m_max_vector_geometry = 8;

    unsigned m_fallback_framebuffer = unsigned(-1);

    TextureLayer* m_texture_layer;

    std::shared_ptr<ShaderProgram> m_shader;
    std::shared_ptr<ShaderProgram> m_fallback_shader;

    std::unique_ptr<Texture> m_acceleration_grid_texture;
    std::vector<std::unique_ptr<Texture>> m_geometry_buffer_texture;
    std::unique_ptr<Texture> m_styles_texture;

    std::unique_ptr<Texture> m_fallback_texture_array_higher;
    std::unique_ptr<Texture> m_fallback_texture_array_lower;

    std::unique_ptr<Texture> m_instanced_zoom;
    std::unique_ptr<Texture> m_instanced_array_index;
    std::unique_ptr<Texture> m_instanced_zoom_fallback;
    std::unique_ptr<Texture> m_instanced_array_index_fallback;

    nucleus::vector_layer::GpuMultiArrayHelper m_gpu_multi_array_helper;

    std::shared_ptr<const nucleus::Raster<glm::u32vec2>> m_styles;

    std::unordered_map<QString, QString> m_defines;

    helpers::ScreenQuadGeometry m_screen_quad_geometry;

    // TODO i think we can move tile id into fallbackmeta
    std::vector<nucleus::tile::Id> m_vector_on_gpu;
    std::vector<FallbackMeta> m_fallback_on_gpu;

    static constexpr int max_fallback_renders_per_frame = 32;
    static constexpr int ms_between_fallback_checks = 500;
    static constexpr int max_ms_for_ortho_wait = 20000;
    static constexpr int max_ortho_zoom_level = 17;

    bool m_fallback_render_possible;
    nucleus::tile::IdMap<uint16_t> m_gpu_fallback_map;

    nucleus::tile::GpuArrayHelper::LayerInfo fallback_layer(nucleus::tile::Id tile_id) const;
    void update_fallback_textures(const std::vector<IdLayer>& tiles_to_render);
    void setup_buffers(std::shared_ptr<ShaderProgram> shader, const std::vector<nucleus::tile::TileBounds>& draw_list, bool draw_fallback) const;
    static std::vector<glm::vec2> gaussian2d(float sigma, int num_points);
    static QString vec2_arr_to_glsl_string(const std::vector<glm::vec2>& points);
};
} // namespace gl_engine
