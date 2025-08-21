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

#include "VectorLayer.h"

#include "ShaderProgram.h"
#include "ShaderRegistry.h"
#include "Texture.h"
#include "TextureLayer.h"
#include "TileGeometry.h"
#include <QOpenGLExtraFunctions>

#include <gl_engine/Framebuffer.h>

#include "nucleus/vector_layer/constants.h"

namespace gl_engine {

using namespace nucleus::vector_layer;

VectorLayer::VectorLayer(unsigned int fallback_resolution, QObject* parent)
    : QObject { parent }
    , m_initialized(false)
    , m_fallback_resolution(fallback_resolution)
    , m_gpu_multi_array_helper()

{
    m_defines = default_defines();
}

std::unordered_map<QString, QString> gl_engine::VectorLayer::default_defines()
{
    std::unordered_map<QString, QString> defines;

    defines[QString("style_bits")] = QString::number(constants::style_bits);
    defines[QString("style_precision")] = QString::number(constants::style_precision);
    defines[QString("max_zoom")] = QString::number(constants::style_zoom_range.y);
    defines[QString("tile_extent")] = QString("float(%1)").arg(constants::tile_extent);
    defines[QString("scale_polygons")] = QString("float(%1)").arg(constants::scale_polygons);
    defines[QString("scale_lines")] = QString("float(%1)").arg(constants::scale_lines);
    defines[QString("grid_size")] = QString("vec2(%1,%1)").arg(constants::grid_size);
    defines[QString("mipmap_levels")] = QString::number(constants::mipmap_levels);

    defines[QString("all_bits")] = QString::number(constants::all_bits);
    defines[QString("coordinate_bits_polygons")] = QString::number(constants::coordinate_bits_polygons);
    defines[QString("coordinate_bits_lines")] = QString::number(constants::coordinate_bits_lines);
    defines[QString("aa_border")] = QString::number(constants::aa_border);

    defines[QString("sampler_offset")] = QString::number(constants::array_helper_all_bits - constants::array_helper_buffer_info_bits);

    return defines;
}

void gl_engine::VectorLayer::set_defines(const std::unordered_map<QString, QString>& defines) { m_defines = defines; }

void gl_engine::VectorLayer::init(ShaderRegistry* shader_registry)
{
    m_screen_quad_geometry = gl_engine::helpers::create_screen_quad_geometry();

    std::vector<QString> defines;
    for (const auto& define : m_defines) {
        defines.push_back("#define " + define.first + " " + define.second);
    }

    std::vector<QString> defines_fallback(defines);
    // defines_fallback.push_back("#define FALLBACK_MODE 1");
    // defines_fallback.push_back("#define SHALLOW_ANGLE_SIGNAL_FREQUENCY 2"); // no angle smoothing
    // defines_fallback.push_back("#define VIEW_MODE 2"); // vector only (for now)

    m_shader = std::make_shared<ShaderProgram>("tile.vert", "vector_layer.frag", ShaderCodeSource::FILE, defines);
    m_fallback_shader = std::make_shared<ShaderProgram>("vector_layer_fallback.vert", "vector_layer_fallback.frag", ShaderCodeSource::FILE, defines_fallback);
    shader_registry->add_shader(m_shader);
    shader_registry->add_shader(m_fallback_shader);

    m_acceleration_grid_texture = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::R32UI);
    m_acceleration_grid_texture->setParams(gl_engine::Texture::Filter::Nearest, gl_engine::Texture::Filter::Nearest);
    // // TODO: might become larger than GL_MAX_ARRAY_TEXTURE_LAYERS
    m_acceleration_grid_texture->allocate_array(constants::grid_size, constants::grid_size, unsigned(m_gpu_multi_array_helper.layer_amount(0)));

    m_instanced_zoom = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::R8UI);
    m_instanced_zoom->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_instanced_array_index = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::RG16UI);
    m_instanced_array_index->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_geometry_buffer_texture.resize(constants::array_layer_tile_amount.size());
    for (uint8_t i = 0; i < constants::array_layer_tile_amount.size(); i++) {
        const auto layer_amount = m_gpu_multi_array_helper.layer_amount(i);

        m_geometry_buffer_texture[i] = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::RG32UI);
        m_geometry_buffer_texture[i]->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);
        m_geometry_buffer_texture[i]->allocate_array(constants::data_size[i], constants::data_size[i], layer_amount);
    }

    m_styles_texture = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::RGBA32UI);
    m_styles_texture->setParams(gl_engine::Texture::Filter::Nearest, gl_engine::Texture::Filter::Nearest);

    m_fallback_texture_array = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::RGBA8);
    m_fallback_texture_array->setParams(Texture::Filter::MipMapLinear, Texture::Filter::Linear, true);
    m_fallback_texture_array->allocate_array(m_fallback_resolution, m_fallback_resolution, unsigned(m_gpu_multi_array_helper.layer_amount(0)));

    m_initialized = true;
    if (m_styles) {
        m_styles_texture->upload(*m_styles);
    }
}

/**
 * returns true if another update step is required (-> we need to call window.update_requested())
 */
bool VectorLayer::check_fallback_textures()
{
    if (!m_fallback_render_possible)
        return false; // there shouldn't be any fallback textures to render -> early exit

    IdLayer fallbacks_to_render;

    for (uint i = 0; i < m_vector_on_gpu.size(); i++) {
        if (m_vector_on_gpu[i] != m_fallback_on_gpu[i]) {
            // TODO also do ortho fotos
            fallbacks_to_render.layer.push_back(i);
            fallbacks_to_render.bounds.push_back({ m_vector_on_gpu[i], {} });

            if (fallbacks_to_render.layer.size() >= max_fallback_renders_per_frame)
                break; // early exit if we are full
        }
    }

    // if we found max amount to render -> we need to look at the next frame to determine if we rendered all
    bool search_next_frame = fallbacks_to_render.layer.size() >= max_fallback_renders_per_frame;
    if (!search_next_frame)
        m_fallback_render_possible = false; // we found all -> we can early exit until we get new tiles

    update_fallback_textures(fallbacks_to_render);

    return search_next_frame;
}

unsigned VectorLayer::tile_count() const { return m_gpu_multi_array_helper.n_occupied(); }

void VectorLayer::setup_buffers(std::shared_ptr<ShaderProgram> shader, const std::vector<nucleus::tile::TileBounds>& draw_list) const
{
    shader->set_uniform("acceleration_grid_sampler", 9);
    m_acceleration_grid_texture->bind(9);

    nucleus::Raster<uint8_t> zoom_level_raster = { glm::uvec2 { 1024, 1 } };
    nucleus::Raster<glm::u16vec2> array_index_raster = { glm::uvec2 { 1024, 4 } };
    for (unsigned i = 0; i < std::min(unsigned(draw_list.size()), 1024u); ++i) {
        // todo -> automatically go to parent ( or parent.parent())
        const auto layer0 = m_gpu_multi_array_helper.layer(draw_list[i].id);
        zoom_level_raster.pixel({ i, 0 }) = layer0.id.zoom_level;
        array_index_raster.pixel({ i, 0 }) = { layer0.index1, layer0.index2 };

        // load the mipmaps -> also make sure that all mipmaps are loaded -> otherwise use parent mipmap
        const auto layer1 = m_gpu_multi_array_helper.layer(layer0.id.parent());
        if (layer1.id.zoom_level < 200) {
            array_index_raster.pixel({ i, 1 }) = { layer1.index1, layer1.index2 };
        } else {
            array_index_raster.pixel({ i, 1 }) = array_index_raster.pixel({ i, 0 });
        }
        const auto layer2 = m_gpu_multi_array_helper.layer(layer0.id.parent().parent());
        if (layer2.id.zoom_level < 200) {
            array_index_raster.pixel({ i, 2 }) = { layer2.index1, layer2.index2 };
        } else {
            array_index_raster.pixel({ i, 2 }) = array_index_raster.pixel({ i, 1 });
        }
        const auto layer3 = m_gpu_multi_array_helper.layer(layer0.id.parent().parent().parent());
        if (layer3.id.zoom_level < 200) {
            array_index_raster.pixel({ i, 3 }) = { layer3.index1, layer3.index2 };
        } else {
            array_index_raster.pixel({ i, 3 }) = array_index_raster.pixel({ i, 2 });
        }
    }

    m_instanced_array_index->bind(10);
    shader->set_uniform("instanced_texture_array_index_sampler_vector", 10);
    m_instanced_array_index->upload(array_index_raster);

    m_instanced_zoom->bind(11);
    shader->set_uniform("instanced_texture_zoom_sampler_vector", 11);
    m_instanced_zoom->upload(zoom_level_raster);

    shader->set_uniform("styles_sampler", 12);
    m_styles_texture->bind(12);

    // upload all geometry buffers
    // binds the buffers to "...buffer_sampler_[0-max]"
    // constexpr uint8_t triangle_vertex_buffer_start = 13;
    // for (uint8_t i = 0; i < constants::array_layer_tile_amount.size(); i++) {
    //     shader->set_uniform("geometry_buffer_sampler_" + std::to_string(i), triangle_vertex_buffer_start + i);
    //     m_geometry_buffer_texture[i]->bind(triangle_vertex_buffer_start + i);
    // }

    shader->set_uniform("geometry_buffer_sampler_0", 13);
    m_geometry_buffer_texture[0]->bind(13);
    shader->set_uniform("geometry_buffer_sampler_1", 14);
    m_geometry_buffer_texture[1]->bind(14);
    shader->set_uniform("geometry_buffer_sampler_2", 15);
    m_geometry_buffer_texture[2]->bind(15);
    shader->set_uniform("geometry_buffer_sampler_3", 16);
    m_geometry_buffer_texture[3]->bind(16);
}

void VectorLayer::draw(const TileGeometry& tile_geometry,
    const TextureLayer& texture_layer,
    const nucleus::camera::Definition& camera,
    const std::vector<nucleus::tile::TileBounds>& draw_list) const
{
    m_shader->bind();

    texture_layer.bind_buffer(m_shader, draw_list); // binds texture_sampler -> for ortho fotos

    m_shader->set_uniform("max_vector_geometry", m_max_vector_geometry);

    setup_buffers(m_shader, draw_list);

    // constexpr uint8_t next_buffer_start = 13 + constants::array_layer_tile_amount.size();
    constexpr uint8_t next_buffer_start = 17;
    m_shader->set_uniform("fallback_texture_array", next_buffer_start);
    m_fallback_texture_array->bind(next_buffer_start);

    tile_geometry.draw(m_shader.get(), camera, draw_list);
}

void VectorLayer::update_gpu_tiles(const std::vector<nucleus::tile::Id>& deleted_tiles, const std::vector<nucleus::tile::GpuVectorLayerTile>& new_tiles)
{
    if (!QOpenGLContext::currentContext()) // can happen during shutdown.
        return;

    for (const auto& quad : deleted_tiles) {
        for (const auto& id : quad.children()) {
            const auto exact_layer = m_gpu_multi_array_helper.exact_layer(0, id);
            if (exact_layer != -1u) {
                m_vector_on_gpu[exact_layer] = {};
                m_fallback_on_gpu[exact_layer] = {};
            }

            m_gpu_multi_array_helper.remove_tile(id);
        }
    }

    for (const auto& tile : new_tiles) {
        // test for validity
        assert(tile.id.zoom_level < 100);
        // if (!tile.acceleration_grid)
        //     continue; // nothing here

        assert(tile.acceleration_grid);
        assert(tile.geometry_buffer);

        // find empty spot and upload texture
        const auto layer_index = m_gpu_multi_array_helper.add_tile(tile.id, tile.buffer_info);

        m_acceleration_grid_texture->upload(*tile.acceleration_grid, layer_index[0]);
        m_geometry_buffer_texture[tile.buffer_info]->upload(*tile.geometry_buffer, layer_index[1]);

        m_vector_on_gpu[layer_index[0]] = tile.id;

        m_fallback_render_possible = true; // there was at least one new tile -> check if we need to render fallbacks
    }
}

void VectorLayer::update_fallback_textures(const IdLayer& id_layer)
{
    // SETUP framebuffer
    QOpenGLExtraFunctions* f = QOpenGLContext::currentContext()->extraFunctions();

    unsigned m_frame_buffer = unsigned(-1);
    f->glGenFramebuffers(1, &m_frame_buffer);
    f->glViewport(0, 0, m_fallback_resolution, m_fallback_resolution);

    f->glDisable(GL_CULL_FACE);
    f->glDisable(GL_DEPTH_TEST);

    f->glBindFramebuffer(GL_FRAMEBUFFER, m_frame_buffer);

    m_fallback_shader->bind();

    m_fallback_shader->set_uniform("min_vector_geometry", m_max_vector_geometry);

    setup_buffers(m_fallback_shader, id_layer.bounds);

    for (unsigned i = 0; i < id_layer.layer.size(); i++) {
        m_fallback_texture_array->bind_layer_to_frame_buffer(0, id_layer.layer[i]);

        const GLfloat clearAlbedoColor[] = { 0.0f, 0.0f, 0.0f, 1.0f };
        f->glClearBufferfv(GL_COLOR, 0, clearAlbedoColor);

        m_fallback_shader->set_uniform("instance_id", int(i));
        m_fallback_shader->set_uniform("tile_zoom", int(id_layer.bounds[i].id.zoom_level));

        // it is possible to render 4 or 8 tiles at once (using different color attachments) -> but not sure if really performance relevant
        m_screen_quad_geometry.draw();
        m_fallback_on_gpu[id_layer.layer[i]] = id_layer.bounds[i].id;
    }

    // generate mipmaps after all the individual layers have been rendered
    m_fallback_texture_array->generate_mipmap();

    f->glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void VectorLayer::set_tile_limit(unsigned int new_limit)
{
    assert(new_limit <= 2048); // array textures with size > 2048 are not supported on all devices
    assert(!m_acceleration_grid_texture);
    assert(m_geometry_buffer_texture.size() == 0);

    m_gpu_multi_array_helper.set_tile_limit(new_limit);

    m_vector_on_gpu.resize(new_limit, {});
    m_fallback_on_gpu.resize(new_limit, {});
}

void VectorLayer::update_style(std::shared_ptr<const nucleus::Raster<glm::u32vec4>> styles)
{
    // at this point we cannot be sure that the texture has been initialized -> save the pointer for now and upload after init
    m_styles = styles;

    if (m_initialized) // only upload if it is already initialized
    {
        m_styles_texture->upload(*m_styles);
    }
}

void VectorLayer::update_max_vector_geometry(unsigned int new_max_vector_geometry)
{
    m_max_vector_geometry = int(new_max_vector_geometry); // uint uniforms do not work ...;
}

} // namespace gl_engine
