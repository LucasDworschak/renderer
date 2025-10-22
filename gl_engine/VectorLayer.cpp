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
    defines[QString("max_offset_levels")] = QString::number(constants::max_offset_levels);
    defines[QString("aa_sample_dist")] = QString("float(%1)").arg(constants::aa_sample_dist);

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
    {
        // fallback meta structure
        QOpenGLExtraFunctions* f = QOpenGLContext::currentContext()->extraFunctions();

        f->glGenFramebuffers(1, &m_fallback_framebuffer);

        m_screen_quad_geometry = gl_engine::helpers::create_screen_quad_geometry();
    }

    std::vector<QString> defines;
    for (const auto& define : m_defines) {
        defines.push_back("#define " + define.first + " " + define.second);
    }

    std::vector<QString> defines_fallback(defines);
    // defines_fallback.push_back("#define FALLBACK_MODE 1");
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

    m_instanced_zoom_fallback = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::R8UI);
    m_instanced_zoom_fallback->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_instanced_array_index_fallback = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::R16UI);
    m_instanced_array_index_fallback->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_geometry_buffer_texture.resize(constants::array_layer_tile_amount.size());
    for (unsigned i = 0; i < constants::array_layer_tile_amount.size(); i++) {
        const auto layer_amount = m_gpu_multi_array_helper.layer_amount(i);

        m_geometry_buffer_texture[i] = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::RG32UI);
        m_geometry_buffer_texture[i]->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);
        m_geometry_buffer_texture[i]->allocate_array(constants::data_size[i], constants::data_size[i], layer_amount);
    }

    m_styles_texture = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::RG32UI);
    m_styles_texture->setParams(gl_engine::Texture::Filter::Nearest, gl_engine::Texture::Filter::Nearest);

    m_fallback_texture_array_higher = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::RGB565);
    m_fallback_texture_array_higher->setParams(Texture::Filter::MipMapLinear, Texture::Filter::Linear, true);
    m_fallback_texture_array_higher->allocate_array(m_fallback_resolution, m_fallback_resolution, unsigned(m_gpu_multi_array_helper.layer_amount(0)));

    m_fallback_texture_array_lower = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::RGB565);
    m_fallback_texture_array_lower->setParams(Texture::Filter::MipMapLinear, Texture::Filter::Linear, true);
    m_fallback_texture_array_lower->allocate_array(
        m_fallback_resolution / 2.0, m_fallback_resolution / 2.0, unsigned(m_gpu_multi_array_helper.layer_amount(0)));

    m_initialized = true;
    if (m_styles) {
        m_styles_texture->upload(*m_styles);
    }
}

void VectorLayer::set_texture_layer(TextureLayer* texture_layer) { m_texture_layer = texture_layer; }

/**
 * returns true to indicate that window.update_requested() should be called
 */
bool VectorLayer::check_fallback_textures()
{
    if (!m_fallback_render_possible)
        return false; // there shouldn't be any fallback textures to render -> early exit

    std::vector<IdLayer> tiles_to_render;

    assert(m_texture_layer);

    unsigned num_unfinished = 0;

    const auto current_time = nucleus::utils::time_since_epoch();

    for (uint i = 0; i < m_vector_on_gpu.size(); i++) {
        // early exit if tile is finished or not even vector data on gpu
        if (m_fallback_on_gpu[i].finished || m_vector_on_gpu[i] == nucleus::tile::Id {})
            continue;

        num_unfinished++;

        // try to render if not finished and diff between last check and now is over threshold
        if (current_time > m_fallback_on_gpu[i].last_check + ms_between_fallback_checks) {

            // ortho uses quads and does not have valid data over certain zoom level
            auto ortho_id = m_vector_on_gpu[i];
            while (ortho_id.zoom_level > max_ortho_zoom_level)
                ortho_id = ortho_id.parent();
            if (ortho_id.zoom_level > 0)
                ortho_id = ortho_id.parent();

            // only finished if ortho tile is loaded
            // additionally we also check if enough time between first_check and now has passed (to allow a worse ortho tile if it cannot be loaded)
            const bool tile_finished = m_texture_layer->is_tile_loaded(ortho_id) || (current_time > m_fallback_on_gpu[i].first_check + max_ms_for_ortho_wait);
            m_fallback_on_gpu[i].finished = tile_finished;
            m_fallback_on_gpu[i].last_check = current_time;

            // only try to render if ortho zoom_level improved from last render
            const auto ortho_zoom_level = m_texture_layer->tile(ortho_id).id.zoom_level;
            if (ortho_zoom_level <= m_fallback_on_gpu[i].last_ortho_zoom)
                continue;
            m_fallback_on_gpu[i].last_ortho_zoom = ortho_zoom_level;

            tiles_to_render.push_back({ { m_vector_on_gpu[i], {} }, i });

            if (tiles_to_render.size() >= max_fallback_renders_per_frame)
                break; // early exit if we are full
        }
    }

    update_fallback_textures(tiles_to_render);

    // we are finished with this method if fallbacks_to_render is not max amount and there aren't any orthos pending
    if (tiles_to_render.size() < max_fallback_renders_per_frame && num_unfinished == tiles_to_render.size()) {
        m_fallback_render_possible = false; // we found all -> we can early exit until we get new tiles
    }

    return true; // always return true here to update
}

unsigned VectorLayer::tile_count() const { return m_gpu_multi_array_helper.n_occupied(); }

nucleus::tile::GpuArrayHelper::LayerInfo VectorLayer::fallback_layer(nucleus::tile::Id tile_id) const
{
    if (tile_id.zoom_level > 200) // zoom may be -1u during startup
        return { {}, 0 };

    while (!m_gpu_fallback_map.contains(tile_id) && tile_id.zoom_level > 0)
        tile_id = tile_id.parent();
    if (!m_gpu_fallback_map.contains(tile_id))
        return { {}, 0 }; // may be empty during startup.
    return { tile_id, m_gpu_fallback_map.at(tile_id) };
}

void VectorLayer::setup_buffers(std::shared_ptr<ShaderProgram> shader, const std::vector<nucleus::tile::TileBounds>& draw_list, bool draw_fallback) const
{
    shader->set_uniform("acceleration_grid_sampler", 9);
    m_acceleration_grid_texture->bind(9);

    // upload all geometry buffers
    // binds the buffers to "...buffer_sampler_[0-max]"
    // NOTE: we moved the texture binding of the geometry_buffer_samples earlier since for some reason the last texture array never bound correctly if done
    // at the end of this function. We however leave the indices to the old ones to make the binding of the other textures easier
    constexpr unsigned triangle_vertex_buffer_start = 15;
    for (unsigned i = 0; i < constants::array_layer_tile_amount.size(); i++) {
        shader->set_uniform("geometry_buffer_sampler_" + std::to_string(i), triangle_vertex_buffer_start + i);
        m_geometry_buffer_texture[i]->bind(triangle_vertex_buffer_start + i);
    }

    nucleus::Raster<uint8_t> zoom_level_raster = { glm::uvec2 { 1024, 1 } };
    nucleus::Raster<glm::u16vec2> array_index_raster = { glm::uvec2 { 1024, 4 } };

    nucleus::Raster<uint8_t> zoom_level_raster_fallback = { glm::uvec2 { 1024, 1 } };
    nucleus::Raster<uint16_t> array_index_raster_fallback = { glm::uvec2 { 1024, 4 } };
    for (unsigned i = 0; i < std::min(unsigned(draw_list.size()), 1024u); ++i) {

        GpuMultiArrayHelper::MultiLayerInfo layer0;

        layer0 = m_gpu_multi_array_helper.layer(draw_list[i].id);

        if (layer0.id.zoom_level > 200)
            continue; // invalid zoom_level -> go to next
        zoom_level_raster.pixel({ i, 0 }) = layer0.id.zoom_level;
        array_index_raster.pixel({ i, 0 }) = { layer0.index1, layer0.index2 };

        // load the level offsets -> also make sure that all level offsets are loaded -> otherwise use parent level
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

        // fallback texture meta
        if (!draw_fallback) {
            // we only want to use the fallback texture when we render the sdf version (when drawing fallback we write to this texture and do not read from it)

            const auto fallback_layer0 = fallback_layer(layer0.id);
            if (fallback_layer0.id.zoom_level > 200)
                continue; // invalid zoom_level -> go to next

            zoom_level_raster_fallback.pixel({ i, 0 }) = fallback_layer0.id.zoom_level;
            array_index_raster_fallback.pixel({ i, 0 }) = fallback_layer0.index;

            const auto fallback_layer1 = fallback_layer(fallback_layer0.id.parent());
            if (fallback_layer1.id.zoom_level < 200) {
                array_index_raster_fallback.pixel({ i, 1 }) = fallback_layer1.index;
            } else {
                array_index_raster_fallback.pixel({ i, 1 }) = array_index_raster_fallback.pixel({ i, 0 });
            }
            const auto fallback_layer2 = fallback_layer(fallback_layer0.id.parent().parent());
            if (fallback_layer2.id.zoom_level < 200) {
                array_index_raster_fallback.pixel({ i, 2 }) = fallback_layer2.index;
            } else {
                array_index_raster_fallback.pixel({ i, 2 }) = array_index_raster_fallback.pixel({ i, 1 });
            }
            const auto fallback_layer3 = fallback_layer(fallback_layer0.id.parent().parent().parent());
            if (fallback_layer3.id.zoom_level < 200) {
                array_index_raster_fallback.pixel({ i, 3 }) = fallback_layer3.index;
            } else {
                array_index_raster_fallback.pixel({ i, 3 }) = array_index_raster_fallback.pixel({ i, 2 });
            }
        }
    }

    m_instanced_array_index->bind(10);
    shader->set_uniform("instanced_texture_array_index_sampler_vector", 10);
    m_instanced_array_index->upload(array_index_raster);

    m_instanced_zoom->bind(11);
    shader->set_uniform("instanced_texture_zoom_sampler_vector", 11);
    m_instanced_zoom->upload(zoom_level_raster);

    m_instanced_array_index_fallback->bind(12);
    shader->set_uniform("instanced_texture_array_index_sampler_vector_fallback", 12);
    m_instanced_array_index_fallback->upload(array_index_raster_fallback);

    m_instanced_zoom_fallback->bind(13);
    shader->set_uniform("instanced_texture_zoom_sampler_vector_fallback", 13);
    m_instanced_zoom_fallback->upload(zoom_level_raster_fallback);

    shader->set_uniform("styles_sampler", 14);
    m_styles_texture->bind(14);
}

void VectorLayer::draw(
    const TileGeometry& tile_geometry, const nucleus::camera::Definition& camera, const std::vector<nucleus::tile::TileBounds>& draw_list) const
{
    m_shader->bind();

    m_texture_layer->bind_buffer(m_shader, draw_list); // binds texture_sampler -> for ortho fotos

    m_shader->set_uniform("max_vector_geometry", m_max_vector_geometry);

    setup_buffers(m_shader, draw_list, false);

    constexpr uint8_t next_buffer_start = 15 + constants::array_layer_tile_amount.size();
    m_shader->set_uniform("fallback_texture_array_higher", next_buffer_start);
    m_fallback_texture_array_higher->bind(next_buffer_start);
    m_shader->set_uniform("fallback_texture_array_lower", next_buffer_start + 1);
    m_fallback_texture_array_lower->bind(next_buffer_start + 1);

    tile_geometry.draw(m_shader.get(), camera, draw_list);
}

void VectorLayer::update_gpu_tiles(const std::vector<nucleus::tile::Id>& deleted_tiles, const std::vector<nucleus::tile::GpuVectorLayerTile>& new_tiles)
{
    if (!QOpenGLContext::currentContext()) // can happen during shutdown.
        return;

    const auto current_time = nucleus::utils::time_since_epoch();

    for (const auto& quad : deleted_tiles) {
        for (const auto& id : quad.children()) {
            const auto exact_layer = m_gpu_multi_array_helper.exact_layer(0, id);
            if (exact_layer < m_vector_on_gpu.size()) { // make sure that we got a valid layer
                m_vector_on_gpu[exact_layer] = {};
                m_fallback_on_gpu[exact_layer] = { false, 0, 0, 0 };
            }

            m_gpu_multi_array_helper.remove_tile(id);
            m_gpu_fallback_map.erase(id);
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
        m_fallback_on_gpu[layer_index[0]] = { false, current_time, current_time - max_ms_for_ortho_wait, 0 };

        m_fallback_render_possible = true; // there was at least one new tile -> check if we need to render fallbacks
    }
}

void VectorLayer::update_fallback_textures(const std::vector<IdLayer>& tiles_to_render)
{
    // SETUP framebuffer
    if (m_fallback_framebuffer == unsigned(-1))
        return; // framebuffer not ready

    QOpenGLExtraFunctions* f = QOpenGLContext::currentContext()->extraFunctions();

    f->glDisable(GL_CULL_FACE);
    f->glDisable(GL_DEPTH_TEST);

    f->glBindFramebuffer(GL_FRAMEBUFFER, m_fallback_framebuffer);

    m_fallback_shader->bind();

    std::vector<nucleus::tile::TileBounds> bounds;

    std::transform(tiles_to_render.begin(), tiles_to_render.end(), std::back_inserter(bounds), [](const IdLayer& tile) { return tile.bounds; });

    m_texture_layer->bind_buffer(m_fallback_shader, bounds); // binds texture_sampler -> for ortho fotos

    m_fallback_shader->set_uniform("min_vector_geometry", m_max_vector_geometry);

    setup_buffers(m_fallback_shader, bounds, true);

    constexpr GLfloat clearAlbedoColor[] = { (242.0f * 0.7f) / 255.0f, (239.0f * 0.7f) / 255.0f, (233.0f * 0.7f) / 255.0f, 255.0f / 255.0f };
    const int mipmap_levels = 1 + static_cast<int>(std::floor(std::log2(m_fallback_resolution)));

    for (unsigned i = 0; i < tiles_to_render.size(); i++) {

        m_fallback_shader->set_uniform("lower_zoom", false);
        m_fallback_shader->set_uniform("instance_id", int(i));
        m_fallback_shader->set_uniform(
            "tile_id", glm::vec3(tiles_to_render[i].bounds.id.coords.x, tiles_to_render[i].bounds.id.coords.y, tiles_to_render[i].bounds.id.zoom_level));

        for (int level = 0; level < mipmap_levels; level++) {
            int mipmapped_resolution = m_fallback_resolution >> level;
            f->glViewport(0, 0, mipmapped_resolution, mipmapped_resolution);
            m_fallback_texture_array_higher->bind_layer_to_frame_buffer(0, tiles_to_render[i].layer, level);

            f->glClearBufferfv(GL_COLOR, 0, clearAlbedoColor);

            // it is possible to render 4 or 8 tiles at once (using different color attachments) -> but not sure if really performance relevant
            m_screen_quad_geometry.draw();
        }

        // lower -> one zoom level lower and with half the resolution of the higher
        m_fallback_shader->set_uniform("lower_zoom", true);
        for (int level = 1; level < mipmap_levels; level++) {
            int mipmapped_resolution = m_fallback_resolution >> level;
            f->glViewport(0, 0, mipmapped_resolution, mipmapped_resolution);
            m_fallback_texture_array_lower->bind_layer_to_frame_buffer(0, tiles_to_render[i].layer, level - 1);

            f->glClearBufferfv(GL_COLOR, 0, clearAlbedoColor);

            // it is possible to render 4 or 8 tiles at once (using different color attachments) -> but not sure if really performance relevant
            m_screen_quad_geometry.draw();
        }

        m_gpu_fallback_map[tiles_to_render[i].bounds.id] = tiles_to_render[i].layer;
    }

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

void VectorLayer::update_style(std::shared_ptr<const nucleus::Raster<glm::u32vec2>> styles)
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
    m_max_vector_geometry = int(new_max_vector_geometry); // uint uniforms do not work ...

    const auto current_time = nucleus::utils::time_since_epoch();

    for (unsigned i = 0; i < m_fallback_on_gpu.size(); i++) {
        m_fallback_on_gpu[i] = { false, current_time, current_time - max_ms_for_ortho_wait, 0 };
    }

    m_fallback_render_possible = true;
}

} // namespace gl_engine
