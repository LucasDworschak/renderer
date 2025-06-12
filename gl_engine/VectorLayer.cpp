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
#include "TileGeometry.h"
#include <QOpenGLExtraFunctions>

#include "nucleus/vector_layer/constants.h"

namespace gl_engine {

using namespace nucleus::vector_layer;

VectorLayer::VectorLayer(QObject* parent)
    : QObject { parent }
    , m_gpu_multi_array_helper()
    , m_initialized(false)

{
}

void gl_engine::VectorLayer::init(ShaderRegistry* shader_registry)
{
    std::vector<QString> defines;
    defines.push_back(QString("#define style_bits %1").arg(constants::style_bits));
    defines.push_back(QString("#define style_precision %1").arg(constants::style_precision));
    defines.push_back(QString("#define max_zoom %1").arg(constants::style_zoom_range.y));
    defines.push_back(QString("#define zoom_blend_steps %1").arg(constants::style_zoom_blend_steps));
    defines.push_back(QString("#define tile_extent %1").arg(constants::tile_extent));
    defines.push_back(QString("#define tile_scale %1").arg(constants::tile_scale));
    defines.push_back(QString("#define grid_size vec2(%1,%1)").arg(constants::grid_size));

    defines.push_back(QString("#define all_bits %1").arg(constants::all_bits));
    defines.push_back(QString("#define cell_bits %1").arg(constants::cell_bits));
    defines.push_back(QString("#define coordinate_bits %1").arg(constants::coordinate_bits));
    defines.push_back(QString("#define aa_border %1").arg(constants::aa_border));

    defines.push_back(QString("#define sampler_offset %1").arg(constants::array_helper_all_bits - constants::array_helper_buffer_info_bits));

    m_shader = std::make_shared<ShaderProgram>("tile.vert", "vector_layer.frag", ShaderCodeSource::FILE, defines);
    shader_registry->add_shader(m_shader);

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

    m_initialized = true;
    if (m_styles) {
        m_styles_texture->upload(*m_styles);
    }
}

unsigned VectorLayer::tile_count() const { return m_gpu_multi_array_helper.n_occupied(); }

void VectorLayer::draw(
    const TileGeometry& tile_geometry, const nucleus::camera::Definition& camera, const std::vector<nucleus::tile::TileBounds>& draw_list) const
{
    m_shader->bind();
    m_shader->set_uniform("acceleration_grid_sampler", 2);
    m_acceleration_grid_texture->bind(2);

    nucleus::Raster<uint8_t> zoom_level_raster = { glm::uvec2 { 1024, 1 } };
    nucleus::Raster<glm::u16vec2> array_index_raster = { glm::uvec2 { 1024, 1 } };
    for (unsigned i = 0; i < std::min(unsigned(draw_list.size()), 1024u); ++i) {
        const auto layer = m_gpu_multi_array_helper.layer(draw_list[i].id);
        zoom_level_raster.pixel({ i, 0 }) = layer.id.zoom_level;
        array_index_raster.pixel({ i, 0 }) = { layer.index1, layer.index2 };
    }

    m_instanced_array_index->bind(5);
    m_shader->set_uniform("instanced_texture_array_index_sampler", 5);
    m_instanced_array_index->upload(array_index_raster);

    m_instanced_zoom->bind(6);
    m_shader->set_uniform("instanced_texture_zoom_sampler", 6);
    m_instanced_zoom->upload(zoom_level_raster);

    m_shader->set_uniform("styles_sampler", 7);
    m_styles_texture->bind(7);

    // upload all geometry buffers
    // binds the buffers to "...buffer_sampler_[0-max]"
    constexpr uint8_t triangle_vertex_buffer_start = 8;
    for (uint8_t i = 0; i < constants::array_layer_tile_amount.size(); i++) {
        m_shader->set_uniform("geometry_buffer_sampler_" + std::to_string(i), triangle_vertex_buffer_start + i);
        m_geometry_buffer_texture[i]->bind(triangle_vertex_buffer_start + i);
    }

    tile_geometry.draw(m_shader.get(), camera, draw_list);
}

void VectorLayer::update_gpu_tiles(const std::vector<nucleus::tile::Id>& deleted_tiles, const std::vector<nucleus::tile::GpuVectorLayerTile>& new_tiles)
{
    if (!QOpenGLContext::currentContext()) // can happen during shutdown.
        return;

    for (const auto& quad : deleted_tiles) {
        for (const auto& id : quad.children()) {
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
    }
}

void VectorLayer::set_tile_limit(unsigned int new_limit)
{
    assert(new_limit <= 2048); // array textures with size > 2048 are not supported on all devices
    assert(!m_acceleration_grid_texture);
    assert(m_geometry_buffer_texture.size() == 0);

    m_gpu_multi_array_helper.set_tile_limit(new_limit);
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

} // namespace gl_engine
