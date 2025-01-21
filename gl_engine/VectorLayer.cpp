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
{
}

void gl_engine::VectorLayer::init(ShaderRegistry* shader_registry)
{
    m_shader = std::make_shared<ShaderProgram>("tile.vert", "vector_layer.frag");
    shader_registry->add_shader(m_shader);

    m_triangle_acceleration_grid_texture = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::R32UI);
    m_triangle_acceleration_grid_texture->setParams(gl_engine::Texture::Filter::Nearest, gl_engine::Texture::Filter::Nearest);
    // // TODO: might become larger than GL_MAX_ARRAY_TEXTURE_LAYERS
    m_triangle_acceleration_grid_texture->allocate_array(constants::grid_size, constants::grid_size, unsigned(m_gpu_multi_array_helper.layer_amount(0)));

    m_tile_id_texture = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::RG32UI);
    m_tile_id_texture->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_array_index_texture = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::RG16UI);
    m_array_index_texture->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_triangle_index_buffer_texture.resize(constants::array_layer_quad_amount.size());
    m_triangle_vertex_buffer_texture.resize(constants::array_layer_quad_amount.size());
    for (uint8_t i = 0; i < constants::array_layer_quad_amount.size(); i++) {
        const auto layer_amount = m_gpu_multi_array_helper.layer_amount(i);

        m_triangle_index_buffer_texture[i] = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::R32UI);
        m_triangle_index_buffer_texture[i]->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);
        m_triangle_index_buffer_texture[i]->allocate_array(constants::data_size[i], constants::data_size[i], layer_amount);

        m_triangle_vertex_buffer_texture[i] = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::R32UI);
        m_triangle_vertex_buffer_texture[i]->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);
        m_triangle_vertex_buffer_texture[i]->allocate_array(constants::data_size[i], constants::data_size[i], layer_amount);
    }

    update_gpu_id_map();
}

void VectorLayer::draw(
    const TileGeometry& tile_geometry, const nucleus::camera::Definition& camera, const nucleus::tile::DrawListGenerator::TileSet& draw_tiles, bool sort_tiles, glm::dvec3 sort_position) const
{
    m_shader->bind();
    m_shader->set_uniform("triangle_acceleration_grid_sampler", 2);
    m_triangle_acceleration_grid_texture->bind(2);

    m_shader->set_uniform("array_index_sampler", 5);
    m_array_index_texture->bind(5);
    m_shader->set_uniform("vector_map_tile_id_sampler", 6);
    m_tile_id_texture->bind(6);

    // upload all index and vertex buffer
    // binds the buffers to "...buffer_sampler_[0-max]"
    constexpr uint8_t triangle_index_buffer_start = 7;
    constexpr uint8_t triangle_vertex_buffer_start = triangle_index_buffer_start + constants::array_layer_quad_amount.size();
    for (uint8_t i = 0; i < constants::array_layer_quad_amount.size(); i++) {
        m_shader->set_uniform("triangle_index_buffer_sampler_" + std::to_string(i), triangle_index_buffer_start + i);
        m_triangle_index_buffer_texture[i]->bind(triangle_index_buffer_start + i);

        m_shader->set_uniform("triangle_vertex_buffer_sampler_" + std::to_string(i), triangle_vertex_buffer_start + i);
        m_triangle_vertex_buffer_texture[i]->bind(triangle_vertex_buffer_start + i);
    }

    tile_geometry.draw(m_shader.get(), camera, draw_tiles, sort_tiles, sort_position);
}

void VectorLayer::update_gpu_quads(const std::vector<nucleus::tile::GpuVectorLayerQuad>& new_quads, const std::vector<nucleus::tile::Id>& deleted_quads)
{
    if (!QOpenGLContext::currentContext()) // can happen during shutdown.
        return;

    for (const auto& quad : deleted_quads) {
        for (const auto& id : quad.children()) {
            m_gpu_multi_array_helper.remove_tile(id);
        }
    }

    for (const auto& quad : new_quads) {
        for (const auto& tile : quad.tiles) {
            // test for validity
            assert(tile.id.zoom_level < 100);
            if (!tile.triangle_acceleration_grid)
                continue; // nothing here

            // assert(tile.triangle_acceleration_grid);
            assert(tile.triangle_index_buffer);
            assert(tile.triangle_vertex_buffer);

            // find empty spot and upload texture
            const auto layer_index = m_gpu_multi_array_helper.add_tile(tile.id, tile.triangle_buffer_info);

            m_triangle_acceleration_grid_texture->upload(*tile.triangle_acceleration_grid, layer_index[0]);

            m_triangle_index_buffer_texture[tile.triangle_buffer_info]->upload(*tile.triangle_index_buffer, layer_index[1]);
            m_triangle_vertex_buffer_texture[tile.triangle_buffer_info]->upload(*tile.triangle_vertex_buffer, layer_index[1]);
        }
    }

    update_gpu_id_map();
}

void VectorLayer::set_quad_limit(unsigned int new_limit) { m_gpu_multi_array_helper.set_quad_limit(new_limit); }
void VectorLayer::update_gpu_id_map()
{
    auto [packed_ids, layers] = m_gpu_multi_array_helper.generate_dictionary();
    m_array_index_texture->upload(layers);
    m_tile_id_texture->upload(packed_ids);
}

} // namespace gl_engine
