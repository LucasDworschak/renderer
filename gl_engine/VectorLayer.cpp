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

VectorLayer::VectorLayer(QObject* parent)
    : QObject { parent }
{
}

void gl_engine::VectorLayer::init(ShaderRegistry* shader_registry)
{
    m_shader = std::make_shared<ShaderProgram>("tile.vert", "vector_layer.frag");
    shader_registry->add_shader(m_shader);

    m_grid_texture = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::R32UI);
    m_grid_texture->setParams(gl_engine::Texture::Filter::Nearest, gl_engine::Texture::Filter::Nearest);
    // // TODO: might become larger than GL_MAX_ARRAY_TEXTURE_LAYERS
    m_grid_texture->allocate_array(nucleus::vector_layer::constants::grid_size, nucleus::vector_layer::constants::grid_size, unsigned(m_gpu_array_helper.size()));

    m_tile_id_texture = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::RG32UI);
    m_tile_id_texture->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_array_index_texture = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::R16UI);
    m_array_index_texture->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_triangle_index_texture = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::R32UI);
    m_triangle_index_texture->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);
    m_triangle_index_texture->allocate_array(nucleus::vector_layer::constants::data_size, nucleus::vector_layer::constants::data_size, unsigned(m_gpu_array_helper.size()));

    m_triangle_data_texture = std::make_unique<Texture>(Texture::Target::_2dArray, Texture::Format::R32UI);
    m_triangle_data_texture->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);
    m_triangle_data_texture->allocate_array(nucleus::vector_layer::constants::data_size, nucleus::vector_layer::constants::data_size, unsigned(m_gpu_array_helper.size()));

    QOpenGLFunctions* f = QOpenGLContext::currentContext()->functions();
    GLint max_layers;
    GLint max_size;
    f->glGetIntegerv(GL_MAX_ARRAY_TEXTURE_LAYERS, &max_layers);
    f->glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max_size);
    std::cout << "max layers: " << max_layers << std::endl;
    std::cout << "max size: " << max_size << std::endl;

    update_gpu_id_map();
}

void VectorLayer::draw(
    const TileGeometry& tile_geometry, const nucleus::camera::Definition& camera, const nucleus::tile::DrawListGenerator::TileSet& draw_tiles, bool sort_tiles, glm::dvec3 sort_position) const
{
    m_shader->bind();
    m_shader->set_uniform("grid_sampler", 2);
    m_grid_texture->bind(2);

    m_shader->set_uniform("grid_index_sampler", 5);
    m_array_index_texture->bind(5);
    m_shader->set_uniform("vector_map_tile_id_sampler", 6);
    m_tile_id_texture->bind(6);

    m_shader->set_uniform("triangle_index_sampler", 7);
    m_triangle_index_texture->bind(7);

    m_shader->set_uniform("triangle_data_sampler", 8);
    m_triangle_data_texture->bind(8);

    tile_geometry.draw(m_shader.get(), camera, draw_tiles, sort_tiles, sort_position);
}

void VectorLayer::update_gpu_quads(const std::vector<nucleus::tile::GpuVectorLayerQuad>& new_quads, const std::vector<nucleus::tile::Id>& deleted_quads)
{
    if (!QOpenGLContext::currentContext()) // can happen during shutdown.
        return;

    for (const auto& quad : deleted_quads) {
        for (const auto& id : quad.children()) {

            if (!m_gpu_array_helper.contains_tile(id))
                continue; // TODO maybe better to move the contains inside m_gpu_array_helper.remove_tile (but currently do not want to change the assert)

            m_gpu_array_helper.remove_tile(id);
        }
    }

    for (const auto& quad : new_quads) {
        for (const auto& tile : quad.tiles) {
            // test for validity
            assert(tile.id.zoom_level < 100);
            if (!tile.grid_triangle)
                continue; // nothing here

            // assert(tile.grid_triangle);
            assert(tile.grid_to_data.size() > 0);
            assert(tile.data_triangle.size() > 0);

            const auto max_layer_amount = std::max(tile.grid_to_data.size(), tile.data_triangle.size());

            // find empty spot and upload texture
            for (size_t i = 0; i < max_layer_amount; i++) {
                const auto layer_index = m_gpu_array_helper.add_tile(tile.id, i);
                if (i == 0) // only one grid texture needed
                    m_grid_texture->upload(*tile.grid_triangle, layer_index);

                // upload index and data if there is still data to upload
                if (i < tile.grid_to_data.size())
                    m_triangle_index_texture->upload(*tile.grid_to_data[i], layer_index);
                if (i < tile.data_triangle.size())
                    m_triangle_data_texture->upload(*tile.data_triangle[i], layer_index);
            }

            // int count = 0;
            // for (auto data : tile.data_triangle->buffer()) {
            //     if (count++ > 15)
            //         break;
            //     std::cout << *reinterpret_cast<float*>(&data) << "\t";
            // }
            // std::cout << std::endl;

            // // int count = 0;
            // for (auto data : tile.grid_to_data->buffer()) {
            //     if (data == 0 || data == -1u)
            //         // count++;
            //         // break;
            //         std::cout << data << "\t";
            // }
            // std::cout << std::endl << std::endl;
        }
    }

    update_gpu_id_map();
}

void VectorLayer::set_quad_limit(unsigned int new_limit) { m_gpu_array_helper.set_quad_limit(new_limit); }

void VectorLayer::update_gpu_id_map()
{
    auto [packed_ids, layers] = m_gpu_array_helper.generate_dictionary();
    m_array_index_texture->upload(layers);
    m_tile_id_texture->upload(packed_ids);
}

} // namespace gl_engine
