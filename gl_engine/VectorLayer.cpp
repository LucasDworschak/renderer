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
#include "nucleus/utils/bit_coding.h"
#include <QOpenGLExtraFunctions>

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
    m_grid_texture->allocate_array(GRID_RESOLUTION, GRID_RESOLUTION, unsigned(m_gpu_array_helper.size()));

    m_tile_id_texture = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::RG32UI);
    m_tile_id_texture->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_meta_texture = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::RG32UI);
    m_meta_texture->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_triangle_index_texture = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::R32UI);
    m_triangle_index_texture->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    m_triangle_data_texture = std::make_unique<Texture>(Texture::Target::_2d, Texture::Format::R32UI);
    m_triangle_data_texture->setParams(Texture::Filter::Nearest, Texture::Filter::Nearest);

    // initial texture upload with default values (but set width and format)
    m_tile_id_texture->upload(nucleus::Raster<glm::u32vec2>({ 256, 256 }, glm::u32vec2(-1, -1)));
    m_meta_texture->upload(nucleus::Raster<glm::u32vec2>({ 256, 256 }, glm::u32vec2(-1, -1)));
    m_triangle_index_texture->upload(nucleus::Raster<uint32_t>({ TEXTURE_RESOLUTION, TEXTURE_RESOLUTION }, -1u));
    m_triangle_data_texture->upload(nucleus::Raster<uint32_t>({ TEXTURE_RESOLUTION, TEXTURE_RESOLUTION }, -1u));
}

void VectorLayer::draw(
    const TileGeometry& tile_geometry, const nucleus::camera::Definition& camera, const nucleus::tile::DrawListGenerator::TileSet& draw_tiles, bool sort_tiles, glm::dvec3 sort_position) const
{
    m_shader->bind();
    m_shader->set_uniform("grid_sampler", 2);
    m_grid_texture->bind(2);

    m_shader->set_uniform("vector_map_meta_sampler", 5);
    m_meta_texture->bind(5);
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
                continue; // TODO maybe better to move the contains inside remove_tile (but currently do not want to change the assert)

            m_gpu_array_helper.remove_tile(id);

            assert(m_id_to_data_bridge.contains(id));
            m_id_to_data_bridge.erase(id);
            assert(m_id_to_triangle_data.contains(id));
            m_id_to_triangle_data.erase(id);
        }
    }

    for (const auto& quad : new_quads) {
        for (const auto& tile : quad.tiles) {
            // test for validity
            assert(tile.id.zoom_level < 100);
            if (!tile.grid_triangle)
                continue; // nothing here

            // assert(tile.grid_triangle);
            assert(tile.grid_to_data);
            assert(tile.data_triangle);

            assert(!m_id_to_data_bridge.contains(tile.id));
            m_id_to_data_bridge.emplace(tile.id, tile.grid_to_data);

            assert(!m_id_to_triangle_data.contains(tile.id));
            m_id_to_triangle_data.emplace(tile.id, tile.data_triangle);

            // find empty spot and upload texture
            const auto layer_index = m_gpu_array_helper.add_tile(tile.id);
            m_grid_texture->upload(*tile.grid_triangle, layer_index);
        }
    }

    update_gpu_data();
}

void VectorLayer::set_quad_limit(unsigned int new_limit) { m_gpu_array_helper.set_quad_limit(new_limit); }

void VectorLayer::update_gpu_data()
{
    // instead of straight up uploading the m_array_index_texture (like TextureLayer),
    // we first want to add the bridge and triangle offsets to the array_index

    // - use m_gpu_array_helper to get layer id and upload the grid as is.
    // - save the data offsets per tile for the data_triangle
    // - combine and reupload data_triangles
    // - combine and reupload grid_to_data
    // - upload m_meta_texture and m_triangle_index_texture and m_triangle_data_texture

    auto [packed_ids, layers] = m_gpu_array_helper.generate_dictionary();
    m_tile_id_texture->reupload(packed_ids);

    // instead of uploading layers directly, we combine the data with global offsets and save it into a meta texture
    // layer -> 16bits
    // triangle_offset -> 24bits
    // data_index_offset -> 24bits
    // = 2*32 bits
    nucleus::Raster<glm::u32vec2> meta({ 256, 256 }, glm::u32vec2(-1, -1)); // RG32UI

    // further more for the index indirection map and the triangle data, we combine the data of all tiles into the following vectors
    // and upload them to their respective textures
    std::vector<uint32_t> triangle_combined_indices;
    triangle_combined_indices.reserve(TEXTURE_RESOLUTION * TEXTURE_RESOLUTION);
    std::vector<uint32_t> triangle_combined_data;
    triangle_combined_data.reserve(TEXTURE_RESOLUTION * TEXTURE_RESOLUTION);

    uint32_t data_index_offset = 0;
    uint32_t triangle_offset = 0;

    const auto hash_to_pixel_255 = [](uint16_t hash) { return glm::uvec2(hash & 255, hash >> 8); };
    for (const auto& [id, data_bridge] : m_id_to_data_bridge) {
        auto hash = nucleus::srs::hash_uint16(id);
        while (meta.pixel(hash_to_pixel_255(hash)) != glm::u32vec2(-1, -1))
            hash++;

        meta.pixel(hash_to_pixel_255(hash)) = nucleus::utils::bit_coding::u16_u24_u24_to_u32_2(layers.pixel(hash_to_pixel_255(hash)), data_index_offset, triangle_offset);

        // add index and triangle data to combined vectors
        triangle_combined_indices.insert(triangle_combined_indices.end(), data_bridge->begin(), data_bridge->end());

        triangle_combined_data.insert(triangle_combined_data.end(), m_id_to_triangle_data.at(id)->begin(), m_id_to_triangle_data.at(id)->end());

        data_index_offset += data_bridge->size();
        triangle_offset += m_id_to_triangle_data.at(id)->size();
    }

    m_meta_texture->reupload(meta);

    // make sure that we have exactly TEXTURE_RESOLUTIONxTEXTURE_RESOLUTION elements
    assert(triangle_combined_indices.size() <= TEXTURE_RESOLUTION * TEXTURE_RESOLUTION);
    assert(triangle_combined_data.size() <= TEXTURE_RESOLUTION * TEXTURE_RESOLUTION);
    triangle_combined_indices.resize(TEXTURE_RESOLUTION * TEXTURE_RESOLUTION, -1u);
    triangle_combined_data.resize(TEXTURE_RESOLUTION * TEXTURE_RESOLUTION, -1u);

    // pack the data into square raster and upload them
    nucleus::Raster<uint32_t> triangle_index(TEXTURE_RESOLUTION, std::move(triangle_combined_indices)); // R32UI
    nucleus::Raster<uint32_t> triangle_data(TEXTURE_RESOLUTION, std::move(triangle_combined_data)); // R32UI

    m_triangle_index_texture->reupload(triangle_index);
    m_triangle_data_texture->reupload(triangle_data);
}

} // namespace gl_engine
