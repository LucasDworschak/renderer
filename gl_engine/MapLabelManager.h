/*****************************************************************************
 * Alpine Terrain Renderer
 * Copyright (C) 2024 Lucas Dworschak
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

#include <QOpenGLBuffer>
#include <QOpenGLTexture>
#include <QOpenGLVertexArrayObject>
#include <unordered_map>

#include "Framebuffer.h"
#include "Texture.h"
#include "nucleus/camera/Definition.h"
#include "nucleus/map_label/FilterDefinitions.h"
#include "nucleus/map_label/LabelFactory.h"

#include "nucleus/tile_scheduler/DrawListGenerator.h"

using namespace nucleus::vector_tile;

namespace gl_engine {
class ShaderProgram;
class ShaderRegistry;

struct GPUVectorTile {
    tile::Id id;
    std::unique_ptr<QOpenGLBuffer> vertex_buffer;
    std::unique_ptr<QOpenGLVertexArrayObject> vao;
    size_t instance_count; // how many characters (+1 for icon)
    glm::dvec3 reference_point = {};
};

class MapLabelManager : public QObject {
    Q_OBJECT

public:
    using TileSet = nucleus::tile_scheduler::DrawListGenerator::TileSet;
    explicit MapLabelManager(const nucleus::tile_scheduler::utils::AabbDecoratorPtr& aabb_decorator, QObject* parent = nullptr);

    void init(ShaderRegistry* shader_registry);
    void draw(Framebuffer* gbuffer, const nucleus::camera::Definition& camera, const TileSet& draw_tiles) const;
    void draw_picker(Framebuffer* gbuffer, const nucleus::camera::Definition& camera, const TileSet& draw_tiles) const;
    TileSet generate_draw_list(const nucleus::camera::Definition& camera) const;

    void update_labels(const std::vector<nucleus::vector_tile::PoiTile>& updated_tiles, const std::vector<tile::Id>& removed_tiles);

private:
    void upload_to_gpu(const tile::Id& id, const PointOfInterestCollection& features);
    void remove_tile(const tile::Id& tile_id);

    std::shared_ptr<ShaderProgram> m_label_shader;
    std::shared_ptr<ShaderProgram> m_picker_shader;

    std::unique_ptr<Texture> m_font_texture;
    std::unique_ptr<Texture> m_icon_texture;

    std::unique_ptr<QOpenGLBuffer> m_index_buffer;
    size_t m_indices_count; // how many vertices per character (most likely 6 since quads)

    nucleus::maplabel::LabelFactory m_mapLabelFactory;

    nucleus::tile_scheduler::DrawListGenerator m_draw_list_generator;
    std::unordered_map<tile::Id, std::shared_ptr<GPUVectorTile>, tile::Id::Hasher> m_gpu_tiles;
};
} // namespace gl_engine
