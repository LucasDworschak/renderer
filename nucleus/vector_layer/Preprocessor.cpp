/*****************************************************************************
 * AlpineMaps.org
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

#include "Preprocessor.h"

#include <QImage>
#include <QString>

#include <CDT.h>

#include "nucleus/Raster.h"
#include <nucleus/utils/bit_coding.h>

#include "nucleus/utils/rasterizer.h"

namespace nucleus::vector_layer {

Preprocessor::Preprocessor(nucleus::tile::Id id)
    : m_x_values_per_y_step(m_grid_size.y, 0)
{
    m_tile_up_direction = (id.scheme == tile::Scheme::Tms) ? 1 : -1;
    m_tile_bounds = { glm::dvec2(0, 0), glm::dvec2 { m_grid_size.x, m_grid_size.y } }; // TODO calculate with srs.tile_bounds method
}

GpuVectorLayerTile Preprocessor::preprocess(const tile::Data data)
{
    GpuVectorLayerTile tile;
    tile.id = data.id;

    // DEBUG polygons
    const std::vector<std::vector<glm::vec2>> triangle_points = { { glm::vec2(10, 30), glm::vec2(30, 10), glm::vec2(50, 50) } };
    const std::vector<unsigned int> style_indices = { 1 };

    // // TODO somehow parse the data to lines and triangles

    // tile.data_line = std::make_shared<VectorLayerLineCollection>(preprocess_lines({}, {}));
    // std::make_shared<VectorLayerTriangleCollection>(preprocess_triangles(triangle_points, style_indices));
    auto processed_triangles = preprocess_triangles(triangle_points, style_indices);

    { // triangles to GpuVectorLayerTile
        std::vector<uint32_t> triangles;
        triangles.reserve(processed_triangles.triangles.size() * 7);

        for (size_t i = 0; i < processed_triangles.triangles.size(); ++i) {
            triangles.insert(triangles.end(), processed_triangles.triangles[i].packed, processed_triangles.triangles[i].packed + 7);
        }

        tile.data_triangle = std::make_shared<const std::vector<uint32_t>>(std::move(triangles));
    }

    { // grid to GpuVectorLayerTile
        std::vector<uint32_t> grid;
        std::vector<uint32_t> grid_to_data;

        uint32_t start_offset = 0;

        for (size_t i = 0; i < processed_triangles.cell_to_data.size(); ++i) {

            // grid is a singular uint32_t value that encodes the start index of the triangle list and the amount of triangles
            if (processed_triangles.cell_to_data[i].size() == 0) {
                grid.push_back(0);
            } else {
                grid.push_back(nucleus::utils::bit_coding::u24_u8_to_u32(start_offset, uint8_t(processed_triangles.cell_to_data[i].size())));
                start_offset += processed_triangles.cell_to_data[i].size();
                grid_to_data.insert(grid_to_data.end(), processed_triangles.cell_to_data[i].begin(), processed_triangles.cell_to_data[i].end());
            }
        }

        tile.grid_triangle = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(m_grid_size.x, std::move(grid)));
        tile.grid_to_data = std::make_shared<const std::vector<uint32_t>>(std::move(grid_to_data));
    }

    return tile;
}

// polygon describe the outer edge of a closed shape
// -> neighbouring vertices form an edge
// last vertex connects to first vertex
VectorLayerTriangleCollection Preprocessor::preprocess_triangles(const std::vector<std::vector<glm::vec2>> polygons, const std::vector<unsigned int> style_indices)
{
    VectorLayerTriangleCollection triangle_collection;
    triangle_collection.cell_to_data = std::vector<std::unordered_set<uint32_t>>(m_grid_size.x * m_grid_size.y, std::unordered_set<uint32_t>());

    float thickness = 5.0f;

    const auto grid_width = m_grid_size.x;
    const auto tile_bounds = m_tile_bounds;

    const auto cell_writer = [&triangle_collection, grid_width, tile_bounds](glm::vec2 pos, int data_index) {
        if (tile_bounds.contains(pos)) {
            triangle_collection.cell_to_data[int(pos.x) + grid_width * int(pos.y)].insert(data_index);
        }
    };

    // create the triangles from polygons
    for (size_t i = 0; i < polygons.size(); ++i) {
        // polygon to ordered triangles
        std::vector<glm::vec2> triangles = nucleus::utils::rasterizer::triangulize(polygons[i]);

        // add ordered triangles to collection
        for (size_t j = 0; j < triangles.size() / 3; ++j) {
            triangle_collection.triangles.push_back({ triangles[j * 3 + 0], triangles[j * 3 + 1], triangles[j * 3 + 2], style_indices[i] });
        }

        nucleus::utils::rasterizer::rasterize_triangle(cell_writer, triangles, thickness);
    }

    return triangle_collection;
}

VectorLayerLineCollection Preprocessor::preprocess_lines(const std::vector<std::vector<glm::vec2>> lines, const std::vector<unsigned int> style_indices)
{
    VectorLayerLineCollection line_collection;
    line_collection.cell_to_data = std::vector<std::unordered_set<uint32_t>>(m_grid_size.x * m_grid_size.y, std::unordered_set<uint32_t>());

    float thickness = 5.0f;

    const auto grid_width = m_grid_size.x;
    const auto tile_bounds = m_tile_bounds;

    const auto cell_writer = [&line_collection, grid_width, tile_bounds](glm::vec2 pos, int data_index) {
        if (tile_bounds.contains(pos)) {
            line_collection.cell_to_data[int(pos.x) + grid_width * int(pos.y)].insert(data_index);
        }
    };

    // create the triangles from polygons
    for (size_t i = 0; i < lines.size(); ++i) {
        create_lines(line_collection, lines[i], style_indices[i]);

        nucleus::utils::rasterizer::rasterize_line(cell_writer, lines[i], thickness);
    }

    return line_collection;
}

void Preprocessor::create_lines(VectorLayerLineCollection& line_collection, const std::vector<glm::vec2> line_points, unsigned int style_index)
{
    // TODO currently line points are duplicated, in the final shader -> ideally we only want to store one large point list for one line.
    line_collection.lines.reserve(line_points.size());
    for (size_t i = 0; i < line_points.size() - 1; i++) {
        line_collection.lines.push_back({ line_points[i], line_points[i + 1], style_index });
    }
}

} // namespace nucleus::vector_layer
