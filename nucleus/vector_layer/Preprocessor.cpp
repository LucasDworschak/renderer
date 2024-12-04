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
    const std::vector<std::vector<glm::vec2>> triangle_points = { { glm::vec2(10.5, 30.5), glm::vec2(30.5, 10.5), glm::vec2(50.5, 50.5) } };
    const std::vector<unsigned int> style_indices = { 1 };

    // // TODO somehow parse the data to lines and triangles

    auto processed_triangles = preprocess_triangles(triangle_points, style_indices);
    tile.data_triangle = std::make_shared<const std::vector<uint32_t>>(std::move(processed_triangles.data));

    condense_data(processed_triangles, processed_triangles, tile);

    return tile;
}

// polygon describe the outer edge of a closed shape
// -> neighbouring vertices form an edge
// last vertex connects to first vertex
VectorLayerCollection Preprocessor::preprocess_triangles(const std::vector<std::vector<glm::vec2>> polygons, const std::vector<unsigned int> style_indices)
{
    VectorLayerCollection triangle_collection;
    triangle_collection.cell_to_temp = std::vector<std::unordered_set<uint32_t>>(m_grid_size.x * m_grid_size.y, std::unordered_set<uint32_t>());

    float thickness = 5.0f;

    const auto grid_width = m_grid_size.x;
    const auto tile_bounds = m_tile_bounds;

    // TODO static assert to make sure that float are indeed the same size as uint32_t -> so as to not remove some important bits

    size_t data_offset = 0;

    // create the triangles from polygons
    for (size_t i = 0; i < polygons.size(); ++i) {

        // polygon to ordered triangles
        std::vector<glm::vec2> triangle_points = nucleus::utils::rasterizer::triangulize(polygons[i]);

        // add ordered triangles to collection
        for (size_t j = 0; j < triangle_points.size() / 3; ++j) {
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 0].x));
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 0].y));
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 1].x));
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 1].y));
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 2].x));
            triangle_collection.data.push_back(*reinterpret_cast<uint32_t*>(&triangle_points[j * 3 + 2].y));
            triangle_collection.data.push_back(style_indices[i]);
        }

        const auto cell_writer = [&triangle_collection, grid_width, tile_bounds, data_offset](glm::vec2 pos, int data_index) {
            if (tile_bounds.contains(pos)) {
                triangle_collection.cell_to_temp[int(pos.x) + grid_width * int(pos.y)].insert(data_index + data_offset);
            }
        };

        nucleus::utils::rasterizer::rasterize_triangle(cell_writer, triangle_points, thickness);

        // add to the data offset for the next polygon
        data_offset += triangle_points.size() / 3;
    }

    return triangle_collection;
}

VectorLayerCollection Preprocessor::preprocess_lines(const std::vector<std::vector<glm::vec2>> lines, const std::vector<unsigned int> style_indices)
{
    VectorLayerCollection line_collection;
    line_collection.cell_to_temp = std::vector<std::unordered_set<uint32_t>>(m_grid_size.x * m_grid_size.y, std::unordered_set<uint32_t>());

    float thickness = 5.0f;

    const auto grid_width = m_grid_size.x;
    const auto tile_bounds = m_tile_bounds;

    size_t data_offset = 0;

    // create the triangles from polygons
    for (size_t i = 0; i < lines.size(); ++i) {
        // create_lines(line_collection, lines[i], style_indices[i]);

        // TODO currently line points are duplicated, in the final shader -> ideally we only want to store one large point list for one line.
        // the problem here is how to also store meta data like style index efficiently (at start or end of line would be ideal but we would need another offset pointer to this location)
        // idea -> maybe store styleindex for lines and triangles separately (-> this would also probably be better for triangles)
        line_collection.data.reserve(lines[i].size() * 3);
        for (size_t j = 0; j < lines[i].size() - 1; j++) {
            glm::vec2 p0 = lines[i][j];
            glm::vec2 p1 = lines[i][j + 1];
            line_collection.data.push_back(*reinterpret_cast<uint32_t*>(&p0.x));
            line_collection.data.push_back(*reinterpret_cast<uint32_t*>(&p0.y));
            line_collection.data.push_back(*reinterpret_cast<uint32_t*>(&p1.x));
            line_collection.data.push_back(*reinterpret_cast<uint32_t*>(&p1.y));
            line_collection.data.push_back(style_indices[i]);
        }

        const auto cell_writer = [&line_collection, grid_width, tile_bounds, data_offset](glm::vec2 pos, int data_index) {
            if (tile_bounds.contains(pos)) {
                line_collection.cell_to_temp[int(pos.x) + grid_width * int(pos.y)].insert(data_index + data_offset);
            }
        };

        nucleus::utils::rasterizer::rasterize_line(cell_writer, lines[i], thickness);

        // add to the data offset for the next polyline
        // TODO test if this is the correct offset (we might be off by one here but it should suffice for now)
        data_offset += lines[i].size();
    }

    return line_collection;
}

/*
 * Function condenses data and fills the GpuVectorLayerTile.
 * condensing:
 *      go over every cell and gather distinct indices
 *      example input(triangle indice for 6 cells): [1], [1], [1,2], [2,3], [2,3], [2,3]
 *      expected output: [[1], [1,2], [2,3]]
 *      Note: it is theoretically possible to further condense the above example to [[1,2,3]] where we can both point to 1, 12 and 23
 *      nevertheless for more complex entries this might be overkill and take more time to compute than it is worth -> only necessary if we need more buffer space
 * simultaneously we also generate the final grid for the tile that stores the offset and size for lookups into the index_bridge
 */
void Preprocessor::condense_data(VectorLayerCollection& triangle_collection, VectorLayerCollection&, GpuVectorLayerTile& tile)
{

    std::unordered_map<std::unordered_set<uint32_t>, uint32_t, Hasher> unique_entries;
    std::vector<uint32_t> grid;
    std::vector<uint32_t> index_bridge;

    uint32_t start_offset = 0;

    for (size_t i = 0; i < triangle_collection.cell_to_temp.size(); ++i) {
        // go through every cell
        if (triangle_collection.cell_to_temp[i].size() == 0) {
            grid.push_back(0); // no data -> only add an emtpy cell
        } else {

            if (unique_entries.contains(triangle_collection.cell_to_temp[i])) {
                // we already have such an entry -> get the offset_size from there
                grid.push_back(unique_entries[triangle_collection.cell_to_temp[i]]);
            } else {
                // we have to add a new element
                const auto offset_size = nucleus::utils::bit_coding::u24_u8_to_u32(start_offset, uint8_t(triangle_collection.cell_to_temp[i].size()));
                unique_entries[triangle_collection.cell_to_temp[i]] = offset_size;

                grid.push_back(offset_size);

                // add the indices to the bridge
                index_bridge.insert(index_bridge.end(), triangle_collection.cell_to_temp[i].begin(), triangle_collection.cell_to_temp[i].end());

                start_offset += triangle_collection.cell_to_temp[i].size();
            }
        }
    }

    tile.grid_triangle = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(m_grid_size.x, std::move(grid)));
    tile.grid_to_data = std::make_shared<const std::vector<uint32_t>>(std::move(index_bridge));

    // TODO do the same for lines
}

} // namespace nucleus::vector_layer
