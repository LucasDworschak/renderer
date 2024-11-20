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

    float dist_to_grid_center = m_tile_bounds.size().x / m_grid_size.x / sqrt(2);
    float thickness = 5.0f + dist_to_grid_center;

    // create the triangles from polygons
    for (size_t i = 0; i < polygons.size(); ++i) {
        create_triangles(triangle_collection, polygons[i], style_indices[i]);
    }

    // fill the grid from the created triangles
    for (size_t i = 0; i < triangle_collection.triangles.size(); ++i) {
        dda_triangle(triangle_collection, i, thickness);
    }

    return triangle_collection;
}

VectorLayerLineCollection Preprocessor::preprocess_lines(const std::vector<std::vector<glm::vec2>> lines, const std::vector<unsigned int> style_indices)
{
    VectorLayerLineCollection line_collection;
    line_collection.cell_to_data = std::vector<std::unordered_set<uint32_t>>(m_grid_size.x * m_grid_size.y, std::unordered_set<uint32_t>());

    float dist_to_grid_center = m_tile_bounds.size().x / m_grid_size.x / sqrt(2);
    float thickness = 5.0f + dist_to_grid_center;

    // create the triangles from polygons
    for (size_t i = 0; i < lines.size(); ++i) {
        create_lines(line_collection, lines[i], style_indices[i]);
    }

    // fill the grid from the created triangles
    for (size_t i = 0; i < line_collection.lines.size(); ++i) {
        auto vertex_start = line_collection.lines[i].data.line_start_vertex;
        auto vertex_end = line_collection.lines[i].data.line_end_vertex;
        glm::vec2 edge = vertex_end - vertex_start;
        glm::vec2 normal = glm::normalize(glm::vec2(-edge.y, edge.x));

        dda_line(line_collection.cell_to_data, vertex_start, edge, normal * thickness, line_collection.lines[i].data.style_index, 0, false);
        // dda_line(vertex_start, edge, normal * 1.0f, 100, 0, false); // DEBUG draw only line without thickness

        add_end_cap(line_collection.cell_to_data, vertex_start, line_collection.lines[i].data.style_index, thickness / 2.0);
        add_end_cap(line_collection.cell_to_data, vertex_end, line_collection.lines[i].data.style_index, thickness / 2.0);
    }

    return line_collection;
}

void Preprocessor::create_triangles(VectorLayerTriangleCollection& triangle_collection, const std::vector<glm::vec2> polygon_points, unsigned int style_index)
{
    std::vector<glm::ivec2> edges;
    { // create the edges
        edges.reserve(polygon_points.size());
        for (size_t i = 0; i < polygon_points.size() - 1; i++) {
            edges.push_back(glm::ivec2(int(i), int(i + 1)));
        }

        // last edge between start and end vertex
        edges.push_back(glm::ivec2(polygon_points.size() - 1, 0));
    }

    // triangulation
    CDT::Triangulation<double> cdt;
    cdt.insertVertices(polygon_points.begin(), polygon_points.end(), [](const glm::vec2& p) { return p.x; }, [](const glm::vec2& p) { return p.y; });
    cdt.insertEdges(edges.begin(), edges.end(), [](const glm::ivec2& p) { return p.x; }, [](const glm::ivec2& p) { return p.y; });
    cdt.eraseOuterTrianglesAndHoles();

    // fill our own data structures
    std::transform(cdt.triangles.begin(), cdt.triangles.end(), std::back_inserter(triangle_collection.triangles), [this, style_index, cdt](CDT::Triangle tri) {
        return create_ordered_triangle(glm::vec2(cdt.vertices[tri.vertices[0]].x, cdt.vertices[tri.vertices[0]].y),
            glm::vec2(cdt.vertices[tri.vertices[1]].x, cdt.vertices[tri.vertices[1]].y),
            glm::vec2(cdt.vertices[tri.vertices[2]].x, cdt.vertices[tri.vertices[2]].y),
            style_index);
    });
}

// index_[a,b,c] are in counter-clockwise (CCW) winding order
// the resulting triangle struct disregards this order and assigns the indices to a top, middle, bottom variable
Triangle Preprocessor::create_ordered_triangle(glm::vec2 triangle_vertex_a, glm::vec2 triangle_vertex_b, glm::vec2 triangle_vertex_c, unsigned int style_index)
{
    Triangle triangle;
    triangle.data.style_index = style_index;

    { // determine y order of all thre triangle points
        // assume bottom/top triangle indices from only the first 2 triangle points
        if (triangle_vertex_a.y < triangle_vertex_b.y) {
            triangle.data.top_vertex = triangle_vertex_a;
            triangle.data.bottom_vertex = triangle_vertex_b;
        } else {
            triangle.data.top_vertex = triangle_vertex_b;
            triangle.data.bottom_vertex = triangle_vertex_a;
        }

        // determine if the actual order of all three triangle points
        if (triangle.data.top_vertex.y < triangle_vertex_c.y) {
            if (triangle.data.bottom_vertex.y > triangle_vertex_c.y) {
                // already ordered correctly
                triangle.data.middle_vertex = triangle_vertex_c;
            } else {
                // third point is actually at the bottom
                triangle.data.middle_vertex = triangle.data.bottom_vertex;
                triangle.data.bottom_vertex = triangle_vertex_c;
            }
        } else {
            // third point is actually at the top
            triangle.data.middle_vertex = triangle.data.top_vertex;
            triangle.data.top_vertex = triangle_vertex_c;
        }
    }

    return triangle;
}

void Preprocessor::create_lines(VectorLayerLineCollection& line_collection, const std::vector<glm::vec2> line_points, unsigned int style_index)
{
    line_collection.lines.reserve(line_points.size());
    for (size_t i = 0; i < line_points.size() - 1; i++) {
        line_collection.lines.push_back({ line_points[i], line_points[i + 1], style_index });
    }
}

std::pair<glm::vec2, int> Preprocessor::calculate_dda_steps(const glm::vec2 line)
{
    int steps = abs(line.x) > abs(line.y) ? abs(line.x) : abs(line.y);

    auto step_size = glm::vec2(line.x / float(steps), line.y / float(steps));

    return std::make_pair(step_size, steps);
}

// writes data_index to m_processed_tiles[tile_id].cell_to_data on a line between top-bottom point
// fill_direction is used to determine if:
//      case 0) it should write to m_x_values_per_y_step the first x position of a new y step
//      casle -1) and case 1) if the algorithm should fill everything between the current x and the x value stored in m_x_values_per_y_step
//      case -1 and case 1 is esentially used to create a filled triangle from a previously created line that had case 0
void Preprocessor::dda_line(VectorLayerGrid& grid, const glm::vec2 origin, const glm::vec2 line, const glm::vec2 thickness_normal, unsigned int data_index, int fill_direction, bool is_triangle)
{
    glm::vec2 thickness_step_size;
    int thickness_steps;
    std::tie(thickness_step_size, thickness_steps) = calculate_dda_steps(thickness_normal);
    thickness_step_size /= 2.0f; // TODO recheck if this is still needed -> without it we might be able to increase the performance significantly
    thickness_steps *= 2;

    glm::vec2 step_size;
    int steps;
    std::tie(step_size, steps) = calculate_dda_steps(line);
    // step_size /= 2.0f;
    // steps *= 2;

    // make sure that at least one step is applied
    thickness_steps = std::max(thickness_steps, 1);
    steps = std::max(steps, 1);

    glm::vec2 current_start_position = origin;

    if (!is_triangle) // visualizing a line means that we want to increase the thickness on both sides of the line (not just outside of a triangle)
        current_start_position -= thickness_normal / 2.0f;

    for (int i = 0; i < thickness_steps; i++) {

        glm::vec2 current_position = current_start_position;

        for (int j = 0; j < steps; j++) {
            write_to_cell(grid, current_position, data_index);

            // add a write step to x/y -> this prevents holes from thickness steps to appear by making the lines thicker
            if (step_size.x < 1)
                write_to_cell(grid, current_position + glm::vec2(1, 0) * glm::sign(step_size.x), data_index);
            if (step_size.y < 1)
                write_to_cell(grid, current_position + glm::vec2(0, 1) * glm::sign(step_size.y), data_index);

            if (i == 0 && is_triangle) // only apply for the inner most layer of a triangle
            {
                if (int(current_position.y + step_size.y) > current_position.y) { // next step would go to next y value
                    if (fill_direction == 0) {
                        // we have to save the x location for a specific y row
                        if (m_tile_bounds.contains(current_position))
                            m_x_values_per_y_step[int(current_position.y)] = int(current_position.x);
                    } else {
                        // we have to fill all values
                        // if()
                        for (int j = current_position.x; j != m_x_values_per_y_step[int(current_position.y)]; j += fill_direction) {
                            write_to_cell(grid, glm::vec2(j, current_position.y), data_index);
                        }
                    }
                }
            }
            current_position += step_size;
        }

        current_start_position += thickness_step_size;
    }
}

void Preprocessor::write_to_cell(VectorLayerGrid& grid, glm::vec2 current_position, unsigned int data_index)
{
    if (m_tile_bounds.contains(current_position)) {
        grid[int(current_position.x) + m_grid_size.x * int(current_position.y)].insert(data_index);
    }
}

void Preprocessor::dda_triangle(VectorLayerTriangleCollection& triangle_collection, unsigned int triangle_index, float thickness)
{
    auto triangle = triangle_collection.triangles[triangle_index];

    auto edge_top_bottom = triangle.data.bottom_vertex - triangle.data.top_vertex;
    auto edge_top_middle = triangle.data.middle_vertex - triangle.data.top_vertex;
    auto edge_middle_bottom = triangle.data.bottom_vertex - triangle.data.middle_vertex;

    auto normal_top_bottom = glm::normalize(glm::vec2(-edge_top_bottom.y, edge_top_bottom.x));
    auto normal_top_middle = glm::normalize(glm::vec2(-edge_top_middle.y, edge_top_middle.x));
    auto normal_middle_bottom = glm::normalize(glm::vec2(-edge_middle_bottom.y, edge_middle_bottom.x));

    { // swap normal direction if they are incorrect (pointing to center)
        glm::vec2 centroid = (triangle.data.top_vertex + triangle.data.middle_vertex + triangle.data.bottom_vertex) / 3.0f;
        if (glm::dot(normal_top_bottom, centroid - triangle.data.top_vertex) > 0) {
            normal_top_bottom *= -1;
        }
        if (glm::dot(normal_top_middle, centroid - triangle.data.top_vertex) > 0) {
            normal_top_middle *= -1;
        }
        if (glm::dot(normal_middle_bottom, centroid - triangle.data.middle_vertex) > 0) {
            normal_middle_bottom *= -1;
        }
    }

    // top bottom line
    dda_line(triangle_collection.cell_to_data, triangle.data.top_vertex, edge_top_bottom, normal_top_bottom * thickness, triangle_index, 0, true);

    // do we have to fill the triangle from right to left(-1) or from left to right(1)
    int fill_direction = (triangle.data.middle_vertex.x < triangle.data.bottom_vertex.x) ? 1 : -1;

    // // top middle line
    dda_line(triangle_collection.cell_to_data, triangle.data.top_vertex, edge_top_middle, normal_top_middle * thickness, triangle_index, fill_direction, true);
    // // middle bottom line
    dda_line(triangle_collection.cell_to_data, triangle.data.middle_vertex, edge_middle_bottom, normal_middle_bottom * thickness, triangle_index, fill_direction, true);

    // DEBUG draw lines with certain index
    // dda_line(triangle.data.top_vertex, edge_top_bottom, normal_top_bottom * 1.0f, 100, 0, true);
    // dda_line(triangle.data.top_vertex, edge_top_middle, normal_top_middle * 1.0f, 100, fill_direction, true);
    // dda_line(triangle.data.middle_vertex, edge_middle_bottom, normal_middle_bottom * 1.0f, 100, fill_direction, true);

    add_end_cap(triangle_collection.cell_to_data, triangle.data.top_vertex, triangle_index, thickness);
    add_end_cap(triangle_collection.cell_to_data, triangle.data.middle_vertex, triangle_index, thickness);
    add_end_cap(triangle_collection.cell_to_data, triangle.data.bottom_vertex, triangle_index, thickness);
}

void Preprocessor::add_end_cap(VectorLayerGrid& grid, const glm::vec2 position, unsigned int data_index, float thickness)
{
    float thickness_test = ceil(thickness);

    for (int i = -thickness_test; i < thickness_test; i++) {
        for (int j = -thickness_test; j < thickness_test; j++) {
            float distance = glm::length(glm::vec2(i, j));
            if (distance <= thickness)
                write_to_cell(grid, position + glm::vec2(i, j), data_index);
        }
    }
}

nucleus::Raster<uint8_t> Preprocessor::visualize_grid(const VectorLayerGrid& processed_grid)
{
    std::vector<uint8_t> output_grid;

    for (size_t i = 0; i < processed_grid.size(); ++i) {

        // grid is a singular uint32_t value that encodes the start index of the triangle list and the amount of triangles
        if (processed_grid[i].size() == 0) {
            output_grid.push_back(0);
        } else {
            output_grid.push_back(255u);
        }
    }

    return nucleus::Raster<uint8_t>(m_grid_size.x, std::move(output_grid));
}

} // namespace nucleus::vector_layer
