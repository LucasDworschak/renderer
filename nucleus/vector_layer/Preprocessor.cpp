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

namespace nucleus::vector_layer {

// TODO current problems:
// - create buffer that will be send to gpu
// - call this class

Preprocessor::Preprocessor()
    : m_x_values_per_y_step(m_grid_size.y * 3, 0)
{
}

// polygon describe the outer edge of a closed shape
// -> neighbouring vertices form an edge
// last vertex connects to first vertex
void Preprocessor::preprocess_triangles(const nucleus::tile::Id tile_id, const std::vector<std::vector<glm::vec2>> polygons, const std::vector<unsigned int> style_indices)
{
    m_processed_tile = VectorLayer();
    m_processed_tile.cell_to_triangle = std::vector<std::unordered_set<uint8_t>>(m_grid_size.x * m_grid_size.y, std::unordered_set<uint8_t>());

    // // define the up direction depending on tile scheme (do we have to add 1 or subtract 1 when looking at above tile)
    m_tile_up_direction = (tile_id.scheme == tile::Scheme::Tms) ? 1 : -1;
    m_tile_bounds = { glm::dvec2(0, 0), glm::dvec2 { 64, 64 } }; // TODO calculate with srs.tile_bounds method

    float dist_to_grid_center = m_tile_bounds.size().x / m_grid_size.x / sqrt(2);
    float thickness = 5.0f + dist_to_grid_center;

    // create the triangles from polygons
    for (size_t i = 0; i < polygons.size(); ++i) {
        create_triangles(polygons[i], style_indices[i]);
    }

    // fill the grid from the created triangles
    for (size_t i = 0; i < m_processed_tile.triangles.size(); ++i) {
        dda_triangle(i, thickness);
    }
}

void Preprocessor::preprocess_lines(const nucleus::tile::Id tile_id, const std::vector<std::vector<glm::vec2>> polygons, const std::vector<unsigned int> style_indices)
{
    m_processed_tile = VectorLayer();
    m_processed_tile.cell_to_triangle = std::vector<std::unordered_set<uint8_t>>(m_grid_size.x * m_grid_size.y, std::unordered_set<uint8_t>());

    // // define the up direction depending on tile scheme (do we have to add 1 or subtract 1 when looking at above tile)
    m_tile_up_direction = (tile_id.scheme == tile::Scheme::Tms) ? 1 : -1;
    m_tile_bounds = { glm::dvec2(0, 0), glm::dvec2 { 64, 64 } }; // TODO calculate with srs.tile_bounds method

    float dist_to_grid_center = m_tile_bounds.size().x / m_grid_size.x / sqrt(2);
    float thickness = 5.0f + dist_to_grid_center;

    // create the triangles from polygons
    for (size_t i = 0; i < polygons.size(); ++i) {
        create_lines(polygons[i], style_indices[i]);
    }

    // fill the grid from the created triangles
    for (size_t i = 0; i < m_processed_tile.lines.size(); ++i) {
        auto vertex_start = m_processed_tile.vertices[m_processed_tile.lines[i].line_start_index];
        auto vertex_end = m_processed_tile.vertices[m_processed_tile.lines[i].line_end_index];
        glm::vec2 edge = vertex_end - vertex_start;
        glm::vec2 normal = glm::normalize(glm::vec2(-edge.y, edge.x));

        dda_line(vertex_start, edge, normal * thickness, m_processed_tile.lines[i].style_index, 0, false);
        // dda_line(vertex_start, edge, normal * 1.0f, 100, 0, false); // DEBUG draw only line without thickness

        add_end_cap(vertex_start, m_processed_tile.lines[i].style_index, thickness / 2.0);
        add_end_cap(vertex_end, m_processed_tile.lines[i].style_index, thickness / 2.0);
    }
}

void Preprocessor::create_triangles(const std::vector<glm::vec2> polygon_points, unsigned int style_index)
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

    // if other polygons have already been parsed, we have to offset the triangles to consider the correct vertices
    unsigned int index_offset = m_processed_tile.vertices.size();

    // fill our own vertice/triangle data structures
    std::transform(cdt.vertices.begin(), cdt.vertices.end(), std::back_inserter(m_processed_tile.vertices), [](CDT::V2d<double> vert) { return glm::vec2 { vert.x, vert.y }; });

    std::transform(cdt.triangles.begin(), cdt.triangles.end(), std::back_inserter(m_processed_tile.triangles), [this, style_index, index_offset](CDT::Triangle tri) {
        return create_ordered_triangle(tri.vertices[0] + index_offset, tri.vertices[1] + index_offset, tri.vertices[2] + index_offset, style_index);
    });
}

// index_[a,b,c] are in counter-clockwise (CCW) winding order
// the resulting triangle struct disregards this order and assigns the indices to a top, middle, bottom variable
Triangle Preprocessor::create_ordered_triangle(unsigned int triangle_index_a, unsigned int triangle_index_b, unsigned int triangle_index_c, unsigned int style_index)
{
    auto vertices = m_processed_tile.vertices;

    Triangle triangle_data;
    triangle_data.style_index = style_index;

    { // determine y order of all thre triangle points
        // assume bottom/top triangle indices from only the first 2 triangle points
        if (vertices[triangle_index_a].y < vertices[triangle_index_b].y) {
            triangle_data.top_index = triangle_index_a;
            triangle_data.bottom_index = triangle_index_b;
        } else {
            triangle_data.top_index = triangle_index_b;
            triangle_data.bottom_index = triangle_index_a;
        }

        // determine if the actual order of all three triangle points
        if (vertices[triangle_data.top_index].y < vertices[triangle_index_c].y) {
            if (vertices[triangle_data.bottom_index].y > vertices[triangle_index_c].y) {
                // already ordered correctly
                triangle_data.middle_index = triangle_index_c;
            } else {
                // third point is actually at the bottom
                triangle_data.middle_index = triangle_data.bottom_index;
                triangle_data.bottom_index = triangle_index_c;
            }
        } else {
            // third point is actually at the top
            triangle_data.middle_index = triangle_data.top_index;
            triangle_data.top_index = triangle_index_c;
        }
    }

    return triangle_data;
}

void Preprocessor::create_lines(const std::vector<glm::vec2> line_points, unsigned int style_index)
{
    unsigned int index_offset = m_processed_tile.vertices.size();

    m_processed_tile.lines.reserve(line_points.size());
    for (size_t i = 0; i < line_points.size() - 1; i++) {
        m_processed_tile.lines.push_back({ index_offset + (unsigned int)(i), index_offset + (unsigned int)(i + 1), style_index });
    }

    m_processed_tile.vertices.insert(m_processed_tile.vertices.end(), line_points.begin(), line_points.end());
}

std::pair<glm::vec2, int> Preprocessor::calculate_dda_steps(const glm::vec2 line)
{
    int steps = abs(line.x) > abs(line.y) ? abs(line.x) : abs(line.y);

    auto step_size = glm::vec2(line.x / float(steps), line.y / float(steps));

    return std::make_pair(step_size, steps);
}

// writes data_index to m_processed_tiles[tile_id].cell_to_triangle on a line between top-bottom point
// fill_direction is used to determine if:
//      case 0) it should write to m_x_values_per_y_step the first x position of a new y step
//      casle -1) and case 1) if the algorithm should fill everything between the current x and the x value stored in m_x_values_per_y_step
//      case -1 and case 1 is esentially used to create a filled triangle from a previously created line that had case 0
void Preprocessor::dda_line(const glm::vec2 origin, const glm::vec2 line, const glm::vec2 thickness_normal, unsigned int data_index, int fill_direction, bool is_triangle)
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

    glm::vec2 currentStartPoint = origin;

    if (!is_triangle) // visualizing a line means that we want to increase the thickness on both sides of the line (not just outside of a triangle)
        currentStartPoint -= thickness_normal / 2.0f;

    for (int i = 0; i < thickness_steps; i++) {

        glm::vec2 currentPoint = currentStartPoint;

        for (int j = 0; j < steps; j++) {
            write_to_cell(currentPoint, data_index);

            // add a write step to x/y -> this prevents holes from thickness steps to appear by making the lines thicker
            if (step_size.x < 1)
                write_to_cell(currentPoint + glm::vec2(1, 0) * glm::sign(step_size.x), data_index);
            if (step_size.y < 1)
                write_to_cell(currentPoint + glm::vec2(0, 1) * glm::sign(step_size.y), data_index);

            if (i == 0 && is_triangle) // only apply for the inner most layer of a triangle
            {
                if (int(currentPoint.y + step_size.y) > currentPoint.y) { // next step would go to next y value
                    if (fill_direction == 0) {
                        // we have to save the x location for a specific y row
                        m_x_values_per_y_step[int(currentPoint.y)] = int(currentPoint.x);
                    } else {
                        // we have to fill all values
                        // if()
                        for (int j = currentPoint.x; j != m_x_values_per_y_step[int(currentPoint.y)]; j += fill_direction) {
                            write_to_cell(glm::vec2(j, currentPoint.y), data_index);
                        }
                    }
                }
            }
            currentPoint += step_size;
        }

        currentStartPoint += thickness_step_size;
    }
}

void Preprocessor::write_to_cell(glm::vec2 current_position, unsigned int data_index)
{
    if (m_tile_bounds.contains(current_position)) {
        m_processed_tile.cell_to_triangle[int(current_position.x) + m_grid_size.x * int(current_position.y)].insert(data_index);
    }
}

void Preprocessor::write_to_cell_sdf(glm::vec2 current_position, std::array<glm::vec2, 3> points, unsigned int data_index, float distance)
{

    // what can we abstract:
    // edge construction
    // sign s
    // dot(e,e)

    glm::vec2 e0 = points[1] - points[0], e1 = points[2] - points[1], e2 = points[0] - points[2];
    glm::vec2 v0 = current_position - points[0], v1 = current_position - points[1], v2 = current_position - points[2];
    glm::vec2 pq0 = v0 - e0 * glm::clamp(glm::dot(v0, e0) / glm::dot(e0, e0), 0.0f, 1.0f);
    glm::vec2 pq1 = v1 - e1 * glm::clamp(glm::dot(v1, e1) / glm::dot(e1, e1), 0.0f, 1.0f);
    glm::vec2 pq2 = v2 - e2 * glm::clamp(glm::dot(v2, e2) / glm::dot(e2, e2), 0.0f, 1.0f);
    float s = glm::sign(e0.x * e2.y - e0.y * e2.x);
    glm::vec2 d
        = min(min(glm::vec2(dot(pq0, pq0), s * (v0.x * e0.y - v0.y * e0.x)), glm::vec2(dot(pq1, pq1), s * (v1.x * e1.y - v1.y * e1.x))), glm::vec2(dot(pq2, pq2), s * (v2.x * e2.y - v2.y * e2.x)));
    if (-sqrt(d.x) * glm::sign(d.y) < distance)
        m_processed_tile.cell_to_triangle[int(current_position.x) + m_grid_size.x * int(current_position.y)].insert(data_index);
}

// also we would probably need a dda for the thickness -> but we can probably only calculated this once at the start? and use those few points as offsets in every step
void Preprocessor::dda_triangle(unsigned int triangle_index, float thickness)
{
    auto vertices = m_processed_tile.vertices;
    auto triangle = m_processed_tile.triangles[triangle_index];

    auto top_vertice = vertices[triangle.top_index];
    auto middle_vertice = vertices[triangle.middle_index];
    auto bottom_vertice = vertices[triangle.bottom_index];

    auto edge_top_bottom = bottom_vertice - top_vertice;
    auto edge_top_middle = middle_vertice - top_vertice;
    auto edge_middle_bottom = bottom_vertice - middle_vertice;

    // TODO normals are not correct -> some might point inwards --> we need to calculate them in create_ordered_triangle
    auto normal_top_bottom = glm::normalize(glm::vec2(-edge_top_bottom.y, edge_top_bottom.x));
    auto normal_top_middle = glm::normalize(glm::vec2(-edge_top_middle.y, edge_top_middle.x));
    auto normal_middle_bottom = glm::normalize(glm::vec2(-edge_middle_bottom.y, edge_middle_bottom.x));

    { // swap normal direction if they are incorrect (pointing to center)
        glm::vec2 centroid = (top_vertice + middle_vertice + bottom_vertice) / 3.0f;
        if (glm::dot(normal_top_bottom, centroid - top_vertice) > 0) {
            normal_top_bottom *= -1;
        }
        if (glm::dot(normal_top_middle, centroid - top_vertice) > 0) {
            normal_top_middle *= -1;
        }
        if (glm::dot(normal_middle_bottom, centroid - middle_vertice) > 0) {
            normal_middle_bottom *= -1;
        }
    }

    // top bottom line
    dda_line(top_vertice, edge_top_bottom, normal_top_bottom * thickness, triangle_index, 0, true);

    // do we have to fill the triangle from right to left(-1) or from left to right(1)
    int fill_direction = (middle_vertice.x < bottom_vertice.x) ? 1 : -1;

    // // top middle line
    dda_line(top_vertice, edge_top_middle, normal_top_middle * thickness, triangle_index, fill_direction, true);
    // // middle bottom line
    dda_line(middle_vertice, edge_middle_bottom, normal_middle_bottom * thickness, triangle_index, fill_direction, true);

    // DEBUG draw lines with certain index
    // dda_line(top_vertice, edge_top_bottom, normal_top_bottom * 1.0f, 100, 0, true);
    // dda_line(top_vertice, edge_top_middle, normal_top_middle * 1.0f, 100, fill_direction, true);
    // dda_line(middle_vertice, edge_middle_bottom, normal_middle_bottom * 1.0f, 100, fill_direction, true);

    add_end_cap(top_vertice, triangle_index, thickness);
    add_end_cap(middle_vertice, triangle_index, thickness);
    add_end_cap(bottom_vertice, triangle_index, thickness);
}

void Preprocessor::sdf_triangle(unsigned int triangle_index)
{
    auto vertices = m_processed_tile.vertices;
    auto triangle = m_processed_tile.triangles[triangle_index];

    auto top_vertice = vertices[triangle.top_index];
    auto middle_vertice = vertices[triangle.middle_index];
    auto bottom_vertice = vertices[triangle.bottom_index];

    float dist_to_grid_center = m_tile_bounds.size().x / m_grid_size.x / sqrt(2);
    float distance = 0 + dist_to_grid_center;

    for (uint8_t i = 0; i < m_grid_size.x; i++) {
        for (uint8_t j = 0; j < m_grid_size.y; j++) {
            auto current_position = glm::vec2 { i, j };

            write_to_cell_sdf(current_position, { top_vertice, middle_vertice, bottom_vertice }, triangle_index, distance);
        }
    }
}

void Preprocessor::add_end_cap(const glm::vec2 position, unsigned int data_index, float thickness)
{
    float thickness_test = ceil(thickness);

    for (int i = -thickness_test; i < thickness_test; i++) {
        for (int j = -thickness_test; j < thickness_test; j++) {
            float distance = glm::length(glm::vec2(i, j));
            if (distance <= thickness)
                write_to_cell(position + glm::vec2(i, j), data_index);
        }
    }
}

void Preprocessor::visualize_grid()
{
    auto raster = nucleus::Raster<uint8_t>(m_grid_size, 0);

    for (size_t i = 0; i < m_processed_tile.cell_to_triangle.size(); ++i) {
        if (m_processed_tile.cell_to_triangle[i].size() > 0) {
            // raster.pixel({ i % m_grid_size.x, i / m_grid_size.x }) = (m_processed_tile.cell_to_triangle[i].size() / 4.0) * 255;

            if (m_processed_tile.cell_to_triangle[i].contains(100))
                raster.pixel({ i % m_grid_size.x, i / m_grid_size.x }) = 100;
            else
                raster.pixel({ i % m_grid_size.x, i / m_grid_size.x }) = 255;
        }
    }

    const auto debug_out = QImage(raster.bytes(), m_grid_size.x, m_grid_size.y, QImage::Format_Grayscale8);
    debug_out.save(QString("grid_test.png"));
}

} // namespace nucleus::vector_layer
