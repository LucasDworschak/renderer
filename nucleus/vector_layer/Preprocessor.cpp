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
// - enlarge the triangle (and also clamp border)
// - only lines
// - create buffer that will be send to gpu
// - call this class

// tileable
// eventually this class will be accessed from transform_and_emit of the specific scheduler.
// this scheduler contains a list of new and deleted quads. It furthermore creates a new gpu_tile where we write our processed data to
// this gpu_tile is than emitted via signals where they are processed and saved for opengl
// so in theory we should be able to get one ID that we should process
// create m_processed_tiles for up to 9 IDs that are finished processing (current tile + up to 8 neighbours)
// and then let the gl class worry about combining the same tile
// having said all this, we still should worry about -> what if the tile is also deleted in this same scheduler call
// and what if new_quads also contains the neighbours
// both might not so bad and might be solved automatically (but with worse performance), if we do nothing but implement it like described above
// alternativeley what if we send the neighbourhood automatically to this function -> create neigbhourhood in scheduler (without deleted_quads)
// send the gpu_tile + neighbours here for filling, combine old gpu quad data with new one here
// we still have to worry about combining quad data on gl side, but it should still be more efficient if we precombine them here.

Preprocessor::Preprocessor()
    : m_x_values_per_y_step(m_grid_size.y * 3, 0)
{
}

// polygon describe the outer edge of a closed shape
// -> neighbouring vertices form an edge
// last vertex connects to first vertex
void Preprocessor::preprocess(const nucleus::tile::Id tile_id, const std::vector<std::vector<glm::vec2>> polygons, const std::vector<unsigned int> style_indices)
{
    m_processed_tile = VectorLayer();
    m_processed_tile.cell_to_triangle = std::vector<std::unordered_set<uint8_t>>(m_grid_size.x * m_grid_size.y, std::unordered_set<uint8_t>());

    // // define the up direction depending on tile scheme (do we have to add 1 or subtract 1 when looking at above tile)
    m_tile_up_direction = (tile_id.scheme == tile::Scheme::Tms) ? 1 : -1;
    m_tile_bounds = { glm::dvec2(0, 0), glm::dvec2 { 64, 64 } }; // TODO calculate with srs.tile_bounds method

    // create the triangles from polygons
    for (size_t i = 0; i < polygons.size(); ++i) {
        create_triangles(polygons[i], style_indices[i]);
    }

    // fill the grid from the created triangles
    for (size_t i = 0; i < m_processed_tile.triangles.size(); ++i) {
        dda_triangle(i);
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

    { // vertex normals
      // TODO calculate after we confirmed sdf triangle shape and then compare if this approach is calculated correctly
      // auto edge0 = vertices[triangle_index_b] - vertices[triangle_index_a];
      // auto edge1 = vertices[triangle_index_c] - vertices[triangle_index_b];
      // auto edge2 = vertices[triangle_index_a] - vertices[triangle_index_c];

        // auto edge_normal0 = glm::normalize(glm::vec2(-edge0.y, edge0.x));
        // auto edge_normal1 = glm::normalize(glm::vec2(-edge1.y, edge1.x));
        // auto edge_normal2 = glm::normalize(glm::vec2(-edge2.y, edge2.x));

        // auto vertex_normal0 = edge_normal0 + edge_normal2;
    }

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

// writes data_index to m_processed_tiles[tile_id].cell_to_triangle on a line between top-bottom point
// fill_direction is used to determine if:
//      case 0) it should write to m_x_values_per_y_step the first x position of a new y step
//      casle -1) and case 1) if the algorithm should fill everything between the current x and the x value stored in m_x_values_per_y_step
//      case -1 and case 1 is esentially used to create a filled triangle from a previously created line that had case 0
void Preprocessor::dda_line(const glm::vec2 top_point, const glm::vec2 bottom_point, unsigned int data_index, int fill_direction)
{
    std::cout << data_index;
    if (top_point == bottom_point)
        return; // both points on same cell --> TODO what to do (also we are currently only checking if they are the same point and not cell)

    glm::vec2 delta = bottom_point - top_point;
    // contains 1,0,-1 depending of the sign of the component (-> do we got left or right in each step)
    glm::vec2 signed_dir = { (delta.x > 0) ? 1 : (delta.x < 0) ? -1 : 0, (delta.y > 0) ? 1 : (delta.y < 0) ? -1 : 0 };

    // calculate the step size for x and y direction (one direction always moves in increments of 1)
    glm::vec2 step_size;
    int steps;
    if (abs(delta.x) > abs(delta.y)) {
        // we increment x by 1/-1 each step and increment y by a fracture
        steps = abs(delta.x);
        step_size = { signed_dir.x, signed_dir.y * (abs(delta.y) / abs(delta.x)) };
    } else {
        // we increment y by 1/-1 each step and increment x by a fracture
        steps = abs(delta.y);
        step_size = { signed_dir.x * abs(delta.x) / abs(delta.y), signed_dir.y };
    }

    // first point
    glm::vec2 currentPoint = top_point;
    glm::vec2 lastPoint = top_point;
    // m_processed_tiles[tile_id].cell_to_triangle[int(currentPoint.x) + m_grid_size.x * int(currentPoint.y)].insert(data_index);
    // write_to_cell(currentPoint, data_index);

    int last_y = currentPoint.y;
    if (fill_direction == 0 && last_y >= 0)
        m_x_values_per_y_step[int(currentPoint.y)] = int(currentPoint.x);

    // draw the points
    for (int i = 0; i < steps; i++) {
        lastPoint = currentPoint;
        currentPoint += step_size;

        // m_processed_tiles[tile_id].cell_to_triangle[int(currentPoint.x) + m_grid_size.x * int(currentPoint.y)].insert(data_index);
        // write_to_cell(currentPoint, data_index);

        if (last_y != int(currentPoint.y)) {
            last_y = currentPoint.y;
            if (last_y < 0)
                continue; // we are outside of current cell -> we dont need to remember

            if (fill_direction == 0)
                m_x_values_per_y_step[currentPoint.y] = currentPoint.x;
            else {
                for (int j = lastPoint.x; j != m_x_values_per_y_step[int(lastPoint.y)]; j += fill_direction) {
                    // write_to_cell(glm::vec2(j, lastPoint.y), data_index);
                }
            }
        }
    }
}

void Preprocessor::write_to_cell(glm::vec2 current_position, std::array<glm::vec2, 3> points, unsigned int data_index, float distance)
{

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

void Preprocessor::dda_triangle(unsigned int triangle_index)
{
    auto vertices = m_processed_tile.vertices;
    auto triangle = m_processed_tile.triangles[triangle_index];

    auto top_vertice = vertices[triangle.top_index];
    auto middle_vertice = vertices[triangle.middle_index];
    auto bottom_vertice = vertices[triangle.bottom_index];

    float dist_to_grid_center = m_tile_bounds.size().x / m_grid_size.x / sqrt(2);
    float distance = 5 + dist_to_grid_center;

    for (uint8_t i = 0; i < m_grid_size.x; i++) {
        for (uint8_t j = 0; j < m_grid_size.y; j++) {
            auto current_position = glm::vec2 { i, j };

            write_to_cell(current_position, { top_vertice, middle_vertice, bottom_vertice }, triangle_index, distance);

            // what can we abstract:
            // edge construction
            // sign s
            // dot(e,e)
        }
    }

    // scale triangle
    {

        // TODO implement

        // move edges along normal by distance
        // edge must be longer
        // in the dda function calculate the distance to the original vertices and calculate if it is in the cell
        //  -> produces rounded triangles
        // edges as origin + distance

        // TODO currently the scaling is done by percentage
        // this means that larger triangles are scaled more than smaller triangles
        // shouldn't we scale it by distance instead? -> so that every triangle is scaled up by xxx meters?
        // solution would probably be something like: calculate the size of the triangle and create an appropriate scale factor from this
        // auto scale_factor = glm::vec2(2.2);
        // auto centroid = (top_vertice + middle_vertice + bottom_vertice) / glm::vec2(3.0);

        // top_vertice = ((top_vertice - centroid) * scale_factor) + centroid;
        // middle_vertice = ((middle_vertice - centroid) * scale_factor) + centroid;
        // bottom_vertice = ((bottom_vertice - centroid) * scale_factor) + centroid;
    }

    // // top bottom line
    // dda_line(top_vertice, bottom_vertice, triangle_index, 0);

    // // do we have to fill the triangle from right to left(-1) or from left to right(1)
    // int fill_direction = (middle_vertice.x < bottom_vertice.x) ? 1 : -1;

    // // // top middle line
    // dda_line(top_vertice, middle_vertice, triangle_index, fill_direction);
    // // // middle bottom line
    // dda_line(middle_vertice, bottom_vertice, triangle_index, fill_direction);
}

void Preprocessor::visualize_grid()
{
    auto raster = nucleus::Raster<uint8_t>(m_grid_size, 0);

    for (size_t i = 0; i < m_processed_tile.cell_to_triangle.size(); ++i) {
        if (m_processed_tile.cell_to_triangle[i].size() > 0) {
            raster.pixel({ i % m_grid_size.x, i / m_grid_size.x }) = 255;
        }
    }

    const auto debug_out = QImage(raster.bytes(), m_grid_size.x, m_grid_size.y, QImage::Format_Grayscale8);
    debug_out.save(QString("grid_test.png"));
}

} // namespace nucleus::vector_layer
