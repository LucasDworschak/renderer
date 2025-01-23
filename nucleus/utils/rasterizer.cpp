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

#include "rasterizer.h"

#include <CDT.h>

namespace nucleus::utils::rasterizer {

std::vector<glm::ivec2> generate_neighbour_edges(size_t num_points, const size_t start_offset)
{
    std::vector<glm::ivec2> edges;
    { // create the edges
        edges.reserve(num_points);
        for (size_t i = 0; i < num_points - 1; i++) {
            edges.push_back(glm::ivec2(start_offset + int(i), start_offset + int(i + 1)));
        }

        // last edge between start and end vertex
        edges.push_back(glm::ivec2(start_offset + num_points - 1, start_offset));
    }

    return edges;
}

std::vector<glm::vec2> triangulize(std::vector<glm::vec2> polygon_points, std::vector<glm::ivec2> edges, bool remove_duplicate_vertices)
{
    std::vector<glm::vec2> processed_triangles;

    // triangulation
    CDT::Triangulation<float> cdt;

    if (remove_duplicate_vertices) {
        CDT::RemoveDuplicatesAndRemapEdges<float>(
            polygon_points,
            [](const glm::vec2& p) { return p.x; },
            [](const glm::vec2& p) { return p.y; },
            edges.begin(),
            edges.end(),
            [](const glm::ivec2& p) { return p.x; },
            [](const glm::ivec2& p) { return p.y; },
            [](CDT::VertInd start, CDT::VertInd end) { return glm::ivec2 { start, end }; });
    }

    cdt.insertVertices(polygon_points.begin(), polygon_points.end(), [](const glm::vec2& p) { return p.x; }, [](const glm::vec2& p) { return p.y; });
    cdt.insertEdges(edges.begin(), edges.end(), [](const glm::ivec2& p) { return p.x; }, [](const glm::ivec2& p) { return p.y; });
    cdt.eraseOuterTrianglesAndHoles();

    processed_triangles.reserve(cdt.triangles.size() * 3);

    // fill our own data structures
    for (size_t i = 0; i < cdt.triangles.size(); ++i) {
        auto tri = cdt.triangles[i];

        const std::array<CDT::V2d<float>, 3> vertices = { cdt.vertices[tri.vertices[0]], cdt.vertices[tri.vertices[1]], cdt.vertices[tri.vertices[2]] };

        const int top_index = (vertices[0].y < vertices[1].y) ? ((vertices[0].y < vertices[2].y) ? 0 : 2) : ((vertices[1].y < vertices[2].y) ? 1 : 2);
        // for middle and bottom index we first initialize them randomly with the values that still need to be tested
        int middle_index;
        int bottom_index;
        if (top_index == 0) {
            middle_index = 1;
            bottom_index = 2;
        } else if (top_index == 1) {
            middle_index = 2;
            bottom_index = 0;
        } else {
            middle_index = 0;
            bottom_index = 1;
        }

        // and now we test if we assigned them correctly
        if (vertices[middle_index].y > vertices[bottom_index].y) {
            // if not we have to interchange them
            int tmp = middle_index;
            middle_index = bottom_index;
            bottom_index = tmp;
        }

        // lastly add the vertices to the vector in the correct order
        processed_triangles.push_back({ vertices[top_index].x, vertices[top_index].y });
        processed_triangles.push_back({ vertices[middle_index].x, vertices[middle_index].y });
        processed_triangles.push_back({ vertices[bottom_index].x, vertices[bottom_index].y });
    }

    return processed_triangles;
}

} // namespace nucleus::utils::rasterizer
