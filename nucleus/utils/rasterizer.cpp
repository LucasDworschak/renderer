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

#ifdef ALP_RASTERIZER_CDT
#include <CDT.h>
#endif

#ifdef ALP_RASTERIZER_EARCUT
#include <earcut.hpp>

// allow vec2 points for earcut input
namespace mapbox {
namespace util {

    template <>
    struct nth<0, glm::vec2> {
        inline static auto get(const glm::vec2& t) { return t.x; };
    };
    template <>
    struct nth<1, glm::vec2> {
        inline static auto get(const glm::vec2& t) { return t.y; };
    };

} // namespace util
} // namespace mapbox
#endif

namespace nucleus::utils::rasterizer {
#ifdef ALP_RASTERIZER_CDT
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

std::vector<glm::vec2> triangulize(std::vector<std::vector<glm::vec2>> polygon_points, bool remove_duplicate_vertices)
{
    std::vector<glm::vec2> processed_triangles;

    // cdt needs to create edges and combine all polygons into one single vector
    std::vector<glm::vec2> combined_polygon_points;
    std::vector<glm::ivec2> edges;
    for (size_t j = 0; j < polygon_points.size(); ++j) {
        auto current_edges = nucleus::utils::rasterizer::generate_neighbour_edges(polygon_points[j].size(), combined_polygon_points.size());
        edges.insert(edges.end(), current_edges.begin(), current_edges.end());
        combined_polygon_points.insert(combined_polygon_points.end(), polygon_points[j].begin(), polygon_points[j].end());
    }

    // triangulation
    CDT::Triangulation<float> cdt;

    if (remove_duplicate_vertices) {
        CDT::RemoveDuplicatesAndRemapEdges<float>(
            combined_polygon_points,
            [](const glm::vec2& p) { return p.x; },
            [](const glm::vec2& p) { return p.y; },
            edges.begin(),
            edges.end(),
            [](const glm::ivec2& p) { return p.x; },
            [](const glm::ivec2& p) { return p.y; },
            [](CDT::VertInd start, CDT::VertInd end) { return glm::ivec2 { start, end }; });
    }

    cdt.insertVertices(combined_polygon_points.begin(), combined_polygon_points.end(), [](const glm::vec2& p) { return p.x; }, [](const glm::vec2& p) { return p.y; });
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

#endif

#ifdef ALP_RASTERIZER_EARCUT

std::pair<uint32_t, uint32_t> get_split_index(uint32_t index, const std::vector<uint32_t>& polygon_sizes)
{
    // first index test different since we use the previous size in the for loop
    if (index < polygon_sizes[0]) {
        return { 0, index };
    }

    for (uint32_t i = 1; i < polygon_sizes.size(); i++) {
        if (index < polygon_sizes[i]) {
            return { i, index - polygon_sizes[i - 1] };
        }
    }

    // should not happen -> the index does not match a valid polygon point
    assert(false);
    return { 0, 0 };
}

std::vector<glm::vec2> triangulize(std::vector<std::vector<glm::vec2>> polygon_points, bool)
{
    std::vector<glm::vec2> processed_triangles;

    auto indices = mapbox::earcut<uint32_t>(polygon_points);

    processed_triangles.reserve(indices.size() * 3);

    // move all polygons sizes to single accumulated array -> so that we can match the index
    std::vector<uint32_t> polygon_sizes;
    uint32_t previous_size = 0;
    for (size_t i = 0; i < polygon_points.size(); ++i) {
        polygon_sizes.push_back(previous_size + polygon_points[i].size());
        previous_size += previous_size + polygon_points[i].size();
    }

    // fill our own data structures
    for (size_t i = 0; i < indices.size() / 3; ++i) {

        auto ind0 = get_split_index(indices[i * 3 + 0], polygon_sizes);
        auto ind1 = get_split_index(indices[i * 3 + 1], polygon_sizes);
        auto ind2 = get_split_index(indices[i * 3 + 2], polygon_sizes);

        const std::array<glm::vec2, 3> vertices = { polygon_points[ind0.first][ind0.second], polygon_points[ind1.first][ind1.second], polygon_points[ind2.first][ind2.second] };

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
#endif

} // namespace nucleus::utils::rasterizer
