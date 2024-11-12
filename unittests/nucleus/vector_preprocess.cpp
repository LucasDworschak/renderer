/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2024 Lucas Dworschak
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

#include <catch2/catch_test_macros.hpp>
#include <glm/glm.hpp>

#include <QImage>
#include <QString>

#include <CDT.h>

#include <radix/tile.h>

#include "nucleus/Raster.h"
#include "nucleus/vector_layer/Preprocessor.h"
#include "nucleus/vector_tile/types.h"

#include <iostream>

using namespace nucleus::vector_tile;
using namespace nucleus::vector_layer;

TEST_CASE("nucleus/vector_preprocess")
{
    SECTION("Triangulation")
    {
        // 5 point polygon
        // basically a square where one side contains an inward facing triangle
        // a triangulation algorithm should be able to discern that 3 triangles are needed to construct this shape
        const std::vector<glm::vec2> points = { glm::vec2(0, 0), glm::vec2(1, 1), glm::vec2(0, 2), glm::vec2(2, 2), glm::vec2(2, 0) };
        const std::vector<glm::ivec2> edges = { glm::ivec2(0, 1), glm::ivec2(1, 2), glm::ivec2(2, 3), glm::ivec2(3, 4), glm::ivec2(4, 0) };

        CDT::Triangulation<double> cdt;
        cdt.insertVertices(points.begin(), points.end(), [](const glm::vec2& p) { return p.x; }, [](const glm::vec2& p) { return p.y; });
        cdt.insertEdges(edges.begin(), edges.end(), [](const glm::ivec2& p) { return p.x; }, [](const glm::ivec2& p) { return p.y; });
        cdt.eraseOuterTrianglesAndHoles();

        auto tri = cdt.triangles;
        auto vert = cdt.vertices;

        // check if only 3 triangles have been found
        CHECK(tri.size() == 3);

        // 1st triangle
        CHECK(vert[tri[0].vertices[0]].x == 0.0);
        CHECK(vert[tri[0].vertices[0]].y == 2.0);
        CHECK(vert[tri[0].vertices[1]].x == 1.0);
        CHECK(vert[tri[0].vertices[1]].y == 1.0);
        CHECK(vert[tri[0].vertices[2]].x == 2.0);
        CHECK(vert[tri[0].vertices[2]].y == 2.0);

        // 2nd triangle
        CHECK(vert[tri[1].vertices[0]].x == 1.0);
        CHECK(vert[tri[1].vertices[0]].y == 1.0);
        CHECK(vert[tri[1].vertices[1]].x == 2.0);
        CHECK(vert[tri[1].vertices[1]].y == 0.0);
        CHECK(vert[tri[1].vertices[2]].x == 2.0);
        CHECK(vert[tri[1].vertices[2]].y == 2.0);

        // 3rd triangle
        CHECK(vert[tri[2].vertices[0]].x == 1.0);
        CHECK(vert[tri[2].vertices[0]].y == 1.0);
        CHECK(vert[tri[2].vertices[1]].x == 0.0);
        CHECK(vert[tri[2].vertices[1]].y == 0.0);
        CHECK(vert[tri[2].vertices[2]].x == 2.0);
        CHECK(vert[tri[2].vertices[2]].y == 0.0);

        // print out all the points of the triangles (to check what might have went wrong)
        // for (std::size_t i = 0; i < tri.size(); i++) {
        //     printf("Triangle points: [[%f, %f], [%f, %f], [%f, %f]]\n",
        //         vert[tri[i].vertices[0]].x, // x0
        //         vert[tri[i].vertices[0]].y, // y0
        //         vert[tri[i].vertices[1]].x, // x1
        //         vert[tri[i].vertices[1]].y, // y1
        //         vert[tri[i].vertices[2]].x, // x2
        //         vert[tri[i].vertices[2]].y // y2
        //     );
        // }
    }

    // SECTION("Triangle y ordering")
    // {
    //     // make sure that the triangle_order function correctly orders the triangle points from lowest y to highest y value
    //     const std::vector<glm::vec2> triangle_points_012 = { { glm::vec2(30, 10), glm::vec2(10, 30), glm::vec2(50, 50) } };
    //     const std::vector<glm::vec2> triangle_points_021 = { glm::vec2(30, 10), glm::vec2(50, 50), glm::vec2(10, 30) };
    //     const std::vector<glm::vec2> triangle_points_102 = { glm::vec2(10, 30), glm::vec2(30, 10), glm::vec2(50, 50) };
    //     const std::vector<glm::vec2> triangle_points_201 = { glm::vec2(10, 30), glm::vec2(50, 50), glm::vec2(30, 10) };
    //     const std::vector<glm::vec2> triangle_points_120 = { glm::vec2(50, 50), glm::vec2(30, 10), glm::vec2(10, 30) };
    //     const std::vector<glm::vec2> triangle_points_210 = { glm::vec2(50, 50), glm::vec2(10, 30), glm::vec2(30, 10) };

    //     const std::vector<std::vector<glm::vec2>> triangle_points = { triangle_points_012, triangle_points_021, triangle_points_102, triangle_points_201, triangle_points_120, triangle_points_210 };

    //     const std::vector<unsigned int> style_indices = { 1, 2, 3, 4, 5, 6 };

    //     Preprocessor p;
    //     p.preprocess(triangle_points, style_indices);
    //     auto t = p.triangles();
    //     CHECK(t.size() == 6);

    //     // std::cout << t[0] << std::endl;
    //     // std::cout << t[1] << std::endl;
    //     // std::cout << t[2] << std::endl;
    //     // std::cout << t[3] << std::endl;
    //     // std::cout << t[4] << std::endl;
    //     // std::cout << t[5] << std::endl;

    //     CHECK(t[0] == Triangle { 0, 1, 2, 1 });
    //     CHECK(t[1] == Triangle { 3, 5, 4, 2 });
    //     CHECK(t[2] == Triangle { 7, 6, 8, 3 });
    //     CHECK(t[3] == Triangle { 11, 9, 10, 4 });
    //     CHECK(t[4] == Triangle { 13, 14, 12, 5 });
    //     CHECK(t[5] == Triangle { 17, 16, 15, 6 });
    // }

    SECTION("Triangle to Grid")
    {
        // we have a basic triangle and want to add an entry on a grid where the triangle overlaps

        // auto raster = nucleus::Raster<uint8_t>({ 64, 64 }, 0);

        const std::vector<glm::vec2> triangle_left_hypo = { glm::vec2(10, 30), glm::vec2(30, 10), glm::vec2(50, 50) };
        const std::vector<glm::vec2> triangle_right_hypo = { glm::vec2(5, 5), glm::vec2(15, 10), glm::vec2(5, 15) };
        // const std::vector<glm::vec2> triangle_test = { glm::vec2(25, 25), glm::vec2(30, 25), glm::vec2(28, 30) };

        const std::vector<std::vector<glm::vec2>> triangle_points = { triangle_left_hypo, triangle_right_hypo };
        const std::vector<std::vector<glm::vec2>> line_points
            = { { glm::vec2(10.5, 40.5), glm::vec2(30.5, 20.5) }, { glm::vec2(10.5, 5.5), glm::vec2(30.5, 5.5) }, { glm::vec2(10.5, 50), glm::vec2(30.5, 50) }

              };
        const std::vector<unsigned int> style_indices = { 1, 2, 3 };

        const auto id = nucleus::tile::Id { .zoom_level = 10, .coords = { 548, 359 }, .scheme = nucleus::tile::Scheme::SlippyMap };

        Preprocessor p;
        p.preprocess_triangles(id, triangle_points, style_indices);
        // p.preprocess_lines(id, line_points, style_indices);
        p.visualize_grid();
    }
}
