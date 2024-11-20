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

#include <QFile>
#include <QImage>
#include <QString>

#include <CDT.h>

#include <nucleus/tile/conversion.h>
#include <nucleus/utils/bit_coding.h>
#include <radix/tile.h>

#include "nucleus/Raster.h"
#include "nucleus/vector_layer/Preprocessor.h"
#include "nucleus/vector_tile/types.h"

using namespace nucleus::vector_tile;
using namespace nucleus::vector_layer;

inline std::ostream& operator<<(std::ostream& os, const glm::uvec3& v) { return os << "{ " << v.x << ", " << v.y << ", " << v.z << " }"; }

QImage example_grid_data_triangles()
{
    auto file = QFile(QString("%1%2").arg(ALP_TEST_DATA_DIR, "vector_layer_grid_triangles.png"));
    file.open(QFile::ReadOnly);
    const auto bytes = file.readAll();
    auto image = QImage::fromData(bytes);
    REQUIRE(!image.isNull());
    return image;
}

QImage example_grid_data_lines()
{
    auto file = QFile(QString("%1%2").arg(ALP_TEST_DATA_DIR, "vector_layer_grid_lines.png"));
    file.open(QFile::ReadOnly);
    const auto bytes = file.readAll();
    auto image = QImage::fromData(bytes);
    REQUIRE(!image.isNull());
    return image;
}

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

        // DEBUG print out all the points of the triangles (to check what might have went wrong)
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

    SECTION("Triangle y ordering")
    {
        // make sure that the triangle_order function correctly orders the triangle points from lowest y to highest y value
        const std::vector<glm::vec2> triangle_points_012 = { { glm::vec2(30, 10), glm::vec2(10, 30), glm::vec2(50, 50) } };
        const std::vector<glm::vec2> triangle_points_021 = { glm::vec2(30, 10), glm::vec2(50, 50), glm::vec2(10, 30) };
        const std::vector<glm::vec2> triangle_points_102 = { glm::vec2(10, 30), glm::vec2(30, 10), glm::vec2(50, 50) };
        const std::vector<glm::vec2> triangle_points_201 = { glm::vec2(10, 30), glm::vec2(50, 50), glm::vec2(30, 10) };
        const std::vector<glm::vec2> triangle_points_120 = { glm::vec2(50, 50), glm::vec2(30, 10), glm::vec2(10, 30) };
        const std::vector<glm::vec2> triangle_points_210 = { glm::vec2(50, 50), glm::vec2(10, 30), glm::vec2(30, 10) };

        const std::vector<std::vector<glm::vec2>> triangle_points = { triangle_points_012, triangle_points_021, triangle_points_102, triangle_points_201, triangle_points_120, triangle_points_210 };

        const std::vector<unsigned int> style_indices = { 1, 1, 1, 1, 1, 1 };

        const auto id = nucleus::tile::Id { .zoom_level = 10, .coords = { 548, 359 }, .scheme = nucleus::tile::Scheme::SlippyMap };

        Preprocessor p(id);
        auto processed = p.preprocess_triangles(triangle_points, style_indices);
        CHECK(processed.triangles.size() == 6);

        Triangle correct(glm::vec2(30, 10), glm::vec2(10, 30), glm::vec2(50, 50), 1);

        CHECK(processed.triangles[0] == correct);
        CHECK(processed.triangles[1] == correct);
        CHECK(processed.triangles[2] == correct);
        CHECK(processed.triangles[3] == correct);
        CHECK(processed.triangles[4] == correct);
        CHECK(processed.triangles[5] == correct);
    }

    SECTION("Triangle to Grid")
    {
        const std::vector<glm::vec2> triangle_left_hypo = { glm::vec2(10, 30), glm::vec2(30, 10), glm::vec2(50, 50) };
        const std::vector<glm::vec2> triangle_right_hypo = { glm::vec2(5, 5), glm::vec2(15, 10), glm::vec2(5, 15) };

        const std::vector<std::vector<glm::vec2>> triangle_points = { triangle_left_hypo, triangle_right_hypo };

        const std::vector<unsigned int> style_indices = { 1, 2 };

        const auto id = nucleus::tile::Id { .zoom_level = 10, .coords = { 548, 359 }, .scheme = nucleus::tile::Scheme::SlippyMap };

        Preprocessor p(id);
        auto processed = p.preprocess_triangles(triangle_points, style_indices);

        auto raster = p.visualize_grid(processed.cell_to_data);
        auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

        auto test_image = example_grid_data_triangles();
        CHECK(image == test_image);

        // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
        image.save(QString("vector_layer_grid_triangles.png"));
    }

    SECTION("Lines to Grid")
    {
        const std::vector<std::vector<glm::vec2>> line_points
            = { { glm::vec2(10.5, 40.5), glm::vec2(30.5, 20.5) }, { glm::vec2(10.5, 5.5), glm::vec2(30.5, 5.5) }, { glm::vec2(10.5, 50), glm::vec2(30.5, 50) } };
        const std::vector<unsigned int> style_indices = { 1, 2, 3 };

        const auto id = nucleus::tile::Id { .zoom_level = 10, .coords = { 548, 359 }, .scheme = nucleus::tile::Scheme::SlippyMap };

        Preprocessor p(id);
        auto processed = p.preprocess_lines(line_points, style_indices);

        auto raster = p.visualize_grid(processed.cell_to_data);

        auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

        // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
        image.save(QString("vector_layer_grid_lines.png"));
    }

    SECTION("vec2 to uint array conversion")
    {
        Triangle t;
        t.data.top_vertex = { 1.3434f, 2.5656f };
        t.data.middle_vertex = { 3.5656f, 4.111f };
        t.data.bottom_vertex = { 5.44432f, 6.0f };
        t.data.style_index = 7u;

        // only copy uint array and check if the union correctly converted the value to glm::vec2
        Triangle t2;
        std::copy(std::begin(t.packed), std::end(t.packed), std::begin(t2.packed));

        CHECK(t2.data.top_vertex == glm::vec2(1.3434f, 2.5656f));
        CHECK(t2.data.middle_vertex == glm::vec2(3.5656f, 4.111f));
        CHECK(t2.data.bottom_vertex == glm::vec2(5.44432f, 6.0f));
        CHECK(t2.data.style_index == 7u);
    }

    SECTION("vec2 to uint array conversion with vector insert")
    {
        auto triangles = std::vector<Triangle> { { { 1.3434f, 2.5656f }, { 3.5656f, 4.111f }, { 5.44432f, 6.0f }, 7u }, { { 10.3434f, 20.5656f }, { 30.5656f, 40.111f }, { 50.44432f, 60.0f }, 70u } };

        std::vector<uint32_t> triangles_ints;
        triangles_ints.reserve(triangles.size() * 7);
        for (size_t i = 0; i < triangles.size(); ++i) {
            triangles_ints.insert(triangles_ints.end(), triangles[i].packed, triangles[i].packed + 7);
        }

        Triangle t1;
        for (size_t i = 0; i < 7; ++i) {
            t1.packed[i] = triangles_ints[i];
        }

        CHECK(t1.data.top_vertex == glm::vec2(1.3434f, 2.5656f));
        CHECK(t1.data.middle_vertex == glm::vec2(3.5656f, 4.111f));
        CHECK(t1.data.bottom_vertex == glm::vec2(5.44432f, 6.0f));
        CHECK(t1.data.style_index == 7u);

        Triangle t2;
        for (size_t i = 7; i < 14; ++i) {
            t2.packed[i - 7] = triangles_ints[i];
        }

        // { 10.3434f, 20.5656f }, { 30.5656f, 40.111f }, { 30.5656f, 40.111f }, 70u }
        CHECK(t2.data.top_vertex == glm::vec2(10.3434f, 20.5656f));
        CHECK(t2.data.middle_vertex == glm::vec2(30.5656f, 40.111f));
        CHECK(t2.data.bottom_vertex == glm::vec2(50.44432f, 60.0f));
        CHECK(t2.data.style_index == 70u);
    }

    SECTION("16/24/24 Bit Data packer")
    {
        {
            auto test = glm::vec<3, uint32_t>(19u, 38u, 44u);
            auto packed = nucleus::utils::bit_coding::u16_u24_u24_to_u32_2(test.x, test.y, test.z);
            auto unpacked = nucleus::utils::bit_coding::u32_2_to_u16_u24_u24(packed);
            CHECK(unpacked == test);

            // std::cout << test << unpacked << std::endl;
        }
        {
            auto test = glm::vec<3, uint32_t>(46757u, 4825733u, 4711465u);
            auto packed = nucleus::utils::bit_coding::u16_u24_u24_to_u32_2(test.x, test.y, test.z);
            auto unpacked = nucleus::utils::bit_coding::u32_2_to_u16_u24_u24(packed);
            CHECK(unpacked == test);

            // std::cout << test << unpacked << std::endl;
        }
    }

    SECTION("24/8 Bit Data packer")
    {
        {
            auto test = glm::vec<2, uint32_t>(19u, 38u);
            auto packed = nucleus::utils::bit_coding::u24_u8_to_u32(test.x, test.y);
            auto unpacked = nucleus::utils::bit_coding::u32_to_u24_u8(packed);
            CHECK(unpacked == test);

            // std::cout << test << unpacked << std::endl;
        }
        {
            auto test = glm::vec<2, uint32_t>(10376621u, 250u);
            auto packed = nucleus::utils::bit_coding::u24_u8_to_u32(test.x, test.y);
            auto unpacked = nucleus::utils::bit_coding::u32_to_u24_u8(packed);
            CHECK(unpacked == test);

            // std::cout << test << unpacked << std::endl;
        }
    }
    Triangle t;
    t.packed[0] = 1106247680u;
    t.packed[1] = 1112014848u;
    t.packed[2] = 1092616192u;
    t.packed[3] = 1112014848u;
    t.packed[4] = -1u;
    t.packed[5] = -1u;
    t.packed[6] = -1u;

    std::cout << t << std::endl;
}
