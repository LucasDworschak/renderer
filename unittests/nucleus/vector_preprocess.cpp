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


#include <nucleus/tile/conversion.h>
#include <nucleus/utils/bit_coding.h>
#include <radix/tile.h>

#include "nucleus/Raster.h"
#include "nucleus/vector_layer/Preprocessor.h"
#include "nucleus/vector_tile/types.h"

using namespace nucleus::vector_tile;
using namespace nucleus::vector_layer;

// helpers for catch2
inline std::ostream& operator<<(std::ostream& os, const glm::uvec3& v) { return os << "{ " << v.x << ", " << v.y << ", " << v.z << " }"; }

inline std::ostream& operator<<(std::ostream& os, const glm::vec2& v) { return os << "{ " << v.x << ", " << v.y << " }"; }

inline std::ostream& operator<<(std::ostream& os, const Triangle& t)
{
    return os << "{ top: " << t.data.top_vertex << ", middle: " << t.data.middle_vertex << ", bottom: " << t.data.bottom_vertex << ", style: " << std::to_string(t.data.style_index) << " }";
}

inline std::ostream& operator<<(std::ostream& os, const Line& t)
{
    return os << "{ start: " << t.data.line_start_vertex << ", end: " << t.data.line_end_vertex << ", style: " << std::to_string(t.data.style_index) << " }";
}

inline bool operator==(const Triangle& t1, const Triangle& t2)
{
    return t1.data.top_vertex == t2.data.top_vertex && t1.data.middle_vertex == t2.data.middle_vertex && t1.data.bottom_vertex == t2.data.bottom_vertex && t1.data.style_index == t2.data.style_index;
}

inline bool operator==(const Line& l1, const Line& l2)
{
    return l1.data.line_start_vertex == l2.data.line_start_vertex && l1.data.line_end_vertex == l2.data.line_end_vertex && l1.data.style_index == l2.data.style_index;
}

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

nucleus::Raster<uint8_t> visualize_grid(const VectorLayerGrid& processed_grid, int grid_width)
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

    return nucleus::Raster<uint8_t>(grid_width, std::move(output_grid));
}

TEST_CASE("nucleus/vector_preprocess")
{

    SECTION("Triangle to Grid")
    {
        const std::vector<glm::vec2> triangle_left_hypo = { glm::vec2(10, 30), glm::vec2(30, 10), glm::vec2(50, 50) };
        const std::vector<glm::vec2> triangle_right_hypo = { glm::vec2(5, 5), glm::vec2(15, 10), glm::vec2(5, 15) };

        const std::vector<std::vector<glm::vec2>> triangle_points = { triangle_left_hypo, triangle_right_hypo };

        const std::vector<unsigned int> style_indices = { 1, 2 };

        const auto id = nucleus::tile::Id { .zoom_level = 10, .coords = { 548, 359 }, .scheme = nucleus::tile::Scheme::SlippyMap };

        Preprocessor p(id);
        auto processed = p.preprocess_triangles(triangle_points, style_indices);

        auto raster = visualize_grid(processed.cell_to_data, 64);
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

        auto raster = visualize_grid(processed.cell_to_data, 64);

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
}
