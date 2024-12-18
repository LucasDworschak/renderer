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

#include <QSignalSpy>
#include <catch2/catch_test_macros.hpp>

#include <QFile>
#include <QImage>
#include <QString>
#include <glm/glm.hpp>

#include <nucleus/tile/TileLoadService.h>
#include <nucleus/tile/conversion.h>
#include <nucleus/tile/utils.h>
#include <nucleus/utils/bit_coding.h>
#include <radix/tile.h>

#include "nucleus/Raster.h"
#include "nucleus/vector_layer/Preprocessor.h"
#include "nucleus/vector_layer/Style.h"

#include "nucleus/utils/rasterizer.h"
#include "nucleus/vector_layer/constants.h"

using namespace nucleus::vector_layer;

// helpers for catch2
inline std::ostream& operator<<(std::ostream& os, const glm::uvec3& v) { return os << "{ " << v.x << ", " << v.y << ", " << v.z << " }"; }

inline std::ostream& operator<<(std::ostream& os, const glm::vec2& v) { return os << "{ " << v.x << ", " << v.y << " }"; }


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

nucleus::Raster<uint8_t> visualize_grid(const nucleus::vector_layer::details::VectorLayerGrid& processed_grid, int grid_width)
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

/*
 *  both raster either have 0 or some value
 *  the value does not have to be the same
 *  -> because raster1 is created by visualizegrid and filled with 255 and raster2 contains actual offset data that varies
 */
bool same_cells_are_filled(const nucleus::Raster<uint8_t>& raster1, std::shared_ptr<const nucleus::Raster<uint32_t>> raster2)
{
    REQUIRE(raster1.size() == raster2->size());

    auto b1 = raster1.buffer();
    auto b2 = raster2->buffer();

    for (size_t i = 0; i < b1.size(); ++i) {

        if ((b1[i] == 0) != (b2[i] == 0))
            return false;
    }

    return true;
}

TEST_CASE("nucleus/vector_preprocess")
{

    SECTION("Tile download basemap")
    {
        // if this fails it is very likely that something on the vector tile server changed
        // manually download the tile from the below link and check if the changes are valid and replace vectortile.mvt with this new file
        // https://osm.cg.tuwien.ac.at/vector_tiles/poi_v1/10/548/359

        const auto id = nucleus::tile::Id { .zoom_level = 10, .coords = { 548, 359 }, .scheme = nucleus::tile::Scheme::SlippyMap };
        nucleus::tile::TileLoadService service("https://mapsneu.wien.gv.at/basemapv/bmapv/3857/tile/", nucleus::tile::TileLoadService::UrlPattern::ZYX_yPointingSouth, ".pbf");

        {
            QSignalSpy spy(&service, &nucleus::tile::TileLoadService::load_finished);
            service.load(id);
            spy.wait(15000);

            REQUIRE(spy.count() == 1);
            QList<QVariant> arguments = spy.takeFirst();
            REQUIRE(arguments.size() == 1);
            nucleus::tile::Data tile = arguments.at(0).value<nucleus::tile::Data>();
            CHECK(tile.id == id);
            CHECK(tile.network_info.status == nucleus::tile::NetworkInfo::Status::Good);
            CHECK(nucleus::utils::time_since_epoch() - tile.network_info.timestamp < 10'000);

            REQUIRE(tile.data->size() > 0);
            CHECK(tile.data->size() > 2000);
        }
    }

    SECTION("Style download basemap")
    {
        Style s("https://mapsneu.wien.gv.at/basemapv/bmapv/3857/resources/styles/");
        QSignalSpy spy(&s, &Style::load_finished);
        s.load();
        spy.wait(10000);
    }

    SECTION("Style parsing color")
    {
        Style s("");
        CHECK(s.parse_color("rgba(134,179,1,0.5)") == 0x86B3017F);
        CHECK(s.parse_color("rgba(134,179,1,1)") == 0x86B301FF);
        CHECK(s.parse_color("rgb(134,179,1)") == 0x86B301FF);
        CHECK(s.parse_color("rgb(50,0,100)") == 0x320064FF);
        CHECK(s.parse_color("#FF55AA") == 0xFF55AAFF);
        CHECK(s.parse_color("#ABCDEF45") == 0xABCDEF45);
    }

    SECTION("Triangle to Grid")
    {
        // make sure that the proposed points are still roughly the same (even after tinkering with grid_size)
        // this would still mean that we have to redo the below data every time we chang the grid scale -> but it isn't too problematic for now
        constexpr float point_scale = 1.0f / 64.0f * nucleus::vector_layer::constants::grid_size;

        const std::vector<glm::vec2> triangle_left_hypo
            = { glm::vec2(10 * point_scale, 30 * point_scale), glm::vec2(30 * point_scale, 5 * point_scale), glm::vec2(50 * point_scale, 50 * point_scale) };
        const std::vector<glm::vec2> triangle_right_hypo = { glm::vec2(5 * point_scale, 5 * point_scale), glm::vec2(25 * point_scale, 10 * point_scale), glm::vec2(5 * point_scale, 15 * point_scale) };

        const std::vector<std::vector<glm::vec2>> triangle_points = { triangle_left_hypo, triangle_right_hypo };

        const std::vector<unsigned int> style_indices = { 1, 1 };

        auto processed = nucleus::vector_layer::details::preprocess_triangles(triangle_points, style_indices);

        auto raster = visualize_grid(processed.cell_to_temp, nucleus::vector_layer::constants::grid_size);
        auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

        auto test_image = example_grid_data_triangles();
        CHECK(image == test_image);

        auto tile = nucleus::vector_layer::details::create_gpu_tile(processed, processed);

        CHECK(same_cells_are_filled(raster, tile.grid_triangle));

        // we provide two triangles that overlap at one point -> we expect four entries [1], [0], ([1], [0])
        auto bridge_data = tile.grid_to_data->buffer();
        REQUIRE(bridge_data.size() == nucleus::vector_layer::constants::data_size * nucleus::vector_layer::constants::data_size);
        CHECK(bridge_data[0] == 1);
        CHECK(bridge_data[1] == 0);
        CHECK(bridge_data[2] == 1);
        CHECK(bridge_data[3] == 0);
        // the rest should be undefined -> -1u
        CHECK(bridge_data[4] == -1u);
        CHECK(bridge_data[5] == -1u);
        CHECK(bridge_data[nucleus::vector_layer::constants::data_size * 1.5] == -1u);


        auto data = tile.data_triangle->buffer();
        REQUIRE(data.size() == nucleus::vector_layer::constants::data_size * nucleus::vector_layer::constants::data_size);
        CHECK(data[0] == 1089470464);
        CHECK(data[1] == 1067450368);
        CHECK(data[2] == 1075838976);
        CHECK(data[3] == 1089470464);
        CHECK(data[4] == 1095237632);
        CHECK(data[5] == 1095237632);
        CHECK(data[6] == 1);
        CHECK(data[7] == 1067450368);
        CHECK(data[8] == 1067450368);
        CHECK(data[9] == 1086849024);
        CHECK(data[10] == 1075838976);
        CHECK(data[11] == 1067450368);
        CHECK(data[12] == 1081081856);
        CHECK(data[13] == 1);
        // the rest should be undefined -> -1u
        CHECK(data[14] == -1u);
        CHECK(data[15] == -1u);
        CHECK(data[nucleus::vector_layer::constants::data_size * 1.7] == -1u);

        // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
        //     std::cout << bridge_data[i] << std::endl;
        // }

        // std::cout << std::endl;
        // for (int i = 0; i < 14; ++i) { // DEBUG expected triangle data
        //     std::cout << data[i] << std::endl;
        // }

        // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
        // image.save(QString("vector_layer_grid_triangles.png"));
    }

    SECTION("Triangle to Grid (outside vertices)")
    {
        // make sure that the proposed points are still roughly the same (even after tinkering with grid_size)
        // this would still mean that we have to redo the below data every time we chang the grid scale -> but it isn't too problematic for now
        constexpr float point_scale = 1.0f / 64.0f * nucleus::vector_layer::constants::grid_size;

        const std::vector<glm::vec2> triangle_left_hypo
            = { glm::vec2(-10 * point_scale, 30 * point_scale), glm::vec2(30 * point_scale, 5 * point_scale), glm::vec2(50 * point_scale, 80 * point_scale) };
        const std::vector<glm::vec2> triangle_right_hypo
            = { glm::vec2(5 * point_scale, -5 * point_scale), glm::vec2(90 * point_scale, 10 * point_scale), glm::vec2(5 * point_scale, 15 * point_scale) };

        const std::vector<std::vector<glm::vec2>> triangle_points = { triangle_left_hypo, triangle_right_hypo };

        const std::vector<unsigned int> style_indices = { 1, 1 };

        auto processed = nucleus::vector_layer::details::preprocess_triangles(triangle_points, style_indices);

        auto raster = visualize_grid(processed.cell_to_temp, nucleus::vector_layer::constants::grid_size);
        auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

        auto tile = nucleus::vector_layer::details::create_gpu_tile(processed, processed);

        CHECK(same_cells_are_filled(raster, tile.grid_triangle));

        // we provide two triangles that overlap at one point -> we expect four entries [1], ([1], [0]), [1]
        auto bridge_data = tile.grid_to_data->buffer();
        REQUIRE(bridge_data.size() == nucleus::vector_layer::constants::data_size * nucleus::vector_layer::constants::data_size);
        CHECK(bridge_data[0] == 1);
        CHECK(bridge_data[1] == 1);
        CHECK(bridge_data[2] == 0);
        CHECK(bridge_data[3] == 0);
        // the rest should be undefined -> -1u
        CHECK(bridge_data[4] == -1u);
        CHECK(bridge_data[5] == -1u);
        CHECK(bridge_data[nucleus::vector_layer::constants::data_size * 1.5] == -1u);

        auto data = tile.data_triangle->buffer();
        REQUIRE(data.size() == nucleus::vector_layer::constants::data_size * nucleus::vector_layer::constants::data_size);
        CHECK(data[0] == 1089470464);
        CHECK(data[1] == 1067450368);
        CHECK(data[2] == 3223322624);
        CHECK(data[3] == 1089470464);
        CHECK(data[4] == 1095237632);
        CHECK(data[5] == 1101004800);
        CHECK(data[6] == 1);
        CHECK(data[7] == 1067450368);
        CHECK(data[8] == 3214934016);
        CHECK(data[9] == 1102315520);
        CHECK(data[10] == 1075838976);
        CHECK(data[11] == 1067450368);
        CHECK(data[12] == 1081081856);
        CHECK(data[13] == 1);
        // the rest should be undefined -> -1u
        CHECK(data[14] == -1u);
        CHECK(data[15] == -1u);
        CHECK(data[nucleus::vector_layer::constants::data_size * 1.7] == -1u);

        // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
        //     std::cout << bridge_data[i] << std::endl;
        // }

        // std::cout << std::endl;
        // for (int i = 0; i < 14; ++i) { // DEBUG expected triangle data
        //     std::cout << data[i] << std::endl;
        // }

        // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
        // image.save(QString("vector_layer_grid_triangles_outside.png"));
    }

    SECTION("Basemap decoding - vectortile Neusiedlersee")
    {
        // TODO test vectortile_neusiedlersee.pbf and see what is wrong with the rasterization of the lake
        // const auto id = nucleus::tile::Id { .zoom_level = 14, .coords = { 4477 * 2, 2850 * 2 + 1 }, .scheme = nucleus::tile::Scheme::SlippyMap };
        const auto id = nucleus::tile::Id { .zoom_level = 13, .coords = { 4477, 2850 }, .scheme = nucleus::tile::Scheme::SlippyMap };

        QByteArray byte_data;

        {
            nucleus::tile::TileLoadService service("https://mapsneu.wien.gv.at/basemapv/bmapv/3857/tile/", nucleus::tile::TileLoadService::UrlPattern::ZYX_yPointingSouth, ".pbf");
            QSignalSpy spy(&service, &nucleus::tile::TileLoadService::load_finished);
            service.load(id);
            spy.wait(15000);

            REQUIRE(spy.count() == 1);
            QList<QVariant> arguments = spy.takeFirst();
            REQUIRE(arguments.size() == 1);
            nucleus::tile::Data tile = arguments.at(0).value<nucleus::tile::Data>();
            byte_data = *tile.data;
        }

        // auto file = QFile(QString("%1%2").arg(ALP_TEST_DATA_DIR, "vectortile_neusiedlersee.pbf"));
        // file.open(QFile::ReadOnly);
        // const auto byte_data = file.readAll();

        Style style("");

        // auto processed = nucleus::vector_layer::preprocess(id, byte_data, style);
        auto polygons = nucleus::vector_layer::details::parse_tile(id, byte_data, style);

        // constexpr float point_scale = 1.0f / 64.0f;

        // const std::vector<glm::vec2> triangle_left_hypo
        //     = { glm::vec2(10 * point_scale, 30 * point_scale), glm::vec2(30 * point_scale, 5 * point_scale), glm::vec2(50 * point_scale, 50 * point_scale) };
        // const std::vector<glm::vec2> triangle_right_hypo = { glm::vec2(5 * point_scale, 5 * point_scale), glm::vec2(25 * point_scale, 10 * point_scale), glm::vec2(5 * point_scale, 15 * point_scale)
        // };

        // const std::vector<std::vector<glm::vec2>> polygons = { triangle_left_hypo, triangle_right_hypo };

        size_t data_offset = 1;

        for (size_t i = 0; i < polygons.size(); ++i) {
            std::vector<glm::vec2> triangle_points = nucleus::utils::rasterizer::triangulize(polygons[i], true);

            std::cout << "o Water" << i << std::endl;
            for (size_t j = 0; j < triangle_points.size() / 3; ++j) {

                std::cout << "v " << triangle_points[j * 3 + 0].x << " 0.0 " << triangle_points[j * 3 + 0].y << std::endl;
                std::cout << "v " << triangle_points[j * 3 + 1].x << " 0.0 " << triangle_points[j * 3 + 1].y << std::endl;
                std::cout << "v " << triangle_points[j * 3 + 2].x << " 0.0 " << triangle_points[j * 3 + 2].y << std::endl;
            }

            for (size_t j = 0; j < triangle_points.size() / 3; ++j) {
                std::cout << "f " << (data_offset + j * 3 + 0) << " " << (data_offset + j * 3 + 1) << " " << (data_offset + j * 3 + 2) << std::endl;
            }
            data_offset += triangle_points.size();
        }

        // TODO here -> visualize the grid first,
        // probably best to execute each step individually and rasterize the resulting polygons with own raster with higher level of detail than grid size

        // auto raster = processed.grid_triangle;
        // auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

        // auto test_image = example_grid_data_triangles();
        // CHECK(image == test_image);

        // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
        //     std::cout << bridge_data[i] << std::endl;
        // }

        // std::cout << std::endl;
        // for (int i = 0; i < 14; ++i) { // DEBUG expected triangle data
        //     std::cout << data[i] << std::endl;
        // }

        // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
        // image.save(QString("vector_layer_grid_triangles.png"));
    }

    SECTION("Lines to Grid")
    {
        const std::vector<std::vector<glm::vec2>> line_points
            = { { glm::vec2(10.5, 40.5), glm::vec2(30.5, 20.5) }, { glm::vec2(10.5, 5.5), glm::vec2(30.5, 5.5) }, { glm::vec2(10.5, 50), glm::vec2(30.5, 50) } };
        const std::vector<unsigned int> style_indices = { 1, 2, 3 };

        // const auto id = nucleus::tile::Id { .zoom_level = 10, .coords = { 548, 359 }, .scheme = nucleus::tile::Scheme::SlippyMap };

        auto processed = nucleus::vector_layer::details::preprocess_lines(line_points, style_indices);

        auto raster = visualize_grid(processed.cell_to_temp, nucleus::vector_layer::constants::grid_size);

        auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

        // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
        // image.save(QString("vector_layer_grid_lines.png"));
    }

    SECTION("float to uint array conversion")
    {
        auto f0 = 1.3434f;
        auto f1 = 3.5656f;
        auto f2 = 5.44432f;

        uint32_t u0 = *reinterpret_cast<uint32_t*>(&f0);
        uint32_t u1 = *reinterpret_cast<uint32_t*>(&f1);
        uint32_t u2 = *reinterpret_cast<uint32_t*>(&f2);

        // convert back to values we can test against
        float t0 = *reinterpret_cast<float*>(&u0);
        float t1 = *reinterpret_cast<float*>(&u1);
        float t2 = *reinterpret_cast<float*>(&u2);

        CHECK(t0 == 1.3434f);
        CHECK(t1 == 3.5656f);
        CHECK(t2 == 5.44432f);
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
