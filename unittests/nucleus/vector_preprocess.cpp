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
#include <catch2/benchmark/catch_benchmark.hpp>
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

#include "nucleus/utils/rasterizer.h"
#include "nucleus/vector_layer/constants.h"

using namespace nucleus::vector_layer;

// helpers for catch2
inline std::ostream& operator<<(std::ostream& os, const glm::uvec3& v) { return os << "{ " << v.x << ", " << v.y << ", " << v.z << " }"; }

inline std::ostream& operator<<(std::ostream& os, const glm::vec2& v) { return os << "{ " << v.x << ", " << v.y << " }"; }
inline std::ostream& operator<<(std::ostream& os, const glm::ivec2& v) { return os << "{ " << v.x << ", " << v.y << " }"; }

QImage example_grid_data_triangles()
{
    auto file = QFile(QString("%1%2").arg(ALP_TEST_DATA_DIR, "vector_layer/grid_triangles.png"));
    file.open(QFile::ReadOnly);
    const auto bytes = file.readAll();
    auto image = QImage::fromData(bytes);
    REQUIRE(!image.isNull());
    return image;
}

QImage example_grid_data_lines()
{
    auto file = QFile(QString("%1%2").arg(ALP_TEST_DATA_DIR, "vector_layer/grid_lines.png"));
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

    SECTION("Triangle to Grid")
    {
        constexpr auto extent = 64u;
        const std::vector<glm::vec2> triangle_left_hypo = { glm::vec2(10, 30), glm::vec2(30, 5), glm::vec2(50, 50) };
        const std::vector<glm::vec2> triangle_right_hypo = { glm::vec2(5, 5), glm::vec2(25, 10), glm::vec2(5, 15) };

        std::vector<nucleus::vector_layer::details::GeometryData> tile_data {
            { triangle_left_hypo, extent, 1, 1, true, nucleus::utils::rasterizer::generate_neighbour_edges(triangle_left_hypo.size()), 0 },
            { triangle_right_hypo, extent, 1, 1, true, nucleus::utils::rasterizer::generate_neighbour_edges(triangle_right_hypo.size()), 0 }
        };

        auto processed = nucleus::vector_layer::details::preprocess_geometry(tile_data);

        auto raster = visualize_grid(processed.acceleration_grid, nucleus::vector_layer::constants::grid_size);
        auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

        auto test_image = example_grid_data_triangles();
        CHECK(image == test_image);

        auto tile = nucleus::vector_layer::details::create_gpu_tile(processed);

        CHECK(same_cells_are_filled(raster, tile.acceleration_grid));

        // we provide two triangles that overlap at one point -> we expect four entries [1], [0], ([1], [0])
        auto bridge_data = tile.index_buffer->buffer();
        REQUIRE(bridge_data.size() == nucleus::vector_layer::constants::data_size[tile.buffer_info] * nucleus::vector_layer::constants::data_size[tile.buffer_info]);
        CHECK(bridge_data[0] == 3); // 1 and 3 values since lowest bit is is_polygon flag
        CHECK(bridge_data[1] == 1);
        CHECK(bridge_data[2] == 1);
        CHECK(bridge_data[3] == 3);
        // the rest should be undefined -> -1u
        CHECK(bridge_data[4] == -1u);
        CHECK(bridge_data[5] == -1u);
        CHECK(bridge_data[nucleus::vector_layer::constants::data_size[tile.buffer_info] * 1.5] == -1u);

        auto data = tile.vertex_buffer->buffer();
        REQUIRE(data.size() == nucleus::vector_layer::constants::data_size[tile.buffer_info] * nucleus::vector_layer::constants::data_size[tile.buffer_info]);
        CHECK(data[0].x == 1081671760);
        CHECK(data[0].y == 1076429280);
        CHECK(data[0].z == 1086915361);
        CHECK(data[1].x == 1075118160);
        CHECK(data[1].y == 1080361120);
        CHECK(data[1].z == 1075118321);
        // the rest should be undefined -> -1u
        CHECK(data[2].x == -1u);
        CHECK(data[3].y == -1u);
        CHECK(data[nucleus::vector_layer::constants::data_size[tile.buffer_info] * 1.7].x == -1u);

        // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
        //     std::cout << bridge_data[i] << std::endl;
        // }

        // std::cout << std::endl;
        // for (int i = 0; i < 14; ++i) { // DEBUG expected triangle data
        //     std::cout << data[i] << std::endl;
        // }

        // // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
        // image.save(QString("vector_layer_grid_triangles.png"));
    }

    SECTION("Triangle to Grid (outside vertices)")
    {
        // make sure that the proposed points are still roughly the same (even after tinkering with grid_size)
        // this would still mean that we have to redo the below data every time we chang the grid scale -> but it isn't too problematic for now
        constexpr auto extent = 64u;

        const std::vector<glm::vec2> triangle_left_hypo = { glm::vec2(-10, 30), glm::vec2(30, 5), glm::vec2(50, 80) };
        const std::vector<glm::vec2> triangle_right_hypo = { glm::vec2(5, -5), glm::vec2(90, 10), glm::vec2(5, 15) };

        std::vector<nucleus::vector_layer::details::GeometryData> tile_data {
            { triangle_left_hypo, extent, 1, 1, true, nucleus::utils::rasterizer::generate_neighbour_edges(triangle_left_hypo.size()), 0 },
            { triangle_right_hypo, extent, 1, 1, true, nucleus::utils::rasterizer::generate_neighbour_edges(triangle_right_hypo.size()), 0 }
        };

        auto processed = nucleus::vector_layer::details::preprocess_geometry(tile_data);

        auto raster = visualize_grid(processed.acceleration_grid, nucleus::vector_layer::constants::grid_size);
        auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

        auto tile = nucleus::vector_layer::details::create_gpu_tile(processed);

        CHECK(same_cells_are_filled(raster, tile.acceleration_grid));

        // we provide two triangles that overlap at one point -> we expect four entries [1], ([1], [0]), [1]
        auto bridge_data = tile.index_buffer->buffer();
        REQUIRE(bridge_data.size() == nucleus::vector_layer::constants::data_size[tile.buffer_info] * nucleus::vector_layer::constants::data_size[tile.buffer_info]);
        CHECK(bridge_data[0] == 3); // 1 and 3 values since lowest bit is is_polygon flag
        CHECK(bridge_data[1] == 1);
        CHECK(bridge_data[2] == 3);
        CHECK(bridge_data[3] == 1);
        // the rest should be undefined -> -1u
        CHECK(bridge_data[4] == -1u);
        CHECK(bridge_data[5] == -1u);
        CHECK(bridge_data[nucleus::vector_layer::constants::data_size[tile.buffer_info] * 1.5] == -1u);

        auto data = tile.vertex_buffer->buffer();
        REQUIRE(data.size() == nucleus::vector_layer::constants::data_size[tile.buffer_info] * nucleus::vector_layer::constants::data_size[tile.buffer_info]);
        CHECK(data[0].x == 1081671760);
        CHECK(data[0].y == 1071186400);
        CHECK(data[0].z == 1086915841);
        CHECK(data[1].x == 1075118000);
        CHECK(data[1].y == 1097400480);
        CHECK(data[1].z == 1075118321);
        // the rest should be undefined -> -1u
        CHECK(data[2].x == -1u);
        CHECK(data[3].y == -1u);

        CHECK(data[nucleus::vector_layer::constants::data_size[tile.buffer_info] * 1.7].x == -1u);

        // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
        //     std::cout << bridge_data[i] << std::endl;
        // }

        // std::cout << std::endl;
        // for (int i = 0; i < 14; ++i) { // DEBUG expected triangle data
        //     std::cout << data[i] << std::endl;
        // }

        // // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
        // image.save(QString("vector_layer_grid_triangles_outside.png"));
    }

    SECTION("Line to Grid")
    {
        constexpr auto extent = 64;
        const std::vector<glm::vec2> line0 = { glm::vec2(10, 30), glm::vec2(30, 50), glm::vec2(50, 30) };

        std::vector<nucleus::vector_layer::details::GeometryData> tile_data { { line0, extent, 1, 1, false, {}, 0 } };

        auto processed = nucleus::vector_layer::details::preprocess_geometry(tile_data);

        auto raster = visualize_grid(processed.acceleration_grid, nucleus::vector_layer::constants::grid_size);
        auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

        // auto test_image = example_grid_data_triangles();
        // CHECK(image == test_image);

        auto tile = nucleus::vector_layer::details::create_gpu_tile(processed);

        CHECK(same_cells_are_filled(raster, tile.acceleration_grid));

        // we provide two line segments that overlap at one point -> we expect four entries [0], [2], ([0], [2])
        auto bridge_data = tile.index_buffer->buffer();
        REQUIRE(bridge_data.size() == nucleus::vector_layer::constants::data_size[tile.buffer_info] * nucleus::vector_layer::constants::data_size[tile.buffer_info]);
        CHECK(bridge_data[0] == 0); // 0 and 2 values since lowest bit is is_polygon flag
        CHECK(bridge_data[1] == 2);
        CHECK(bridge_data[2] == 0);
        CHECK(bridge_data[3] == 2);
        // the rest should be undefined -> -1u
        CHECK(bridge_data[4] == -1u);
        CHECK(bridge_data[5] == -1u);
        CHECK(bridge_data[nucleus::vector_layer::constants::data_size[tile.buffer_info] * 1.5] == -1u);

        auto data = tile.vertex_buffer->buffer();
        REQUIRE(data.size() == nucleus::vector_layer::constants::data_size[tile.buffer_info] * nucleus::vector_layer::constants::data_size[tile.buffer_info]);
        CHECK(data[0].x == 1076429280);
        CHECK(data[0].y == 1081672480);
        CHECK(data[0].z == 1073807361);
        CHECK(data[1].x == 1081672480);
        CHECK(data[1].y == 1086915040);
        CHECK(data[1].z == 1073807361);
        // the rest should be undefined -> -1u
        CHECK(data[2].x == -1u);
        CHECK(data[3].y == -1u);

        CHECK(data[nucleus::vector_layer::constants::data_size[tile.buffer_info] * 1.7].x == -1u);

        // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
        //     std::cout << bridge_data[i] << std::endl;
        // }

        // std::cout << std::endl;
        // for (int i = 0; i < 14; ++i) { // DEBUG expected triangle data
        //     std::cout << data[i] << std::endl;
        // }

        // // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
        image.save(QString("vector_layer_grid_lines.png"));
    }

    // SECTION("Tile exploration") // section mostly used for tile debugging -> not a real test
    // {
    //     // const auto id = nucleus::tile::Id { .zoom_level = 14, .coords = { 4477 * 2, 2850 * 2 + 1 }, .scheme = nucleus::tile::Scheme::SlippyMap };
    //     // const auto id = nucleus::tile::Id { .zoom_level = 13, .coords = { 4477, 2850 }, .scheme = nucleus::tile::Scheme::SlippyMap };
    //     // const auto id = nucleus::tile::Id { .zoom_level = 15, .coords = { (4477 * 2 + 1) * 2, (2850 * 2 + 1) * 2 }, .scheme = nucleus::tile::Scheme::SlippyMap };
    //     const auto id = nucleus::tile::Id { .zoom_level = 19, .coords = { 285987, 181795 }, .scheme = nucleus::tile::Scheme::SlippyMap }; // bodensee

    //     QByteArray byte_data;

    //     {
    //         nucleus::tile::TileLoadService service("https://osm.cg.tuwien.ac.at/vector_tiles/vector_layer_v1/", nucleus::tile::TileLoadService::UrlPattern::ZXY_yPointingSouth, "");
    //         QSignalSpy spy(&service, &nucleus::tile::TileLoadService::load_finished);
    //         service.load(id);
    //         spy.wait(15000);

    //         REQUIRE(spy.count() == 1);
    //         QList<QVariant> arguments = spy.takeFirst();
    //         REQUIRE(arguments.size() == 1);
    //         nucleus::tile::Data tile = arguments.at(0).value<nucleus::tile::Data>();
    //         byte_data = *tile.data;
    //     }

    //     Style s(":/vectorlayerstyles/openstreetmap.json");
    //     s.load();

    //     auto tile_data = nucleus::vector_layer::details::parse_tile(id, byte_data, s);

    //     size_t data_offset = 1;
    //     auto acceleration_grid = std::vector<std::set<uint32_t>>(nucleus::vector_layer::constants::grid_size * nucleus::vector_layer::constants::grid_size, std::set<uint32_t>());

    //     for (size_t i = 0; i < tile_data.size(); ++i) {
    //         if (tile_data[i].is_polygon) {

    //             if (tile_data[i].style != 88u)
    //                 continue;

    //             std::vector<glm::vec2> triangle_points = nucleus::utils::rasterizer::triangulize(tile_data[i].vertices, tile_data[i].edges, true);

    //             std::cout << std::endl << "o pedestrian_area" << i << std::endl;
    //             for (size_t j = 0; j < triangle_points.size() / 3; ++j) {
    //                 std::cout << "v " << triangle_points[j * 3 + 0].x << " 0 " << triangle_points[j * 3 + 0].y << std::endl;
    //                 std::cout << "v " << triangle_points[j * 3 + 1].x << " 0 " << triangle_points[j * 3 + 1].y << std::endl;
    //                 std::cout << "v " << triangle_points[j * 3 + 2].x << " 0 " << triangle_points[j * 3 + 2].y << std::endl;

    //                 auto f_ind = ((j + data_offset) * 3);
    //                 std::cout << "f " << f_ind << " " << (f_ind + 1) << " " << (f_ind + 2) << std::endl;
    //             }

    //             data_offset += triangle_points.size() / 3;

    //         } else {
    //             continue;
    //             // qDebug() << tile_data[i].style;
    //             if (tile_data[i].style != 92u)
    //                 continue;
    //             const auto cell_writer = [&acceleration_grid, data_offset](glm::vec2 pos, int data_index) {
    //                 // if in grid_size bounds and not already present -> than add index to vector
    //                 if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(nucleus::vector_layer::constants::grid_size), pos))) {
    //                     // last bit of index indicates that this is a line
    //                     // !! IMPORTANT !! we have to use the lowest bit since set orders the input depending on key -> lines and polygons have to stay intermixed
    //                     const auto index = ((data_index + data_offset) << 1); // | 0u;
    //                     acceleration_grid[int(pos.x) + nucleus::vector_layer::constants::grid_size * int(pos.y)].insert(index);
    //                 }
    //             };

    //             std::cout << std::endl << "o road" << i << std::endl;
    //             for (size_t j = 0; j < tile_data[i].vertices.size(); ++j) {
    //                 std::cout << "v " << tile_data[i].vertices[j].x << " 0 " << tile_data[i].vertices[j].y << std::endl;
    //             }
    //             std::cout << "l ";
    //             for (size_t j = 0; j < tile_data[i].vertices.size(); ++j) {
    //                 std::cout << (j + data_offset) << " ";
    //             }

    //             std::cout << std::endl << std::endl;

    //             const auto scale = float(nucleus::vector_layer::constants::grid_size) / float(tile_data[i].extent);
    //             nucleus::utils::rasterizer::rasterize_line(cell_writer, tile_data[i].vertices, tile_data[i].line_width * scale, scale);

    //             data_offset += tile_data[i].vertices.size();
    //         }
    //     }

    //     auto raster = visualize_grid(acceleration_grid, nucleus::vector_layer::constants::grid_size);
    //     auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);
    //     image.save(QString("vector_layer_debuggggg.png"));

    //     // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
    //     //     std::cout << bridge_data[i] << std::endl;
    //     // }

    //     // std::cout << std::endl;
    //     // for (int i = 0; i < 14; ++i) { // DEBUG expected triangle data
    //     //     std::cout << data[i] << std::endl;
    //     // }
    // }

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

    SECTION("Triangle Data packer")
    {
        {
            // auto a = glm::ivec2(-4095, 0b000000000000);
            // auto b = glm::ivec2(4095, 0b000000000000);
            // auto c = glm::ivec2(-1, 0b000000000000);
            // uint32_t style = 0b000000000000;

            auto a = glm::ivec2(0b110101010100, 0b111101110100);
            auto b = glm::ivec2(0b001111000101, 0b001111101101);
            auto c = glm::ivec2(0b010110011011, 0b100111001110);
            uint32_t style = 0b110001110001;

            auto packed = nucleus::vector_layer::details::pack_triangle_data(a, b, c, style);
            auto unpacked = nucleus::vector_layer::details::unpack_triangle_data(packed);

            CHECK(a == glm::ivec2(std::get<0>(unpacked)));
            CHECK(b == glm::ivec2(std::get<1>(unpacked)));
            CHECK(c == glm::ivec2(std::get<2>(unpacked)));
            CHECK(style == std::get<3>(unpacked));

            // std::cout << packed << std::endl;
            // std::cout << a.x << " " << a.y << std::endl;
            // std::cout << b.x << " " << b.y << std::endl;
            // std::cout << c.x << " " << c.y << std::endl;
            // std::cout << std::get<3>(unpacked) << std::endl;
            // std::cout << int32_t(std::get<0>(unpacked).x) << " " << int32_t(std::get<0>(unpacked).y) << std::endl;
            // std::cout << int32_t(std::get<1>(unpacked).x) << " " << int32_t(std::get<1>(unpacked).y) << std::endl;
            // std::cout << int32_t(std::get<2>(unpacked).x) << " " << int32_t(std::get<2>(unpacked).y) << std::endl;
            // std::cout << std::get<3>(unpacked) << std::endl;
        }

        {
            // test with negative numbers
            auto a = glm::ivec2(-nucleus::vector_layer::constants::tile_extent / 3.0, 0);
            auto c = glm::ivec2(0, -nucleus::vector_layer::constants::tile_extent / 4.0);
            auto b = glm::ivec2(nucleus::vector_layer::constants::tile_extent + -nucleus::vector_layer::constants::tile_extent / 5.0, -1);
            uint32_t style = 0;

            auto packed = nucleus::vector_layer::details::pack_triangle_data(a, b, c, style);
            auto unpacked = nucleus::vector_layer::details::unpack_triangle_data(packed);

            CHECK(a == glm::ivec2(std::get<0>(unpacked)));
            CHECK(b == glm::ivec2(std::get<1>(unpacked)));
            CHECK(c == glm::ivec2(std::get<2>(unpacked)));
            CHECK(style == std::get<3>(unpacked));

            // std::cout << packed << std::endl;
            // std::cout << a.x << " " << a.y << std::endl;
            // std::cout << b.x << " " << b.y << std::endl;
            // std::cout << c.x << " " << c.y << std::endl;
            // std::cout << std::get<3>(unpacked) << std::endl;
            // std::cout << int32_t(std::get<0>(unpacked).x) << " " << int32_t(std::get<0>(unpacked).y) << std::endl;
            // std::cout << int32_t(std::get<1>(unpacked).x) << " " << int32_t(std::get<1>(unpacked).y) << std::endl;
            // std::cout << int32_t(std::get<2>(unpacked).x) << " " << int32_t(std::get<2>(unpacked).y) << std::endl;
            // std::cout << std::get<3>(unpacked) << std::endl;
        }

        {
            // test with negative numbers
            auto a = glm::ivec2(0b010101010100 * -1, 0b011101110100);
            auto b = glm::ivec2(0b001111000101, 0b001111101101 * -1);
            auto c = glm::ivec2(0b010110011011 * -1, 0b010111001110 * -1);
            uint32_t style = 0b00100111000101;

            auto packed = nucleus::vector_layer::details::pack_triangle_data(a, b, c, style);
            auto unpacked = nucleus::vector_layer::details::unpack_triangle_data(packed);

            CHECK(a == glm::ivec2(std::get<0>(unpacked)));
            CHECK(b == glm::ivec2(std::get<1>(unpacked)));
            CHECK(c == glm::ivec2(std::get<2>(unpacked)));
            CHECK(style == std::get<3>(unpacked));

            // std::cout << packed << std::endl;
            // std::cout << a.x << " " << a.y << std::endl;
            // std::cout << b.x << " " << b.y << std::endl;
            // std::cout << c.x << " " << c.y << std::endl;
            // std::cout << std::get<3>(unpacked) << std::endl;
            // std::cout << int32_t(std::get<0>(unpacked).x) << " " << int32_t(std::get<0>(unpacked).y) << std::endl;
            // std::cout << int32_t(std::get<1>(unpacked).x) << " " << int32_t(std::get<1>(unpacked).y) << std::endl;
            // std::cout << int32_t(std::get<2>(unpacked).x) << " " << int32_t(std::get<2>(unpacked).y) << std::endl;
            // std::cout << std::get<3>(unpacked) << std::endl;
        }
    }
}

TEST_CASE("nucleus/vector_preprocess benchmarks")
{
    // BENCHMARK("triangulize polygons")
    // {
    //     const std::vector<glm::vec2> polygon_points = { glm::vec2(10.5, 10.5), glm::vec2(30.5, 10.5), glm::vec2(50.5, 50.5), glm::vec2(10.5, 30.5) };
    //     const auto edges = nucleus::utils::rasterizer::generate_neighbour_edges(polygon_points);
    //     nucleus::utils::rasterizer::triangulize(polygon_points, edges);
    // };
}
