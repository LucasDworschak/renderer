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

        std::cout << service.build_tile_url(id).toStdString();

        {
            QSignalSpy spy(&service, &nucleus::tile::TileLoadService::load_finished);
            service.load(id);
            spy.wait(10000);

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
        const std::vector<glm::vec2> triangle_left_hypo = { glm::vec2(10, 30), glm::vec2(30, 10), glm::vec2(50, 50) };
        const std::vector<glm::vec2> triangle_right_hypo = { glm::vec2(5, 5), glm::vec2(15, 10), glm::vec2(5, 15) };

        const std::vector<std::vector<glm::vec2>> triangle_points = { triangle_left_hypo, triangle_right_hypo };

        const std::vector<unsigned int> style_indices = { 1, 2 };

        auto processed = nucleus::vector_layer::details::preprocess_triangles(triangle_points, style_indices);

        auto raster = visualize_grid(processed.cell_to_temp, 64);
        auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

        auto test_image = example_grid_data_triangles();
        CHECK(image == test_image);

        auto tile = nucleus::vector_layer::details::create_gpu_tile(processed, processed);

        CHECK(same_cells_are_filled(raster, tile.grid_triangle));

        // we provide two triangles that overlap at one point -> we expect four entries ([1], [0], ([1], [0])
        REQUIRE(tile.grid_to_data->size() == 4);
        CHECK(tile.grid_to_data->at(0) == 1);
        CHECK(tile.grid_to_data->at(1) == 0);
        CHECK(tile.grid_to_data->at(2) == 1);
        CHECK(tile.grid_to_data->at(3) == 0);

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

        auto raster = visualize_grid(processed.cell_to_temp, 64);

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
