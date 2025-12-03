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


#include <CDT.h>

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

    template <>
    struct nth<0, nucleus::vector_layer::ClipperPoint> {
        inline static auto get(const nucleus::vector_layer::ClipperPoint& t) { return t.x; };
    };
    template <>
    struct nth<1, nucleus::vector_layer::ClipperPoint> {
        inline static auto get(const nucleus::vector_layer::ClipperPoint& t) { return t.y; };
    };

} // namespace util
} // namespace mapbox

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

void visualize_acceleration_grid(std::shared_ptr<const nucleus::Raster<uint32_t>> raster)
{
    std::vector<uint8_t> output_grid;

    const auto& grid = raster->buffer();

    for (size_t i = 0; i < grid.size(); ++i) {

        const auto data = nucleus::utils::bit_coding::u32_to_u24_u8(grid[i]);

        // grid is a singular uint32_t value that encodes the start index of the triangle list and the amount of triangles
        if (data.y == 0) {
            output_grid.push_back(0);
        } else if (data.y == 1) {
            output_grid.push_back(60u);
        } else if (data.y == 2) {
            output_grid.push_back(120u);
        } else if (data.y == 3) {
            output_grid.push_back(180u);
        } else {
            output_grid.push_back(255u);
        }
    }

    auto output_raster = nucleus::Raster<uint8_t>(nucleus::vector_layer::constants::grid_size, std::move(output_grid));
    auto image = nucleus::tile::conversion::u8raster_to_qimage(output_raster);
    image.save(QString("vector_layer_amounts.png"));
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

TEST_CASE("nucleus/vector_preprocess/clipping")
{
    SECTION("Clip to rect if all outside")
    {
        ClipperPaths shapes = { Clipper2Lib::MakePath<ClipperResolution>({ -10, -10, -10, 10, 10, 10, 10, -10 }) };

        ClipperRect rect = ClipperRect(-5, -5, 5, 5);
        ClipperPaths solution = RectClipFunc<ClipperResolution>(rect, shapes);

        CHECK(solution.size() == 1);
        CHECK(solution[0].size() == 4);
        CHECK(solution[0][0] == ClipperPoint { -5, -5 });
        CHECK(solution[0][1] == ClipperPoint { -5, 5 });
        CHECK(solution[0][2] == ClipperPoint { 5, 5 });
        CHECK(solution[0][3] == ClipperPoint { 5, -5 });
    }

    SECTION("Clip to rect no overlap")
    {
        ClipperPaths shapes = { Clipper2Lib::MakePath<ClipperResolution>({ -20, -20, -20, -10, -10, -10, -10, -20 }) };

        ClipperRect rect = ClipperRect(-5, -5, 5, 5);
        ClipperPaths solution = RectClipFunc<ClipperResolution>(rect, shapes);

        CHECK(solution.size() == 0);
    }

    SECTION("Clip diamond with rect")
    { // checks clipping against each edge
        ClipperPaths shapes = { Clipper2Lib::MakePath<ClipperResolution>({ 0, -7, -7, 0, 0, 7, 7, 0 }) };

        ClipperRect rect = ClipperRect(-5, -5, 5, 5);
        ClipperPaths solution = RectClipFunc<ClipperResolution>(rect, shapes);

        CHECK(solution.size() == 1);
        CHECK(solution[0].size() == 8);
        CHECK(solution[0][0] == ClipperPoint { 5, 2 });
        CHECK(solution[0][1] == ClipperPoint { 5, -2 });
        CHECK(solution[0][2] == ClipperPoint { 2, -5 });
        CHECK(solution[0][3] == ClipperPoint { -2, -5 });
        CHECK(solution[0][4] == ClipperPoint { -5, -2 });
        CHECK(solution[0][5] == ClipperPoint { -5, 2 });
        CHECK(solution[0][6] == ClipperPoint { -2, 5 });
        CHECK(solution[0][7] == ClipperPoint { 2, 5 });
    }

    SECTION("Clip poly with hole")
    {
        // polygon and hole are outside
        // -> will be cliped to same shape
        ClipperPaths shapes = { Clipper2Lib::MakePath<ClipperResolution>({ -15, -15, 15, -15, 15, 15, -15, 15 }),
            Clipper2Lib::MakePath<ClipperResolution>({ -10, -10, -10, 10, 10, 10, 10, -10 }) };

        ClipperRect rect = ClipperRect(-5, -5, 5, 5);
        ClipperPaths solution = RectClipFunc<ClipperResolution>(rect, shapes);

        // winding order is kept (points might be a bit rearanged though)
        CHECK(solution.size() == 2);
        CHECK(solution[0].size() == 4);
        CHECK(solution[0][0] == ClipperPoint { -5, 5 });
        CHECK(solution[0][1] == ClipperPoint { -5, -5 });
        CHECK(solution[0][2] == ClipperPoint { 5, -5 });
        CHECK(solution[0][3] == ClipperPoint { 5, 5 });
        CHECK(solution[1].size() == 4);
        CHECK(solution[1][0] == ClipperPoint { -5, -5 });
        CHECK(solution[1][1] == ClipperPoint { -5, 5 });
        CHECK(solution[1][2] == ClipperPoint { 5, 5 });
        CHECK(solution[1][3] == ClipperPoint { 5, -5 });

        // both polygons are visualizing the same. after rasterization we expect that no polygon is rasterized
        // -> since we are fully in a hole

        // convert clipping solution to vector that triangulize method can use
        std::vector<std::vector<glm::vec2>> clipped_poly;
        for (size_t i = 0; i < solution.size(); i++) {
            std::vector<glm::vec2> p;
            for (size_t j = 0; j < solution[i].size(); j++) {
                p.emplace_back(solution[i][j].x, solution[i][j].y);
            }
            clipped_poly.push_back(p);
        }

        std::vector<glm::vec2> triangle_points = nucleus::utils::rasterizer::triangulize(clipped_poly, true);

        // no triangles were created
        CHECK(triangle_points.size() == 0);
    }

    SECTION("Clip poly with diamond hole")
    {
        // polygon and hole are outside
        // -> will be cliped to same shape
        ClipperPaths shapes = { Clipper2Lib::MakePath<ClipperResolution>({ -15, -15, 15, -15, 15, 15, -15, 15 }),
            Clipper2Lib::MakePath<ClipperResolution>({ 0, -7, -7, 0, 0, 7, 7, 0 }) };

        ClipperRect rect = ClipperRect(-5, -5, 5, 5);
        ClipperPaths solution = RectClipFunc<ClipperResolution>(rect, shapes);

        // winding order is kept (points might be a bit rearanged though)
        CHECK(solution.size() == 2);
        CHECK(solution[0].size() == 4);
        CHECK(solution[0][0] == ClipperPoint { -5, 5 });
        CHECK(solution[0][1] == ClipperPoint { -5, -5 });
        CHECK(solution[0][2] == ClipperPoint { 5, -5 });
        CHECK(solution[0][3] == ClipperPoint { 5, 5 });
        CHECK(solution[1].size() == 8);
        CHECK(solution[1][0] == ClipperPoint { 5, 2 });
        CHECK(solution[1][1] == ClipperPoint { 5, -2 });
        CHECK(solution[1][2] == ClipperPoint { 2, -5 });
        CHECK(solution[1][3] == ClipperPoint { -2, -5 });
        CHECK(solution[1][4] == ClipperPoint { -5, -2 });
        CHECK(solution[1][5] == ClipperPoint { -5, 2 });
        CHECK(solution[1][6] == ClipperPoint { -2, 5 });
        CHECK(solution[1][7] == ClipperPoint { 2, 5 });

        // we expect that the rasterizer sees 4 triangles at the corners of the clip -> the diamond remains a hole

        // convert clipping solution to vector that triangulize method can use
        std::vector<std::vector<glm::vec2>> clipped_poly;
        for (size_t i = 0; i < solution.size(); i++) {
            std::vector<glm::vec2> p;
            for (size_t j = 0; j < solution[i].size(); j++) {
                p.emplace_back(solution[i][j].x, solution[i][j].y);
            }
            clipped_poly.push_back(p);
        }

        std::vector<glm::vec2> triangle_points = nucleus::utils::rasterizer::triangulize(clipped_poly, true);

        // 4 triangles with 3 points each
        CHECK(triangle_points.size() == 4 * 3);
        CHECK(triangle_points[0] == glm::vec2 { -5, -5 });
        CHECK(triangle_points[1] == glm::vec2 { -2, -5 });
        CHECK(triangle_points[2] == glm::vec2 { -5, -2 });

        CHECK(triangle_points[3] == glm::vec2 { -5, 2 });
        CHECK(triangle_points[4] == glm::vec2 { -2, 5 });
        CHECK(triangle_points[5] == glm::vec2 { -5, 5 });

        CHECK(triangle_points[6] == glm::vec2 { 5, -5 });
        CHECK(triangle_points[7] == glm::vec2 { 2, -5 });
        CHECK(triangle_points[8] == glm::vec2 { 5, -2 });

        CHECK(triangle_points[9] == glm::vec2 { 5, 2 });
        CHECK(triangle_points[10] == glm::vec2 { 5, 5 });
        CHECK(triangle_points[11] == glm::vec2 { 2, 5 });
    }

    SECTION("clip lines")
    {
        // 4 points forming 3 line segments
        // first line clipped on both ends
        // second line completely outside
        // third line outside to middle
        ClipperPaths in = { Clipper2Lib::MakePath<ClipperResolution>({ -15, -15, 15, 15, 30, 30, 0, 0 }) };
        ClipperRect rect = ClipperRect(-5, -5, 5, 5);

        auto clipper = RectClipLines({ rect });
        ClipperPaths solution = clipper.Execute(in);

        CHECK(solution.size() == 2);
        CHECK(solution[0].size() == 2);
        CHECK(solution[1].size() == 2);

        // first line segment
        CHECK(solution[0][0].x == -5);
        CHECK(solution[0][0].y == -5);
        CHECK(solution[0][1].x == 5);
        CHECK(solution[0][1].y == 5);

        // second line segment
        CHECK(solution[1][0].x == 5);
        CHECK(solution[1][0].y == 5);
        CHECK(solution[1][1].x == 0);
        CHECK(solution[1][1].y == 0);
    }

    SECTION("fully covers cell")
    {

        ClipperPaths shapes1 = { Clipper2Lib::MakePath<ClipperResolution>({ -10, -10, -10, 10, 10, 10, 10, -10 }) };
        ClipperPaths shapes2 = { Clipper2Lib::MakePath<ClipperResolution>({ -10, -10, -10, 10, 10, 10, 10, -10 }) };
        ClipperPaths shapes3
            = { Clipper2Lib::MakePath<ClipperResolution>({ 0, -7, -7, 0, 0, 7, 7, 0 }) }; // all are outside but in a diamond shape -> not fully
        ClipperPaths shapes4 = { Clipper2Lib::MakePath<ClipperResolution>({ -10, -10, -10, 10, 5, 15, 10, 10, 10, -10 }) };

        ClipperRect rect = ClipperRect(-5, -5, 5, 5);
        ClipperPaths solution1 = RectClipFunc<ClipperResolution>(rect, shapes1);
        ClipperPaths solution2 = RectClipFunc<ClipperResolution>(rect, shapes2);
        ClipperPaths solution3 = RectClipFunc<ClipperResolution>(rect, shapes3);
        ClipperPaths solution4 = RectClipFunc<ClipperResolution>(rect, shapes4);

        CHECK(nucleus::vector_layer::Preprocessor::fully_covers(solution1, rect));
        CHECK(nucleus::vector_layer::Preprocessor::fully_covers(solution2, rect));
        CHECK(!nucleus::vector_layer::Preprocessor::fully_covers(solution3, rect));
        CHECK(nucleus::vector_layer::Preprocessor::fully_covers(solution4, rect));
    }

    SECTION("line fully covers cell")
    {
        ClipperPaths shapes1 = { Clipper2Lib::MakePath<ClipperResolution>({ 0, -5, 0, 5 }) };
        ClipperPaths shapes2 = { Clipper2Lib::MakePath<ClipperResolution>({ 1, -5, 1, 5 }) }; // barely not fully covered with 7.5 width
        ClipperPaths shapes3 = { Clipper2Lib::MakePath<ClipperResolution>({ 0, -15, 0, -5 }) }; // one end point is barely at the edge of the rect

        ClipperRect rect = ClipperRect(-5, -5, 5, 5);

        // auto clipper = Clipper2Lib::RectClipLines64({ rect });
        // ClipperPaths solution1 = clipper.Execute(shapes1);
        // ClipperPaths solution2 = clipper.Execute(shapes2);

        // line_width is too small -> we do not even try
        CHECK(nucleus::vector_layer::Preprocessor::line_fully_covers(shapes1, 5.0, rect) == -1u);
        // barely fully covers
        CHECK(nucleus::vector_layer::Preprocessor::line_fully_covers(shapes1, 7.5, rect) == 0);
        // barely not fully covered since translated by 1 unit
        CHECK(nucleus::vector_layer::Preprocessor::line_fully_covers(shapes2, 7.5, rect) == -1u);

        CHECK(nucleus::vector_layer::Preprocessor::line_fully_covers(shapes3, 7.1 + 5.0, rect) == 0); // barely full cover
        CHECK(nucleus::vector_layer::Preprocessor::line_fully_covers(shapes3, 7.0 + 5.0, rect) == -1u); // barely outside full cover

        // CHECK(nucleus::vector_layer::Preprocessor::line_fully_covers(solution1, 5.0, rect));
    }

    SECTION("clipping vector tile to cell")
    { // real example
        constexpr size_t expected_process_amount = 134940;

        Style style(":/vectorlayerstyles/openstreetmap.json"); // 13
        // Style style(":/vectorlayerstyles/qwant.json"); // 9
        // Style style(":/vectorlayerstyles/osm-bright.json"); // 10
        style.load();

        auto id = nucleus::tile::Id { .zoom_level = 14, .coords = { 8936, 5681 }, .scheme = nucleus::tile::Scheme::SlippyMap };
        auto file = QFile(QString("%1%2").arg(ALP_TEST_DATA_DIR, "vector_layer/vectortile_benchmark_14_8936_5681.pbf"));
        file.open(QFile::ReadOnly);
        const auto bytes = file.readAll();

        Preprocessor preprocessor(std::move(style));

        auto tile_data = preprocessor.parse_tile(id, bytes);
        preprocessor.preprocess_geometry(tile_data, id.zoom_level);
        auto tile = preprocessor.create_gpu_tile();

        CHECK(preprocessor.processed_amount() == expected_process_amount);

        BENCHMARK("parse tile")
        {
            auto output = preprocessor.parse_tile(id, bytes);
            return output;
        };

        BENCHMARK("preprocess geometry")
        {
            preprocessor.preprocess_geometry(tile_data, id.zoom_level);
            return;
        };

        // BENCHMARK("create gpu tile")
        // {
        //     auto output = preprocessor.create_gpu_tile(temp_data);
        //     return output;
        // };

        BENCHMARK("complete preprocess")
        {
            auto output = preprocessor.preprocess(id, bytes);
            return output;
        };

        // for (int i = 0; i < 500; i++) {
        //     auto output = preprocessor.preprocess(id, bytes);
        // }
    }
}

TEST_CASE("nucleus/vector_preprocess")
{

    SECTION("Simplify styles")
    {
        {
            // only draw first found style (front to back rendering)
            std::vector<std::vector<glm::u32vec2>> styles;

            styles.push_back({
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
            });

            styles.push_back({ LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment() });

            auto style_buffer = Style::create_style_buffer_data(styles);

            std::vector<uint32_t> style_indices { 0, 1 };
            const auto simplified = nucleus::vector_layer::Style::simplify_styles(&style_indices, 15, style_buffer);

            CHECK(simplified.size() == 1);
            CHECK(simplified[0] == 1);
        }
        {

            // draw both styles
            std::vector<std::vector<glm::u32vec2>> styles;

            styles.push_back({ LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment() });

            styles.push_back({
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 200, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
            });

            auto style_buffer = Style::create_style_buffer_data(styles);

            std::vector<uint32_t> style_indices { 0, 1 };
            const auto simplified = nucleus::vector_layer::Style::simplify_styles(&style_indices, 15, style_buffer);

            CHECK(simplified.size() == 2);
            CHECK(simplified[0] == 1); // layer 1 first
        }

        {
            // width changed -> draw 1 than 3

            std::vector<std::vector<glm::u32vec2>> styles;

            styles.push_back({
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 10, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
            });

            styles.push_back({
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
            });

            styles.push_back({ LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment(),
                LayerStyle { 255, 0, 1 * nucleus::vector_layer::constants::style_precision, 1, false }.buffer_alignment() });

            auto style_buffer = Style::create_style_buffer_data(styles);
            std::vector<uint32_t> style_indices { 0, 1, 2 };
            const auto simplified = nucleus::vector_layer::Style::simplify_styles(&style_indices, 15, style_buffer);

            CHECK(simplified.size() == 2);
            CHECK(simplified[0] == 2);
            // style with index 1 is not used as style with index 0 fully covered it
            CHECK(simplified[1] == 0);
        }
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

    SECTION("Triangle Data packer")
    {

        auto a = glm::i64vec2(0b010100, 0b110100);
        auto b = glm::i64vec2(0b100101, 0b101101);
        auto c = glm::i64vec2(0b011011, 0b001110);

        uint16_t style = 343u;

        auto packed = nucleus::vector_layer::Preprocessor::pack_shader_data({ a, b, c, glm::bvec3(), style, true });
        auto unpacked = nucleus::vector_layer::Preprocessor::unpack_shader_data(packed);

        CHECK(a == glm::i64vec2(unpacked.a));
        CHECK(b == glm::i64vec2(unpacked.b));
        CHECK(c == glm::i64vec2(unpacked.c));
        CHECK(style == unpacked.style_index);
    }

    SECTION("Line Data packer")
    {

        auto a = glm::ivec2(0b010100, 0b110100);
        auto b = glm::ivec2(0b100101, 0b101101);
        // auto c = glm::ivec2(0, 0);

        uint16_t style = 646u;

        auto packed = nucleus::vector_layer::Preprocessor::pack_shader_data({ a, b, b, glm::bvec3(), style, false });
        auto unpacked = nucleus::vector_layer::Preprocessor::unpack_shader_data(packed);

        CHECK(a == glm::ivec2(unpacked.a));
        CHECK(b == glm::ivec2(unpacked.b));
        // CHECK(c == glm::ivec2(unpacked.c));
        CHECK(style == unpacked.style_index);
    }

    SECTION("std::map order behaviour")
    {
        // make sure that the map behaviour is consistent across platforms
        // mainly the keys are sorted correctly

        std::map<uint32_t, uint32_t> map;
        map[10] = 10;
        map[20] = 30;
        map[15] = 60;

        auto values = std::vector<uint32_t>();
        std::transform(map.begin(), map.end(), std::back_inserter(values), [](std::pair<uint32_t, uint32_t> pair) { return pair.second; });

        CHECK(values.size() == 3);
        CHECK(values[0] == 10);
        CHECK(values[1] == 60); // is inserted on second position
        CHECK(values[2] == 30);
    }

    SECTION("Polygon area calculation - clockwise polygon")
    {
        auto vertices = Clipper2Lib::MakePath<ClipperResolution>({ 2560, 3584, 3072, 3584, 3072, 4096 });

        CHECK(vertices.size() == 3);
        auto area = Preprocessor::polygon_area(vertices);
        CHECK(area > 0.0);
    }

    SECTION("Polygon area calculation - counter-clockwise polygon")
    {
        auto vertices = Clipper2Lib::MakePath<ClipperResolution>({ 3072, 3584, 2560, 3584, 3072, 4096 });

        CHECK(vertices.size() == 3);
        auto area = Preprocessor::polygon_area(vertices);
        CHECK(area < 0.0);
    }
}
