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

// nucleus::Raster<uint8_t> visualize_grid(const nucleus::vector_layer::VectorLayerGrid& processed_grid, int grid_width)
// {
//     std::vector<uint8_t> output_grid;

//     for (size_t i = 0; i < processed_grid.size(); ++i) {

//         // grid is a singular uint32_t value that encodes the start index of the triangle list and the amount of triangles
//         if (processed_grid[i].size() == 0) {
//             output_grid.push_back(0);
//         } else {
//             output_grid.push_back(255u);
//         }
//     }

//     return nucleus::Raster<uint8_t>(grid_width, std::move(output_grid));
// }

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

// TODOs

// fully filled cell
// - vector<Clipper2Lib::RectClip64> grid to vector<struct{is_done, Clipper2Lib::RectClip64}> grid
// - is_done is set to true, if polygon fully covers cell
// - beware of transparancies

// polylines
// - aabb needs to expand by 1/2 * line width in all directions
// - clipper rects need to expand by 1/2 + line width in all directions
// - do clipper rects need to be uniformely sized? -> in theory i think not since the origin stays the same (and we use negative numbers)
// -- but we would need to create new rectclip objects for each individual line width (which is expensive)
// -> therefore one oversized rectclip with some kind of upper boundary of how large lines could be
// -> but the larger this is, the more cells will be filled unintentionally
// - alternative -> create two boundary polylines (Clipper2Lib::InflatePaths ?)
// - only if the boundary line falls into the inner clip, the outer clip will be executed

// triangles/polylines coordinates are anchored to cell origin

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

    // SECTION("clip test")
    // {
    //     // clipper2 result that is outside rect
    //     ClipperPaths shapes = { Clipper2Lib::MakePath<ClipperResolution>({ 128, -64, 91, 150, 77, 150, 114, -64 }) };

    //     ClipperRect rect = ClipperRect(64, 64, 127, 127);
    //     ClipperPaths solution = RectClip(rect, shapes);

    //     qDebug() << "points:";
    //     for (size_t i = 0; i < solution.size(); ++i) {
    //         for (size_t j = 0; j < solution[i].size(); ++j) {
    //             qDebug() << solution[i][j].x << solution[i][j].y;
    //         }
    //     }

    //     // output
    //     // 91 63
    //     // 105 63
    //     // 94 127
    //     // 80 127
    // }

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

    // SECTION("preprocess simple tile")
    // {
    //     // the main reason for this test case is to check if the amount of data within a certain cell is correct
    //     // we look at a tile on a mountain at a cell where only one polygon should be present

    //     Style style(":/vectorlayerstyles/openstreetmap.json");
    //     style.load();

    //     auto id = nucleus::tile::Id { .zoom_level = 14, .coords = { 8781, 5760 }, .scheme = nucleus::tile::Scheme::SlippyMap };
    //     auto file = QFile(QString("%1%2").arg(ALP_TEST_DATA_DIR, "vector_layer/vectortile_mountain_14_8781_5760.pbf"));
    //     file.open(QFile::ReadOnly);
    //     const auto bytes = file.readAll();

    //     // const auto style_buffer = style.styles()->buffer();

    //     Preprocessor preprocessor(std::move(style));

    //     auto tile_data = preprocessor.parse_tile(id, bytes);
    //     preprocessor.preprocess_geometry(tile_data);
    //     auto tile = preprocessor.create_gpu_tile();

    //     const auto& pixel = tile.acceleration_grid->pixel({ 40, 40 });
    //     const auto data = nucleus::utils::bit_coding::u32_to_u24_u8(pixel);

    //     // investigate geometry in detail
    //     // const auto geom1 = Preprocessor::unpack_data(tile.geometry_buffer->buffer()[data.x]);
    //     // const auto geom2 = nucleus::vector_layer::Preprocessor::unpack_data(tile.geometry_buffer->buffer()[data.x + 1]);
    //     // qDebug() << geom1.a.x << geom1.a.y << geom1.b.x << geom1.b.y << geom1.c.x << geom1.c.y << geom1.style_index;
    //     // qDebug() << geom2.a.x << geom2.a.y << geom2.b.x << geom2.b.y << geom2.c.x << geom2.c.y << geom2.style_index;

    //     // make sure that the amount of geometry at the cell is exactly 2
    //     // 2 because we triangulate a square to two triangles
    //     CHECK(data.y == 2);

    //     // visualize_acceleration_grid(tile.acceleration_grid);
    // }

    SECTION("clipping vector tile to cell")
    { // real example
        constexpr size_t expected_process_amount = 128006;

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
        preprocessor.preprocess_geometry(tile_data);
        auto tile = preprocessor.create_gpu_tile();

        CHECK(preprocessor.processed_amount() == expected_process_amount);

        BENCHMARK("parse tile")
        {
            auto output = preprocessor.parse_tile(id, bytes);
            return output;
        };

        BENCHMARK("preprocess geometry")
        {
            preprocessor.preprocess_geometry(tile_data);
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

    // SECTION("Tile download basemap")
    // {
    //     // if this fails it is very likely that something on the vector tile server changed
    //     // manually download the tile from the below link and check if the changes are valid and replace vectortile.mvt with this new file
    //     // https://osm.cg.tuwien.ac.at/vector_tiles/poi_v1/10/548/359

    //     const auto id = nucleus::tile::Id { .zoom_level = 10, .coords = { 548, 359 }, .scheme = nucleus::tile::Scheme::SlippyMap };
    //     nucleus::tile::TileLoadService service("https://mapsneu.wien.gv.at/basemapv/bmapv/3857/tile/",
    //     nucleus::tile::TileLoadService::UrlPattern::ZYX_yPointingSouth, ".pbf");

    //     {
    //         QSignalSpy spy(&service, &nucleus::tile::TileLoadService::load_finished);
    //         service.load(id);
    //         spy.wait(15000);

    //         REQUIRE(spy.count() == 1);
    //         QList<QVariant> arguments = spy.takeFirst();
    //         REQUIRE(arguments.size() == 1);
    //         nucleus::tile::Data tile = arguments.at(0).value<nucleus::tile::Data>();
    //         CHECK(tile.id == id);
    //         CHECK(tile.network_info.status == nucleus::tile::NetworkInfo::Status::Good);
    //         CHECK(nucleus::utils::time_since_epoch() - tile.network_info.timestamp < 10'000);

    //         REQUIRE(tile.data->size() > 0);
    //         CHECK(tile.data->size() > 2000);
    //     }
    // }

    // SECTION("Triangle to Grid")
    // { // TODO redo after refactor
    //     const std::vector<glm::u32vec4> style_buffer = { glm::u32vec4(0, 0, 0, 0) };

    //     constexpr auto extent = 64u;
    //     const std::vector<std::vector<glm::vec2>> triangle_left_hypo = { { glm::vec2(10, 30), glm::vec2(30, 5), glm::vec2(50, 50) } };
    //     const std::vector<std::vector<glm::vec2>> triangle_right_hypo = { { glm::vec2(5, 5), glm::vec2(25, 10), glm::vec2(5, 15) } };

    //     std::vector<GeometryData> tile_data { { triangle_left_hypo, extent, { { 0, 0 } }, true }, { triangle_right_hypo,
    //     extent, { { 0, 1 } }, true } };

    //     auto processed = nucleus::vector_layer::details::preprocess_geometry(tile_data, style_buffer);

    //     auto raster = visualize_grid(processed, nucleus::vector_layer::constants::grid_size);
    //     auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

    //     auto test_image = example_grid_data_triangles();
    //     CHECK(image == test_image);

    //     auto tile = nucleus::vector_layer::details::create_gpu_tile(processed);

    //     CHECK(same_cells_are_filled(raster, tile.acceleration_grid));

    //     auto data = tile.vertex_buffer->buffer();
    //     REQUIRE(data.size() == nucleus::vector_layer::constants::data_size[tile.buffer_info] *
    //     nucleus::vector_layer::constants::data_size[tile.buffer_info]); CHECK(data[0].x == 1081671760); CHECK(data[0].y == 1076429280); CHECK(data[0].z ==
    //     1086915360); CHECK(data[1].x == 1075118160); CHECK(data[1].y == 1080361120); CHECK(data[1].z == 1075118320);
    //     // the rest should be undefined -> -1u
    //     CHECK(data[2].x == -1u);
    //     CHECK(data[3].y == -1u);
    //     CHECK(data[nucleus::vector_layer::constants::data_size[tile.buffer_info] * 1.7].x == -1u);

    //     // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
    //     //     std::cout << index_buffer[i] << std::endl;
    //     // }

    //     // std::cout << std::endl;
    //     // for (int i = 0; i < 14; ++i) { // DEBUG expected triangle data
    //     //     std::cout << data[i] << std::endl;
    //     // }

    //     // // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
    //     // image.save(QString("vector_layer_grid_triangles.png"));
    // }

    // SECTION("Triangle to Grid (outside vertices)")
    // { // TODO redo after refactor
    //     // make sure that the proposed points are still roughly the same (even after tinkering with grid_size)
    //     // this would still mean that we have to redo the below data every time we chang the grid scale -> but it isn't too problematic for now
    //     constexpr auto extent = 64u;

    //     const std::vector<glm::u32vec4> style_buffer = { glm::u32vec4(0, 0, 0, 0) };

    //     const std::vector<std::vector<glm::vec2>> triangle_left_hypo = { { glm::vec2(-10, 30), glm::vec2(30, 5), glm::vec2(50, 80) } };
    //     const std::vector<std::vector<glm::vec2>> triangle_right_hypo = { { glm::vec2(5, -5), glm::vec2(90, 10), glm::vec2(5, 15) } };

    //     std::vector<nucleus::vector_layer::details::GeometryData> tile_data { { triangle_left_hypo, extent, { { 0, 0 } }, true }, { triangle_right_hypo,
    //     extent, { { 0, 1 } }, true } };

    //     auto processed = nucleus::vector_layer::details::preprocess_geometry(tile_data, style_buffer);

    //     auto raster = visualize_grid(processed, nucleus::vector_layer::constants::grid_size);
    //     auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

    //     auto tile = nucleus::vector_layer::details::create_gpu_tile(processed);

    //     CHECK(same_cells_are_filled(raster, tile.acceleration_grid));

    //     auto data = tile.vertex_buffer->buffer();
    //     REQUIRE(data.size() == nucleus::vector_layer::constants::data_size[tile.buffer_info] *
    //     nucleus::vector_layer::constants::data_size[tile.buffer_info]); CHECK(data[0].x == 1081671760); CHECK(data[0].y == 1071186400); CHECK(data[0].z ==
    //     1086915840); CHECK(data[1].x == 1075118000); CHECK(data[1].y == 1097400480); CHECK(data[1].z == 1075118320);
    //     // the rest should be undefined -> -1u
    //     CHECK(data[2].x == -1u);
    //     CHECK(data[3].y == -1u);

    //     CHECK(data[nucleus::vector_layer::constants::data_size[tile.buffer_info] * 1.7].x == -1u);

    //     // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
    //     //     std::cout << index_buffer[i] << std::endl;
    //     // }

    //     // std::cout << std::endl;
    //     // for (int i = 0; i < 14; ++i) { // DEBUG expected triangle data
    //     //     std::cout << data[i] << std::endl;
    //     // }

    //     // // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
    //     // image.save(QString("vector_layer_grid_triangles_outside.png"));
    // }

    // SECTION("Line to Grid")
    // { // TODO redo after refactor
    //     const std::vector<glm::u32vec4> style_buffer = { glm::u32vec4(0, 0, 0, 0), glm::u32vec4(1, 0, 0, 0) };

    //     constexpr auto extent = 64;
    //     const std::vector<std::vector<glm::vec2>> line0 = { { glm::vec2(10, 30), glm::vec2(30, 50), glm::vec2(50, 30) } };

    //     std::vector<nucleus::vector_layer::details::GeometryData> tile_data { { line0, extent, { { 0, 0 }, { 1, 1 } }, false } };

    //     auto processed = nucleus::vector_layer::details::preprocess_geometry(tile_data, style_buffer);

    //     auto raster = visualize_grid(processed, nucleus::vector_layer::constants::grid_size);
    //     auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);

    //     // auto test_image = example_grid_data_triangles();
    //     // CHECK(image == test_image);

    //     auto tile = nucleus::vector_layer::details::create_gpu_tile(processed);

    //     CHECK(same_cells_are_filled(raster, tile.acceleration_grid));

    //     // we provide two line segments that overlap at one point.
    //     // the line has two "different" styles on two different layers
    //     // we expect that the index buffer has indices 1) with only start 2) with only end 3) with both start and end indices
    //     // furthermore we expect that layer with higher layer indices are stored earlier

    //     auto data = tile.vertex_buffer->buffer();
    //     REQUIRE(data.size() == nucleus::vector_layer::constants::data_size[tile.buffer_info] *
    //     nucleus::vector_layer::constants::data_size[tile.buffer_info]); CHECK(data[0].x == 1076429280); CHECK(data[0].y == 1081672480); CHECK(data[0].z ==
    //     1073807360); CHECK(data[1].x == 1081672480); CHECK(data[1].y == 1086915040); CHECK(data[1].z == 1073807360);
    //     // the rest should be undefined -> -1u
    //     CHECK(data[2].x == -1u);
    //     CHECK(data[3].y == -1u);
    //     CHECK(data[nucleus::vector_layer::constants::data_size[tile.buffer_info] * 1.7].x == -1u);

    //     // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
    //     //     std::cout << index_buffer[i] << std::endl;
    //     // }

    //     // std::cout << std::endl;
    //     // for (int i = 0; i < 14; ++i) { // DEBUG expected triangle data
    //     //     std::cout << data[i] << std::endl;
    //     // }

    //     // // DEBUG: save image (image saved to build/Desktop-Profile/unittests/nucleus)
    //     image.save(QString("vector_layer_grid_lines.png"));
    // }
    SECTION("Simplify styles")
    {
        {
            // only draw second style (passed by style_indices)
            // style width index 0 and 2 in style buffer are here for blending between zoom steps
            std::vector<glm::u32vec4> style_buffer { { 200, 0, 0, 0 }, { 200, 0, 0, 0 }, { 255, 0, 0, 0 }, { 255, 0, 0, 0 } };
            std::vector<StyleLayerIndex> style_indices { { 1, 1 }, { 3, 3 } };
            const auto simplified = nucleus::vector_layer::Preprocessor::simplify_styles(&style_indices, style_buffer);

            CHECK(simplified.size() == 1);
            CHECK(simplified[0].style_index == 3);
        }
        {
            // draw both styles
            // style width index 0 in style buffer is here for blending between zoom steps
            std::vector<glm::u32vec4> style_buffer { { 200, 0, 0, 0 }, { 200, 0, 0, 0 }, { 200, 0, 0, 0 } };
            std::vector<StyleLayerIndex> style_indices { { 1, 1 }, { 2, 2 } };
            const auto simplified = nucleus::vector_layer::Preprocessor::simplify_styles(&style_indices, style_buffer);

            CHECK(simplified.size() == 2);
            CHECK(simplified[0].style_index == 2); // but layer 1 first
        }

        {
            // width changed -> draw 3 than 1
            // styles width index 0 and 2 in style buffer are here for blending between zoom steps
            std::vector<glm::u32vec4> style_buffer { { 255, 0, 10, 0 }, { 255, 0, 10, 0 }, { 255, 0, 0, 0 }, { 255, 0, 0, 0 }, { 255, 0, 0, 0 } };
            std::vector<StyleLayerIndex> style_indices { { 1, 1 }, { 3, 3 }, { 4, 4 } };
            const auto simplified = nucleus::vector_layer::Preprocessor::simplify_styles(&style_indices, style_buffer);

            CHECK(simplified.size() == 2);
            CHECK(simplified[0].style_index == 4);
            // style with index 3 is not used as style with index 4 fully covered it
            CHECK(simplified[1].style_index == 1);
        }
    }

    // SECTION("Tile exploration") // section mostly used for tile debugging -> not a real test
    // {
    //     const auto id = nucleus::tile::Id { .zoom_level = 14, .coords = { 4477 * 2, 2850 * 2 + 1 }, .scheme = nucleus::tile::Scheme::SlippyMap };
    //     // const auto id = nucleus::tile::Id { .zoom_level = 13, .coords = { 4477, 2850 }, .scheme = nucleus::tile::Scheme::SlippyMap };
    //     // const auto id = nucleus::tile::Id { .zoom_level = 15, .coords = { (4477 * 2 + 1) * 2, (2850 * 2 + 1) * 2 }, .scheme = nucleus::tile::Scheme::SlippyMap };
    //     // const auto id = nucleus::tile::Id { .zoom_level = 19, .coords = { 285987, 181795 }, .scheme = nucleus::tile::Scheme::SlippyMap };

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

    //     size_t data_offset = 0;
    //     auto acceleration_grid = std::vector<std::set<uint32_t>>(nucleus::vector_layer::constants::grid_size * nucleus::vector_layer::constants::grid_size, std::set<uint32_t>());

    //     for (size_t i = 0; i < tile_data.size(); ++i) {
    //         if (tile_data[i].is_polygon) {

    //             // if (tile_data[i].style != 88u)
    //             //     continue;

    //             std::vector<glm::vec2> triangle_points = nucleus::utils::rasterizer::triangulize(tile_data[i].vertices, true);

    //             std::cout << std::endl << "o pedestrian_area" << i << std::endl;
    //             for (size_t j = 0; j < triangle_points.size() / 3; ++j) {
    //                 std::cout << "v " << triangle_points[j * 3 + 0].x << " 0 " << triangle_points[j * 3 + 0].y << std::endl;
    //                 std::cout << "v " << triangle_points[j * 3 + 1].x << " 0 " << triangle_points[j * 3 + 1].y << std::endl;
    //                 std::cout << "v " << triangle_points[j * 3 + 2].x << " 0 " << triangle_points[j * 3 + 2].y << std::endl;

    //                 auto f_ind = ((j + data_offset) * 3) + 1;
    //                 std::cout << "f " << f_ind << " " << (f_ind + 1) << " " << (f_ind + 2) << std::endl;
    //             }

    //             data_offset += triangle_points.size() / 3;

    //         } else {
    //             // continue;
    //             // qDebug() << tile_data[i].style;
    //             // if (tile_data[i].style != 92u)
    //             //     continue;
    //             // const auto cell_writer = [&acceleration_grid, data_offset](glm::vec2 pos, int data_index) {
    //             //     // if in grid_size bounds and not already present -> than add index to vector
    //             //     if (glm::all(glm::lessThanEqual({ 0, 0 }, pos)) && glm::all(glm::greaterThan(glm::vec2(nucleus::vector_layer::constants::grid_size), pos))) {
    //             //         // last bit of index indicates that this is a line
    //             //         // !! IMPORTANT !! we have to use the lowest bit since set orders the input depending on key -> lines and polygons have to stay intermixed
    //             //         const auto index = ((data_index + data_offset) << 1); // | 0u;
    //             //         acceleration_grid[int(pos.x) + nucleus::vector_layer::constants::grid_size * int(pos.y)].insert(index);
    //             //     }
    //             // };

    //             // std::cout << std::endl << "o road" << i << std::endl;
    //             // for (size_t j = 0; j < tile_data[i].vertices.size(); ++j) {
    //             //     std::cout << "v " << tile_data[i].vertices[j].x << " 0 " << tile_data[i].vertices[j].y << std::endl;
    //             // }
    //             // std::cout << "l ";
    //             // for (size_t j = 0; j < tile_data[i].vertices.size(); ++j) {
    //             //     std::cout << (j + data_offset) << " ";
    //             // }

    //             // std::cout << std::endl << std::endl;

    //             // const auto scale = float(nucleus::vector_layer::constants::grid_size) / float(tile_data[i].extent);
    //             // nucleus::utils::rasterizer::rasterize_line(cell_writer, tile_data[i].vertices, tile_data[i].line_width * scale, scale);

    //             // data_offset += tile_data[i].vertices.size();
    //         }
    //     }

    //     auto raster = visualize_grid(acceleration_grid, nucleus::vector_layer::constants::grid_size);
    //     auto image = nucleus::tile::conversion::u8raster_to_qimage(raster);
    //     image.save(QString("vector_layer_debuggggg.png"));

    //     // for (int i = 0; i < 10; ++i) { // DEBUG expected bridge data
    //     //     std::cout << index_buffer[i] << std::endl;
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

        auto a = glm::i64vec2(0b010100, 0b110100);
        auto b = glm::i64vec2(0b100101, 0b101101);
        auto c = glm::i64vec2(0b011011, 0b001110);

        uint16_t style = 343u;

        auto packed = nucleus::vector_layer::Preprocessor::pack_triangle_data({ a, b, c, style, true });
        auto unpacked = nucleus::vector_layer::Preprocessor::unpack_data(packed);

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

        auto packed = nucleus::vector_layer::Preprocessor::pack_line_data(a, b, style);
        auto unpacked = nucleus::vector_layer::Preprocessor::unpack_data(packed);

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
}

// TEST_CASE("nucleus/vector_preprocess benchmarks")
// {
//     // load tile data
//     // the normal tile is a tile in a small city, with more than half of the tile consisting of a mountain.
//     // the zoom level of 13 and 14 are the tiles with the most amount of data present, but since most of this tile is in the countryside, it only contains
//     68kb
//     // of data
//     // TODO this vector tile uses basemap -> it is not comparable for the benchmark
//     auto id_normal = nucleus::tile::Id { .zoom_level = 13, .coords = { 4412, 2893 }, .scheme = nucleus::tile::Scheme::SlippyMap };
//     auto file_normal = QFile(QString("%1%2").arg(ALP_TEST_DATA_DIR, "vector_layer/vectortile_13_4412_2893.pbf"));
//     file_normal.open(QFile::ReadOnly);
//     const auto bytes_normal = file_normal.readAll();

//     // the "worst" tile is one of the largest tile with 868kb (right over vienna -> lots of buildings and other details)
//     // so the output of the benchmark can be regarded as a worst case approximation
//     auto id_worst = nucleus::tile::Id { .zoom_level = 14, .coords = { 8936, 5681 }, .scheme = nucleus::tile::Scheme::SlippyMap };
//     auto file_worst = QFile(QString("%1%2").arg(ALP_TEST_DATA_DIR, "vector_layer/vectortile_benchmark_14_8936_5681.pbf"));
//     file_worst.open(QFile::ReadOnly);
//     const auto bytes_worst = file_worst.readAll();

//     // load style
//     Style style(":/vectorlayerstyles/openstreetmap.json");
//     style.load();

//     BENCHMARK("preprocess normal tile") { nucleus::vector_layer::preprocess(id_normal, bytes_normal, style); };
//     BENCHMARK("preprocess worse case tile") { nucleus::vector_layer::preprocess(id_worst, bytes_worst, style); };

//     // BENCHMARK("triangulize polygons")
//     // {
//     //     const std::vector<glm::vec2> polygon_points = { glm::vec2(10.5, 10.5), glm::vec2(30.5, 10.5), glm::vec2(50.5, 50.5), glm::vec2(10.5, 30.5) };
//     //     const auto edges = nucleus::utils::rasterizer::generate_neighbour_edges(polygon_points);
//     //     nucleus::utils::rasterizer::triangulize(polygon_points, edges);
//     // };
// }
