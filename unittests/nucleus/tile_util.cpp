/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2023 Adam Celarek
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

#include <QBuffer>
#include <QFile>
#include <QImage>
#include <QThread>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nucleus/camera/Definition.h>

#include "nucleus/camera/PositionStorage.h"
#include "nucleus/tile/utils.h"
#include "radix/quad_tree.h"

using Catch::Approx;
using namespace nucleus::tile;
namespace quad_tree = radix::quad_tree;
using nucleus::tile::utils::AabbDecorator;
using radix::TileHeights;

TEST_CASE("nucleus/tile/utils/make_bounds altitude correction")
{
    SECTION("root tile") {
        const auto bounds = utils::make_bounds(Id { 0, { 0, 0 } }, 100, 1000);
        CHECK(bounds.min.z == Approx(100 - 0.5).scale(100));
        CHECK(bounds.max.z > 1000);
    }
    SECTION("tile from equator") {
        auto bounds = utils::make_bounds(Id { 1, { 0, 1 } }, 100, 1000); // top left
        CHECK(bounds.min.z == Approx(100 - 0.5).scale(100));
        CHECK(bounds.max.z > 1000);
        bounds = utils::make_bounds(Id { 1, { 0, 0 } }, 100, 1000); // bottom left
        CHECK(bounds.min.z == Approx(100 - 0.5).scale(100));
        CHECK(bounds.max.z > 1000);
    }
}

TEST_CASE("tile/utils/refine_functor")
{
    // todo: optimise / benchmark refine functor
    // todo: optimise / benchmark draw list generator
    auto camera = nucleus::camera::stored_positions::stephansdom_closeup();

    QFile file(":/map/height_data.atb");
    const auto open = file.open(QIODeviceBase::OpenModeFlag::ReadOnly);
    assert(open);
    Q_UNUSED(open);
    const QByteArray data = file.readAll();
    const auto decorator = AabbDecorator::make(radix::TileHeights::deserialise(data));

    const auto refine_functor = utils::refineFunctor(camera, decorator, 256, 18);
    const auto all_leaves = quad_tree::onTheFlyTraverse(Id { 0, { 0, 0 } }, [](const Id& v) { return v.zoom_level < 6; }, [](const Id& v) { return v.children(); });

    BENCHMARK("refine functor double")
    {
        auto retval = true;
        for (const auto &id : all_leaves) {
            retval = retval && !refine_functor(id);
        }
        return retval;
    };

    const auto refine_functor_float = utils::refine_functor_float(camera, decorator, 256, 18);
    BENCHMARK("refine functor float")
    {
        auto retval = true;
        for (const auto &id : all_leaves) {
            retval = retval && !refine_functor_float(id);
        }
        return retval;
    };
}

TEST_CASE("tile/utils/refine_functor_max")
{
    // tests if the method to determine the max position in an aabb is correct
    auto camera = nucleus::camera::stored_positions::stephansdom_closeup();
    QFile file(":/map/height_data.atb");
    const auto open = file.open(QIODeviceBase::OpenModeFlag::ReadOnly);
    assert(open);
    Q_UNUSED(open);
    const QByteArray data = file.readAll();
    const auto decorator = AabbDecorator::make(radix::TileHeights::deserialise(data));

    // const auto refine_functor = utils::refineFunctor(camera, decorator, 1.0);
    const auto all_leaves = quad_tree::onTheFlyTraverse(Id { 0, { 0, 0 } }, [](const Id& v) { return v.zoom_level < 6; }, [](const Id& v) { return v.children(); });

    for (const auto& id : all_leaves) {

        const auto aabb = decorator->aabb(id);

        // code from refine_functor START
        const auto camera_position = camera.position();
        glm::dvec3 max_corner = { camera_position.x < (aabb.min.x + aabb.max.x) / 2.0 ? aabb.max.x : aabb.min.x,
            camera_position.y < (aabb.min.y + aabb.max.y) / 2.0 ? aabb.max.y : aabb.min.y,
            camera_position.z < (aabb.min.z + aabb.max.z) / 2.0 ? aabb.max.z : aabb.min.z };
        const auto delta = max_corner - camera_position;
        const auto distance = float(std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z));
        // code from refine_functor END

        std::array<glm::dvec3, 8> all_corners = { glm::dvec3(aabb.min.x, aabb.min.y, aabb.min.z),
            glm::dvec3(aabb.min.x, aabb.min.y, aabb.max.z),
            glm::dvec3(aabb.min.x, aabb.max.y, aabb.min.z),
            glm::dvec3(aabb.min.x, aabb.max.y, aabb.max.z),
            glm::dvec3(aabb.max.x, aabb.min.y, aabb.min.z),
            glm::dvec3(aabb.max.x, aabb.min.y, aabb.max.z),
            glm::dvec3(aabb.max.x, aabb.max.y, aabb.min.z),
            glm::dvec3(aabb.max.x, aabb.max.y, aabb.max.z) };

        for (const auto& corner : all_corners) {
            const auto curr_delta = corner - camera_position;
            const auto curr_distance = float(sqrt(curr_delta.x * curr_delta.x + curr_delta.y * curr_delta.y + curr_delta.z * curr_delta.z));

            CHECK(curr_distance <= distance);
        }
    }
}

TEST_CASE("tile/utils/camera_frustum_contains_tile")
{
    using nucleus::tile::utils::camera_frustum_contains_tile;
    using nucleus::tile::utils::camera_frustum_contains_tile_old;
    QFile file(":/map/height_data.atb");
    const auto open = file.open(QIODeviceBase::OpenModeFlag::ReadOnly);
    assert(open);
    Q_UNUSED(open);
    const QByteArray data = file.readAll();
    const auto decorator = AabbDecorator::make(radix::TileHeights::deserialise(data));

    SECTION("case 1")
    {
        auto cam = nucleus::camera::Definition(glm::dvec3 { 0, 0, 0 }, glm::dvec3 { 0, 1, 0 });
        cam.set_viewport_size({ 100, 100 });
        cam.set_field_of_view(90);
        cam.set_near_plane(0.01f);

        CHECK(camera_frustum_contains_tile(cam.frustum(), SrsAndHeightBounds { { -1., 9., -1. }, { 1., 10., 1. } }));
        CHECK(camera_frustum_contains_tile(cam.frustum(), SrsAndHeightBounds { { 0., 0., 0. }, { 1., 1., 1. } }));
        CHECK(camera_frustum_contains_tile(cam.frustum(), SrsAndHeightBounds { { -10., -10., -10. }, { 10., 1., 10. } }));
        CHECK(!camera_frustum_contains_tile(cam.frustum(), SrsAndHeightBounds { { -10., -10., -10. }, { 10., -1., 10. } }));
        CHECK(!camera_frustum_contains_tile(cam.frustum(), SrsAndHeightBounds { { -10., 0., -10. }, { -9., 1., -9. } }));
    }
    SECTION("case 2")
    {
        auto cam = nucleus::camera::stored_positions::wien();
        cam.set_viewport_size({ 1920, 1080 });
        cam.set_field_of_view(60);

        CHECK(!camera_frustum_contains_tile_old(cam.frustum(), SrsAndHeightBounds { { 1878516.4071364924, 5635549.221409474, 0.0 }, { 2504688.5428486564, 6261721.357121637, 1157.4707507122087 } }));
        CHECK(!camera_frustum_contains_tile(cam.frustum(), SrsAndHeightBounds { { 1878516.4071364924, 5635549.221409474, 0.0 }, { 2504688.5428486564, 6261721.357121637, 1157.4707507122087 } }));

        CHECK(camera_frustum_contains_tile(cam.frustum(), decorator->aabb({ 3, { 4, 4 } })) == camera_frustum_contains_tile_old(cam.frustum(), decorator->aabb({ 3, { 4, 4 } })));
    }
    SECTION("case 3")
    {
        auto cam = nucleus::camera::Definition({ 1.76106e+06, 6.07163e+06, 2510.08 },
            glm::dvec3 { 1.76106e+06, 6.07163e+06, 2510.08 } - glm::dvec3 { 0.9759, 0.19518, 0.09759 });
        cam.set_viewport_size({ 2561, 1369 });
        cam.set_field_of_view(60);
        CHECK(camera_frustum_contains_tile(cam.frustum(), decorator->aabb({ 10, { 557, 667 } })) == camera_frustum_contains_tile_old(cam.frustum(), decorator->aabb({ 10, { 557, 667 } })));
    }

    SECTION("many automated test cases")
    {
        std::vector<Id> tile_ids;
        quad_tree::onTheFlyTraverse(
            Id { 0, { 0, 0 } },
            [](const Id& v) { return v.zoom_level < 7; },
            [&tile_ids](const Id& v) {
                tile_ids.push_back(v);
                return v.children();
            });

        QFile file(":/map/height_data.atb");
        const auto open = file.open(QIODeviceBase::OpenModeFlag::ReadOnly);
        assert(open);
        Q_UNUSED(open);
        const QByteArray data = file.readAll();
        const auto decorator = AabbDecorator::make(TileHeights::deserialise(data));

        const auto camera_positions = std::vector {
            nucleus::camera::stored_positions::karwendel(),
            nucleus::camera::stored_positions::grossglockner(),
            nucleus::camera::stored_positions::oestl_hochgrubach_spitze(),
            nucleus::camera::stored_positions::schneeberg(),
            nucleus::camera::stored_positions::wien(),
            nucleus::camera::stored_positions::stephansdom(),
        };
        for (const auto& camera_position : camera_positions) {
            quad_tree::onTheFlyTraverse(Id { 0, { 0, 0 } }, nucleus::tile::utils::refineFunctor(camera_position, decorator, 256, 18), [&tile_ids](const Id& v) {
                tile_ids.push_back(v);
                return v.children();
            });
        }
        for (const auto& camera : camera_positions) {
            const auto camera_frustum = camera.frustum();
            for (const auto& tile_id : tile_ids) {
                const auto aabb = decorator->aabb(tile_id);
                CHECK(camera_frustum_contains_tile(camera_frustum, aabb) == camera_frustum_contains_tile_old(camera_frustum, aabb));
            }
        }

        BENCHMARK("camera_frustum_contains_tile")
        {
            bool retval = false;
            for (const auto& camera : camera_positions) {
                const auto camera_frustum = camera.frustum();
                for (const auto& tile_id : tile_ids) {
                    const auto aabb = decorator->aabb(tile_id);
                    retval = retval != camera_frustum_contains_tile(camera_frustum, aabb);
                }
            }
            return retval;
        };

        BENCHMARK("camera_frustum_contains_tile_old")
        {
            bool retval = false;
            for (const auto& camera : camera_positions) {
                const auto camera_frustum = camera.frustum();
                for (const auto& tile_id : tile_ids) {
                    const auto aabb = decorator->aabb(tile_id);
                    retval = retval != camera_frustum_contains_tile_old(camera_frustum, aabb);
                }
            }
            return retval;
        };
    }
}
