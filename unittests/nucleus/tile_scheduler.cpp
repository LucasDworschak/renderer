/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2023 Adam Celarek
 * Copyright (C) 2024 Gerald Kimmersdorfer
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

#include <unordered_set>

#include "nucleus/utils/Stopwatch.h"
#include "test_helpers.h"
#include <QImage>
#include <QSignalSpy>
#include <QThread>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nucleus/camera/PositionStorage.h>
#include <nucleus/tile/SchedulerDirector.h>
#include <nucleus/tile/TextureScheduler.h>
#include <nucleus/tile/conversion.h>
#include <nucleus/tile/types.h>
#include <nucleus/tile/utils.h>
#include <nucleus/utils/image_loader.h>
#include <radix/TileHeights.h>
#include <radix/tile.h>

using nucleus::tile::utils::AabbDecorator;
using radix::TileHeights;
using namespace nucleus::tile;

namespace {

#ifdef _MSC_VER
constexpr auto timing_multiplicator = 20;
#elif defined __ANDROID__
constexpr auto timing_multiplicator = 1;
#else
constexpr auto timing_multiplicator = 1;
#endif
std::unique_ptr<Scheduler> scheduler_with_true_heights()
{
    auto scheduler = std::make_unique<TextureScheduler>(Scheduler::Settings {});
    QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);

    QFile file(":/map/height_data.atb");
    const auto open = file.open(QIODeviceBase::OpenModeFlag::ReadOnly);
    assert(open);
    Q_UNUSED(open);
    const QByteArray data = file.readAll();
    const auto decorator = AabbDecorator::make(TileHeights::deserialise(data));
    scheduler->set_aabb_decorator(decorator);
    scheduler->set_update_timeout(0);
    scheduler->set_enabled(true);
    spy.wait(5 * timing_multiplicator); // wait for quad requests triggered by set_enabled
    REQUIRE(spy.size() == 1);
    // 'disable' timer updates. most tests trigger them automatically, the others will reset it.
    scheduler->set_update_timeout(1'000'000);
    return scheduler;
}

std::unique_ptr<TextureScheduler> default_scheduler()
{
    auto scheduler = std::make_unique<TextureScheduler>(Scheduler::Settings {});
    scheduler->set_name("test");
    QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);
    TileHeights h;
    h.emplace({ 0, { 0, 0 } }, { 100, 4000 });
    scheduler->set_aabb_decorator(AabbDecorator::make(std::move(h)));
    scheduler->set_update_timeout(1);
    scheduler->set_enabled(true);
    spy.wait(2 * timing_multiplicator); // wait for quad requests triggered by set_enabled
    REQUIRE(spy.size() == 1);
    // 'disable' timer updates. most tests trigger them automatically, the others will reset it.
    scheduler->set_update_timeout(1'000'000);
    return scheduler;
}

std::unique_ptr<Scheduler> scheduler_with_disk_cache()
{
    auto sch = default_scheduler();
    CHECK(sch->read_disk_cache());
    return sch;
}

std::unique_ptr<Scheduler> scheduler_with_aabb()
{
    auto scheduler = std::make_unique<TextureScheduler>(Scheduler::Settings {});
    std::filesystem::remove_all(scheduler->disk_cache_path());
    TileHeights h;
    h.emplace({ 0, { 0, 0 } }, { 100, 4000 });
    scheduler->set_aabb_decorator(AabbDecorator::make(std::move(h)));
    return scheduler;
}

QByteArray example_tile_data()
{
    auto height_file = QFile(QString("%1%2").arg(ALP_TEST_DATA_DIR, "test-tile_ortho.jpeg"));
    height_file.open(QFile::ReadOnly);
    const auto height_bytes = height_file.readAll();
    REQUIRE(!QImage::fromData(height_bytes).isNull());
    return height_bytes;
}

nucleus::tile::DataQuad example_tile_quad_for(const Id& id, unsigned n_children = 4, NetworkInfo::Status status = NetworkInfo::Status::Good)
{
    const auto children = id.children();
    REQUIRE(n_children <= 4);
    nucleus::tile::DataQuad cpu_quad;
    cpu_quad.id = id;
    cpu_quad.n_tiles = n_children;
    static const auto example_data = example_tile_data();
    for (unsigned i = 0; i < n_children; ++i) {
        cpu_quad.tiles[i].id = children[i];
        cpu_quad.tiles[i].data = std::make_shared<QByteArray>(example_data);
        cpu_quad.tiles[i].network_info.status = status;
        cpu_quad.tiles[i].network_info.timestamp = nucleus::utils::time_since_epoch();
    }
    return cpu_quad;
}

std::vector<nucleus::tile::DataQuad> example_quads_for(const std::vector<Id>& ids)
{
    std::vector<nucleus::tile::DataQuad> quads;
    quads.reserve(ids.size());
    for (const auto& id : ids) {
        quads.push_back(example_tile_quad_for(id));
    }

    return quads;
}

std::vector<nucleus::tile::DataQuad> example_quads_for_steffl_and_gg()
{
    static std::vector<nucleus::tile::DataQuad> retval = example_quads_for({
        Id { 0, { 0, 0 } },
        Id { 1, { 1, 1 } },
        Id { 2, { 2, 2 } },
        Id { 3, { 4, 5 } },
        Id { 4, { 8, 10 } },
        Id { 5, { 17, 20 } },
        Id { 6, { 34, 41 } },
        Id { 7, { 69, 83 } }, // stephans dom
        Id { 8, { 139, 167 } },
        Id { 9, { 279, 334 } },
        Id { 10, { 558, 668 } },
        Id { 10, { 558, 669 } },
        Id { 11, { 1117, 1337 } },
        Id { 11, { 1117, 1338 } },
        Id { 11, { 1116, 1337 } },
        Id { 11, { 1116, 1338 } },
        Id { 12, { 2234, 2675 } },
        Id { 7, { 68, 83 } }, // grossglockner
        Id { 7, { 68, 82 } },
        Id { 8, { 136, 166 } },
        Id { 8, { 137, 166 } },
        Id { 8, { 136, 165 } },
        Id { 8, { 137, 165 } },
        Id { 9, { 273, 332 } },
        Id { 9, { 274, 332 } },
        Id { 9, { 273, 331 } },
        Id { 9, { 274, 331 } },
        Id { 10, { 547, 664 } },
        Id { 10, { 548, 664 } },
        Id { 11, { 1095, 1328 } },
        Id { 11, { 1096, 1328 } },
        Id { 12, { 2191, 2657 } },
        Id { 12, { 2192, 2657 } },
        Id { 12, { 2191, 2656 } },
        Id { 12, { 2192, 2656 } },
        Id { 13, { 4384, 5313 } },
        Id { 13, { 4385, 5313 } },
        Id { 13, { 4384, 5312 } },
        Id { 13, { 4385, 5312 } },
    });
    return retval;
}
std::vector<nucleus::tile::DataQuad> example_quads_many()
{
    static std::vector<nucleus::tile::DataQuad> retval = []() {
        auto scheduler = default_scheduler();
        QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);
        auto camera = nucleus::camera::stored_positions::grossglockner();
        camera.set_viewport_size({ 3840, 2160 });
        scheduler->update_camera(camera);
        scheduler->send_quad_requests();
        REQUIRE(spy.size() == 1);
        return example_quads_for(spy.front().front().value<std::vector<Id>>());
    }();
    return retval;
}

} // namespace

TEST_CASE("nucleus/tile/Scheduler")
{
    SECTION("enable / disable tile scheduler")
    {
        auto scheduler = scheduler_with_aabb();
        // disabled by default (enable by signal, in order to prevent a race with opengl)
        scheduler->set_update_timeout(1);

        QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        spy.wait(2 * timing_multiplicator);
        CHECK(spy.empty());

        scheduler->set_enabled(true);
        spy.wait(2 * timing_multiplicator);
        CHECK(!spy.empty());
    }

#ifndef __EMSCRIPTEN__
    SECTION("update timeout & camera updates are collected")
    {
        auto scheduler = default_scheduler();
        scheduler->set_enabled(true);
        scheduler->set_update_timeout(5 * timing_multiplicator);

        QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        scheduler->update_camera(nucleus::camera::stored_positions::oestl_hochgrubach_spitze());
        test_helpers::process_events_for(2 * timing_multiplicator);
        CHECK(spy.empty());

        test_helpers::process_events_for(5 * timing_multiplicator);
        CHECK(spy.size() == 1);
    }

    SECTION("timer is not restarted on camera updates && stops if nothing happens")
    {
        auto scheduler = default_scheduler();
        scheduler->set_update_timeout(6 * timing_multiplicator);
        QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        test_helpers::process_events_for(2 * timing_multiplicator);
        scheduler->update_camera(nucleus::camera::stored_positions::oestl_hochgrubach_spitze());
        test_helpers::process_events_for(2 * timing_multiplicator);
        CHECK(spy.empty());
        scheduler->update_camera(nucleus::camera::stored_positions::grossglockner());
        test_helpers::process_events_for(3 * timing_multiplicator);
        CHECK(spy.size() == 1);
        test_helpers::process_events_for(7 * timing_multiplicator);
        CHECK(spy.size() == 1);
    }
#endif

    SECTION("quads are being requested")
    {
        auto scheduler = default_scheduler();
        QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        scheduler->send_quad_requests();
        REQUIRE(spy.size() == 1);
        const auto quads = spy.constFirst().constFirst().value<std::vector<Id>>();
        REQUIRE(quads.size() >= 5);
        // high level tiles that contain stephansdom
        // according to https://www.maptiler.com/google-maps-coordinates-tile-bounds-projection/#4/6.45/50.74
        // at the time of writing:
        // the corresponding tile list generation code is tested only by visual inspection (of the rendered scene).
        // seems to work. the following are sanity checks.
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 0, { 0, 0 } }) != quads.end());
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 1, { 1, 1 } }) != quads.end());
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 2, { 2, 2 } }) != quads.end());
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 3, { 4, 5 } }) != quads.end());
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 4, { 8, 10 } }) != quads.end());

        // currently we have tiles up until zoom level 18 for the depth tiles
        // so max quad tile level is 17.
        CHECK(std::find_if(quads.cbegin(), quads.cend(), [](const Id& id) { return id.zoom_level == 17; }) != quads.end());
        CHECK(std::find_if(quads.cbegin(), quads.cend(), [](const Id& id) { return id.zoom_level == 18; }) == quads.end());
    }

    SECTION("quads are not requested if there is no network")
    {
        auto scheduler = default_scheduler();
        scheduler->set_network_reachability(QNetworkInformation::Reachability::Disconnected);
        QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        scheduler->send_quad_requests();
        REQUIRE(spy.size() == 0);
        scheduler->set_network_reachability(QNetworkInformation::Reachability::Online);
        scheduler->send_quad_requests();
        REQUIRE(spy.size() == 1);
    }

    SECTION("delivered quads are not requested again")
    {
        auto scheduler = default_scheduler();
        scheduler->receive_quad(example_tile_quad_for(Id { 0, { 0, 0 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 1, { 1, 1 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 2, { 2, 2 } }));
        QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        scheduler->send_quad_requests();
        REQUIRE(spy.size() == 1);
        const auto quads = spy.constFirst().constFirst().value<std::vector<Id>>();
        REQUIRE(quads.size() >= 5);
        // high level tiles that contain stephansdom
        // according to https://www.maptiler.com/google-maps-coordinates-tile-bounds-projection/#4/6.45/50.74
        // at the time of writing:
        // the corresponding tile list generation code is tested only by visual inspection (of the rendered scene).
        // seems to work. the following are sanity checks.
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 0, { 0, 0 } }) == quads.end());
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 1, { 1, 1 } }) == quads.end());
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 2, { 2, 2 } }) == quads.end());
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 3, { 4, 5 } }) != quads.end());
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 4, { 8, 10 } }) != quads.end());
    }

    SECTION("network failed tiles are ignored, not found tiles are not ignored")
    {
        auto scheduler = default_scheduler();
        scheduler->receive_quad(example_tile_quad_for(Id { 0, { 0, 0 } }, 4, NetworkInfo::Status::NotFound));
        scheduler->receive_quad(example_tile_quad_for(Id { 1, { 1, 1 } }, 4, NetworkInfo::Status::NetworkError));
        scheduler->receive_quad(example_tile_quad_for(Id { 2, { 2, 2 } }, 4, NetworkInfo::Status::NotFound));
        scheduler->receive_quad(example_tile_quad_for(Id { 3, { 4, 5 } }, 4, NetworkInfo::Status::NetworkError));
        QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        scheduler->send_quad_requests();
        REQUIRE(spy.size() == 1);
        const auto quads = spy.constFirst().constFirst().value<std::vector<Id>>();
        REQUIRE(quads.size() >= 5);
        // high level tiles that contain stephansdom
        // according to https://www.maptiler.com/google-maps-coordinates-tile-bounds-projection/#4/6.45/50.74
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 0, { 0, 0 } }) == quads.end());
#ifndef __EMSCRIPTEN__
        // fails because webassembly doesn't always see 404 and therefore we need a workaround in scheduler->receive_quad, which fails the test.
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 1, { 1, 1 } }) != quads.end());
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 3, { 4, 5 } }) != quads.end());
#endif
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 2, { 2, 2 } }) == quads.end());
        CHECK(std::find(quads.cbegin(), quads.cend(), Id { 4, { 8, 10 } }) != quads.end());
    }

    SECTION("delivered tiles are requested again after they get too old")
    {
        auto scheduler = default_scheduler();
        scheduler->set_retirement_age_for_tile_cache(5 * timing_multiplicator);
        scheduler->receive_quad(example_tile_quad_for(Id { 0, { 0, 0 } }, 4, NetworkInfo::Status::Good));
        scheduler->receive_quad(example_tile_quad_for(Id { 1, { 1, 1 } }, 4, NetworkInfo::Status::NotFound));

        QSignalSpy spy(scheduler.get(), &Scheduler::quads_requested);
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        scheduler->send_quad_requests();
        {
            REQUIRE(spy.size() == 1);
            const auto quads = spy.constFirst().constFirst().value<std::vector<Id>>();
            // not found, as already in cache
            CHECK(std::find(quads.cbegin(), quads.cend(), Id { 0, { 0, 0 } }) == quads.end());
            CHECK(std::find(quads.cbegin(), quads.cend(), Id { 1, { 1, 1 } }) == quads.end());
        }

        QThread::msleep(10 * timing_multiplicator);
        scheduler->send_quad_requests();
        {
            REQUIRE(spy.size() == 2);
            const auto quads = spy.constLast().constFirst().value<std::vector<Id>>();
            // found, as cache version is to old and should be updated
            CHECK(std::find(quads.cbegin(), quads.cend(), Id { 0, { 0, 0 } }) != quads.end());
            CHECK(std::find(quads.cbegin(), quads.cend(), Id { 1, { 1, 1 } }) != quads.end());
        }
    }

    SECTION("delivered quads are sent on to the gpu (with no repeat, only the ones in the tree)")
    {
        auto scheduler = default_scheduler();
        QSignalSpy spy(scheduler.get(), &TextureScheduler::gpu_tiles_updated);
        scheduler->receive_quad(example_tile_quad_for(Id { 0, { 0, 0 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 1, { 1, 1 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 2, { 2, 2 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 4, { 8, 10 } }));
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        scheduler->update_gpu_quads();
        REQUIRE(spy.size() == 1);
        const auto gpu_tiles = spy[0][1].value<std::vector<nucleus::tile::GpuTextureTile>>();
        REQUIRE(gpu_tiles.size() == 3);
        CHECK(gpu_tiles[0].id == Id { 0, { 0, 0 } }); // order does not matter
        CHECK(gpu_tiles[1].id == Id { 1, { 1, 1 } });
        CHECK(gpu_tiles[2].id == Id { 2, { 2, 2 } });

        scheduler->receive_quad(example_tile_quad_for(Id { 3, { 4, 5 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 5, { 17, 20 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 6, { 34, 41 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 7, { 69, 83 } }));
        scheduler->update_gpu_quads();
        REQUIRE(spy.size() == 2);
        const auto new_gpu_tiles = spy[1][1].value<std::vector<nucleus::tile::GpuTextureTile>>();
        REQUIRE(new_gpu_tiles.size() == 5);
        CHECK(new_gpu_tiles[0].id == Id { 3, { 4, 5 } }); // order does not matter
        CHECK(new_gpu_tiles[1].id == Id { 4, { 8, 10 } });
        CHECK(new_gpu_tiles[2].id == Id { 5, { 17, 20 } });
        CHECK(new_gpu_tiles[3].id == Id { 6, { 34, 41 } });
        CHECK(new_gpu_tiles[4].id == Id { 7, { 69, 83 } });
    }

    SECTION("incomplete tiles are replaced with default ones, when sending to gpu")
    {
        auto scheduler = default_scheduler();
        QSignalSpy spy(scheduler.get(), &TextureScheduler::gpu_tiles_updated);
        auto quad = example_tile_quad_for({ 0, { 0, 0 } }, 4);
        quad.tiles[2].data->resize(0);

        scheduler->receive_quad(quad);
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        scheduler->update_gpu_quads();
        REQUIRE(spy.size() == 1);
        const auto gpu_tiles = spy[0][1].value<std::vector<nucleus::tile::GpuTextureTile>>();
        REQUIRE(gpu_tiles.size() == 1);
        REQUIRE(gpu_tiles[0].texture);
        CHECK(gpu_tiles[0].texture->at(0).width() == 512);
        CHECK(gpu_tiles[0].texture->at(0).height() == 512);
    }

    SECTION("gpu quads are updated when serving from cache")
    {
        auto scheduler = default_scheduler();
        QSignalSpy spy(scheduler.get(), &TextureScheduler::gpu_tiles_updated);
        for (const auto& q : example_quads_for_steffl_and_gg())
            scheduler->receive_quad(q);

        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        scheduler->update_gpu_quads();
        CHECK(spy.size() == 1);

        scheduler->update_camera(nucleus::camera::stored_positions::grossglockner());
        scheduler->update_gpu_quads();
        CHECK(spy.size() == 2);
    }

    SECTION("number of gpu quads doesn't exceede the limit")
    {
        constexpr auto tested_gpu_quad_limit = 32;
        auto scheduler = default_scheduler();
        scheduler->set_gpu_quad_limit(tested_gpu_quad_limit);
        QSignalSpy spy(scheduler.get(), &TextureScheduler::gpu_tiles_updated);

        std::unordered_set<Id, Id::Hasher> cached_tiles;

        {
            scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
            for (const auto& q : example_quads_for(scheduler->missing_quads_for_current_camera()))
                scheduler->receive_quad(q);

            scheduler->update_gpu_quads();
            REQUIRE(spy.size() == 1);
            const auto new_quads = spy[0][1].value<std::vector<nucleus::tile::GpuTextureTile>>();
            const auto deleted_quads = spy[0][0].value<std::vector<Id>>();
            CHECK(new_quads.size() == tested_gpu_quad_limit);
            CHECK(deleted_quads.empty());

            nucleus::tile::Cache<nucleus::tile::GpuTextureTile> test_cache;
            for (const auto& q : new_quads)
                test_cache.insert(q);

            auto n_tiles = 0;
            test_cache.visit([&n_tiles](const auto&) {
                n_tiles++;
                return true;
            });
            CHECK(n_tiles == tested_gpu_quad_limit); // check that all are reachable
            for (const auto& tile : new_quads) {
                cached_tiles.insert(tile.id);
            }
        }

        {
            scheduler->update_camera(nucleus::camera::stored_positions::grossglockner());
            for (const auto& q : example_quads_for(scheduler->missing_quads_for_current_camera()))
                scheduler->receive_quad(q);

            scheduler->update_gpu_quads();
            REQUIRE(spy.size() == 2);
            const auto new_quads = spy[1][1].value<std::vector<nucleus::tile::GpuTextureTile>>();
            const auto deleted_quads = spy[1][0].value<std::vector<Id>>();
            CHECK(new_quads.size() == deleted_quads.size());
            CHECK(!new_quads.empty());

            { // check that all are reachable
                for (const auto& id : deleted_quads) {
                    CHECK(cached_tiles.contains(id));
                    cached_tiles.erase(id);
                }
                for (const auto& tile : new_quads) {
                    CHECK(!cached_tiles.contains(tile.id));
                    cached_tiles.insert(tile.id);
                }
                CHECK(cached_tiles.size() == tested_gpu_quad_limit);
                nucleus::tile::Cache<nucleus::tile::GpuCacheInfo> test_cache;
                for (const auto& id : cached_tiles) {
                    test_cache.insert({ { id } });
                }
                auto n_tiles = 0;
                test_cache.visit([&n_tiles](const auto&) {
                    n_tiles++;
                    return true;
                });
                CHECK(n_tiles == tested_gpu_quad_limit);
            }
            { // check for double entries
                std::unordered_set<Id, Id::Hasher> deleted_quads_set(deleted_quads.cbegin(), deleted_quads.cend());
                for (const auto& new_quad : new_quads) {
                    CHECK(!deleted_quads_set.contains(new_quad.id));
                }
            }
        }
        {
            scheduler->update_camera(nucleus::camera::stored_positions::grossglockner_shadow());
            for (const auto& q : example_quads_for(scheduler->missing_quads_for_current_camera()))
                scheduler->receive_quad(q);

            scheduler->update_gpu_quads();
            REQUIRE(spy.size() == 3);
            const auto new_quads = spy[2][1].value<std::vector<nucleus::tile::GpuTextureTile>>();
            const auto deleted_quads = spy[2][0].value<std::vector<Id>>();
            CHECK(new_quads.size() == deleted_quads.size());
            CHECK(!new_quads.empty());

            { // check that all are reachable
                for (const auto& id : deleted_quads) {
                    CHECK(cached_tiles.contains(id));
                    cached_tiles.erase(id);
                }
                for (const auto& tile : new_quads) {
                    CHECK(!cached_tiles.contains(tile.id));
                    cached_tiles.insert(tile.id);
                }
                CHECK(cached_tiles.size() == tested_gpu_quad_limit);
                nucleus::tile::Cache<nucleus::tile::GpuCacheInfo> test_cache;
                for (const auto& id : cached_tiles) {
                    test_cache.insert({ { id } });
                }
                auto n_tiles = 0;
                test_cache.visit([&n_tiles](const auto&) {
                    n_tiles++;
                    return true;
                });
                CHECK(n_tiles == tested_gpu_quad_limit);
            }
            { // check for double entries
                std::unordered_set<Id, Id::Hasher> deleted_quads_set(deleted_quads.cbegin(), deleted_quads.cend());
                for (const auto& new_quad : new_quads) {
                    CHECK(!deleted_quads_set.contains(new_quad.id));
                }
            }
        }
    }

    SECTION("gpu tiles are optimised for the current camera position")
    {
        auto scheduler = default_scheduler();
        scheduler->set_gpu_quad_limit(17);
        QSignalSpy spy(scheduler.get(), &TextureScheduler::gpu_tiles_updated);
        for (const auto& q : example_quads_for_steffl_and_gg())
            scheduler->receive_quad(q);

        std::unordered_set<Id, Id::Hasher> cached_tiles;
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        scheduler->update_gpu_quads();
        REQUIRE(spy.size() == 1);
        const auto new_quads = spy[0][1].value<std::vector<nucleus::tile::GpuTextureTile>>();
        for (const auto& tile : new_quads) {
            cached_tiles.insert(tile.id);
        }
        CHECK(cached_tiles.contains({ 11, { 1117, 1337 } }));
        CHECK(cached_tiles.contains({ 11, { 1117, 1338 } }));
        CHECK(cached_tiles.contains({ 11, { 1116, 1337 } }));
        CHECK(cached_tiles.contains({ 11, { 1116, 1338 } }));
        CHECK(cached_tiles.contains({ 12, { 2234, 2675 } }));
    }

    SECTION("ram tiles are purged")
    {
        auto scheduler = default_scheduler();
        scheduler->set_ram_quad_limit(17);
        scheduler->set_purge_timeout(1);
        for (const auto& q : example_quads_for_steffl_and_gg())
            scheduler->receive_quad(q);
        test_helpers::process_events_for(2 * timing_multiplicator);
        CHECK(scheduler->ram_cache().n_cached_objects() == 17);
    }

    SECTION("purging tiles based on camera")
    {
        auto scheduler = default_scheduler();
        scheduler->set_ram_quad_limit(17);
        scheduler->set_purge_timeout(2 * timing_multiplicator);
        scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
        for (const auto& q : example_quads_for_steffl_and_gg())
            scheduler->receive_quad(q);
        test_helpers::process_events_for(3 * timing_multiplicator);
        CHECK(scheduler->ram_cache().n_cached_objects() == 17);
        CHECK(scheduler->ram_cache().contains({ 11, { 1117, 1337 } }));
        CHECK(scheduler->ram_cache().contains({ 11, { 1117, 1338 } }));
        CHECK(scheduler->ram_cache().contains({ 11, { 1116, 1337 } }));
        CHECK(scheduler->ram_cache().contains({ 11, { 1116, 1338 } }));
        CHECK(scheduler->ram_cache().contains({ 12, { 2234, 2675 } }));
    }

    SECTION("purging ram tiles with tolerance")
    {
        auto scheduler = default_scheduler();
        // example_quads_for_steffl_and_gg().size() == 39
        const unsigned limit = 38;
        REQUIRE(example_quads_for_steffl_and_gg().size() > limit);
        REQUIRE(example_quads_for_steffl_and_gg().size() == 39);
        scheduler->set_ram_quad_limit(limit);
        for (const auto& q : example_quads_for_steffl_and_gg())
            scheduler->receive_quad(q);
        scheduler->purge_ram_cache();
        CHECK(scheduler->ram_cache().n_cached_objects() == example_quads_for_steffl_and_gg().size());
        scheduler->receive_quad(example_tile_quad_for(Id { 10, { 0, 0 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 11, { 1, 1 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 12, { 2, 2 } }));
        scheduler->purge_ram_cache();
        CHECK(scheduler->ram_cache().n_cached_objects() == limit);
    }

#ifndef __EMSCRIPTEN__
    SECTION("purging happens with a delay (collects purge events) and the timer is not restarted on tile delivery")
    {
        nucleus::utils::Stopwatch sw;
        auto scheduler = default_scheduler();
        scheduler->set_purge_timeout(9 * timing_multiplicator);
        scheduler->set_update_timeout(100 * timing_multiplicator); // sending quads to gpu takes long and makes tight timing impossible in debug mode
        scheduler->set_persist_timeout(100 * timing_multiplicator);
        scheduler->set_ram_quad_limit(2);
        scheduler->receive_quad(example_tile_quad_for(Id { 0, { 0, 0 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 1, { 1, 1 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 2, { 2, 2 } }));
        CHECK(scheduler->ram_cache().n_cached_objects() == 3);
        test_helpers::process_events_for(3 * timing_multiplicator);
        CHECK(scheduler->ram_cache().n_cached_objects() == 3);
        scheduler->receive_quad(example_tile_quad_for(Id { 1, { 0, 0 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 1, { 1, 0 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 2, { 2, 1 } }));

        CHECK(scheduler->ram_cache().n_cached_objects() == 6);
        test_helpers::process_events_for(3 * timing_multiplicator);
        CHECK(scheduler->ram_cache().n_cached_objects() == 6);

        scheduler->receive_quad(example_tile_quad_for(Id { 1, { 0, 1 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 2, { 1, 1 } }));
        scheduler->receive_quad(example_tile_quad_for(Id { 2, { 1, 2 } }));
        CHECK(scheduler->ram_cache().n_cached_objects() == 9);
        test_helpers::process_events_for(4 * timing_multiplicator);
        CHECK(scheduler->ram_cache().n_cached_objects() == 2);
    }
#endif

    const auto check_persisted_tile = [](const auto& scheduler, const Id& id) {
        const auto example_quad = example_tile_quad_for(id);
        REQUIRE(scheduler->ram_cache().contains(id));
        REQUIRE(scheduler->ram_cache().peak_at(id).n_tiles == example_quad.n_tiles);
        REQUIRE(scheduler->ram_cache().peak_at(id).id == id);
        REQUIRE(id == example_quad.id);
        const auto children = id.children();
        for (unsigned i = 0; i < 4; ++i) {
            const nucleus::tile::Data& child_tile = scheduler->ram_cache().peak_at(id).tiles[i];
            CHECK(child_tile.id == children[i]);
            CHECK(*child_tile.data == *example_quad.tiles[i].data);
        }
    };

    const auto check_persited_tiles = [&check_persisted_tile](const auto& scheduler, const auto& ids) {
        for (const auto& id : ids) {
            check_persisted_tile(scheduler, id);
        }
    };

    SECTION("persisting data does error on unnamed schedulers")
    {

        auto scheduler = std::make_unique<TextureScheduler>(Scheduler::Settings {});
        CHECK(!scheduler->persist_tiles());
        CHECK(!scheduler->read_disk_cache());
    }

    SECTION("persisting data works")
    {
        {
            auto scheduler = default_scheduler();
            scheduler->receive_quad(example_tile_quad_for(Id { 0, { 0, 0 } }));
            scheduler->receive_quad(example_tile_quad_for(Id { 1, { 1, 1 } }));
            scheduler->receive_quad(example_tile_quad_for(Id { 2, { 2, 2 } }));
            CHECK(scheduler->persist_tiles());
        }
        auto scheduler = scheduler_with_disk_cache();
        CHECK(scheduler->ram_cache().n_cached_objects() == 3);

        check_persited_tiles(scheduler, std::vector { Id { 0, { 0, 0 } }, Id { 1, { 1, 1 } }, Id { 2, { 2, 2 } } });

        CHECK(scheduler->ram_cache().peak_at(Id { 0, { 0, 0 } }).tiles[0].id == Id { 1, { 0, 0 } }); // order does not matter!
        CHECK(scheduler->ram_cache().peak_at(Id { 0, { 0, 0 } }).tiles[1].id == Id { 1, { 1, 0 } });
        CHECK(scheduler->ram_cache().peak_at(Id { 0, { 0, 0 } }).tiles[2].id == Id { 1, { 0, 1 } });
        CHECK(scheduler->ram_cache().peak_at(Id { 0, { 0, 0 } }).tiles[3].id == Id { 1, { 1, 1 } });

        CHECK(scheduler->ram_cache().peak_at(Id { 1, { 1, 1 } }).tiles[0].id == Id { 2, { 2, 2 } });
        CHECK(scheduler->ram_cache().peak_at(Id { 1, { 1, 1 } }).tiles[1].id == Id { 2, { 3, 2 } });
        CHECK(scheduler->ram_cache().peak_at(Id { 1, { 1, 1 } }).tiles[2].id == Id { 2, { 2, 3 } });
        CHECK(scheduler->ram_cache().peak_at(Id { 1, { 1, 1 } }).tiles[3].id == Id { 2, { 3, 3 } });

        CHECK(scheduler->ram_cache().peak_at(Id { 2, { 2, 2 } }).tiles[0].id == Id { 3, { 4, 4 } });
        CHECK(scheduler->ram_cache().peak_at(Id { 2, { 2, 2 } }).tiles[1].id == Id { 3, { 5, 4 } });
        CHECK(scheduler->ram_cache().peak_at(Id { 2, { 2, 2 } }).tiles[2].id == Id { 3, { 4, 5 } });
        CHECK(scheduler->ram_cache().peak_at(Id { 2, { 2, 2 } }).tiles[3].id == Id { 3, { 5, 5 } });
        std::filesystem::remove_all(scheduler->disk_cache_path());
    }

    SECTION("persisting data works als with itterative updates")
    {
        {
            auto scheduler = default_scheduler();
            scheduler->receive_quad(example_tile_quad_for(Id { 0, { 0, 0 } }));
            scheduler->receive_quad(example_tile_quad_for(Id { 1, { 1, 1 } }));
            scheduler->receive_quad(example_tile_quad_for(Id { 2, { 2, 2 } }));
            scheduler->persist_tiles();
        }
        {
            auto scheduler = scheduler_with_disk_cache();
            CHECK(scheduler->ram_cache().n_cached_objects() == 3);
            check_persited_tiles(scheduler, std::vector { Id { 0, { 0, 0 } }, Id { 1, { 1, 1 } }, Id { 2, { 2, 2 } } });
            scheduler->receive_quad(example_tile_quad_for(Id { 3, { 0, 0 } }));
            scheduler->receive_quad(example_tile_quad_for(Id { 4, { 0, 1 } }));
            scheduler->persist_tiles();
        }
        {
            auto scheduler = scheduler_with_disk_cache();
            CHECK(scheduler->ram_cache().n_cached_objects() == 5);
            check_persited_tiles(scheduler, std::vector { Id { 0, { 0, 0 } }, Id { 1, { 1, 1 } }, Id { 2, { 2, 2 } }, Id { 3, { 0, 0 } }, Id { 4, { 0, 1 } } });
            scheduler->set_ram_quad_limit(3);
            scheduler->update_camera(nucleus::camera::stored_positions::stephansdom());
            scheduler->update_gpu_quads();
            scheduler->purge_ram_cache();
            CHECK(scheduler->ram_cache().n_cached_objects() == 3);
            check_persited_tiles(scheduler, std::vector { Id { 0, { 0, 0 } }, Id { 1, { 1, 1 } }, Id { 2, { 2, 2 } } });
            scheduler->persist_tiles();
        }
        {
            auto scheduler = scheduler_with_disk_cache();
            CHECK(scheduler->ram_cache().n_cached_objects() == 3);
            check_persited_tiles(scheduler, std::vector { Id { 0, { 0, 0 } }, Id { 1, { 1, 1 } }, Id { 2, { 2, 2 } } });
            std::filesystem::remove_all(scheduler->disk_cache_path());
        }
    }

    SECTION("notification, when a tile is received")
    {
        auto scheduler = default_scheduler();
        QSignalSpy spy(scheduler.get(), &Scheduler::quad_received);
        scheduler->receive_quad(example_tile_quad_for(Id { 0, { 0, 0 } }));
        REQUIRE(spy.size() == 1);
        scheduler->receive_quad(example_tile_quad_for(Id { 1, { 1, 1 } }));
        REQUIRE(spy.size() == 2);
        CHECK(spy[0][0].value<Id>() == Id { 0, { 0, 0 } });
        CHECK(spy[1][0].value<Id>() == Id { 1, { 1, 1 } });
    }
}

TEST_CASE("nucleus/tile/Scheduler benchmarks")
{
    auto camera = nucleus::camera::stored_positions::grossglockner();
    camera.set_viewport_size({ 3840, 2160 });

    BENCHMARK("construct")
    {
        return scheduler_with_true_heights();
    };

    {
        auto scheduler = scheduler_with_true_heights();
        BENCHMARK("request quads")
        {
            scheduler->update_camera(camera);
            scheduler->send_quad_requests();
        };
    }

    BENCHMARK("receive " + std::to_string(example_quads_many().size()) + " quads")
    {
        auto scheduler = default_scheduler();
        for (const auto& q : example_quads_many())
            scheduler->receive_quad(q);
    };

    BENCHMARK("receive " + std::to_string(example_quads_for_steffl_and_gg().size()) + " quads + update_gpu_quads")
    {
        auto scheduler = default_scheduler();
        scheduler->update_camera(camera);
        for (const auto& q : example_quads_for_steffl_and_gg())
            scheduler->receive_quad(q);
        // unpacking byte arrays takes long, hence only the smaller dataset
        scheduler->update_gpu_quads();
    };

    BENCHMARK("receive " + std::to_string(example_quads_many().size()) + " quads + purge_ram_cache")
    {
        auto scheduler = default_scheduler();
        scheduler->update_camera(camera);
        for (const auto& q : example_quads_for_steffl_and_gg())
            scheduler->receive_quad(q);
        scheduler->purge_ram_cache();
    };

    {
        auto scheduler = default_scheduler();
        scheduler->receive_quad({example_tile_quad_for({0, {0, 0}}),});
        for (unsigned i = 1; i < 100; ++i) {
            scheduler->receive_quad(example_tile_quad_for({ i, { 0, 0 } }));
            scheduler->receive_quad(example_tile_quad_for({ i, { 1, 0 } }));
            scheduler->receive_quad(example_tile_quad_for({ i, { 1, 1 } }));
            scheduler->receive_quad(example_tile_quad_for({ i, { 0, 1 } }));
        }
        BENCHMARK("write cache to disk") {
            scheduler->persist_tiles();
        };
    }

    BENCHMARK("read cache from disk") {
        auto scheduler = scheduler_with_disk_cache();
    };
    auto scheduler = scheduler_with_disk_cache();
    std::filesystem::remove_all(scheduler->disk_cache_path());
}

TEST_CASE("nucleus/tile/TextureScheduler")
{
    SECTION("to_raster")
    {
        nucleus::tile::DataQuad quad;
        quad.id = radix::tile::Id { 6, { 34, 41 } };
        unsigned idx = 1;
        for (const auto& c : quad.id.children()) {
            quad.tiles[idx].id = c;
            auto ba = test_helpers::load_test_file(QString("quad/%1_%2_%3.jpg").arg(c.zoom_level).arg(c.coords.x).arg(c.coords.y));
            quad.tiles[idx].data = std::make_shared<QByteArray>(std::move(ba));
            quad.tiles[idx].network_info = { NetworkInfo::Status::Good, 12345 };
            idx = (idx + 1) % 4;
        }
        quad.n_tiles = 4;
        const auto joined = TextureScheduler::to_raster(quad, { { 256, 256 }, glm::u8vec4 { 255, 255, 255, 255 } });
        const auto correct = nucleus::utils::image_loader::rgba8(test_helpers::load_test_file(QString("quad/merged.jpg")));
        REQUIRE(joined.width() == correct->width());
        REQUIRE(joined.height() == correct->height());
        auto correct_iter = correct->begin();
        auto rmse = 0.0;
        for (const auto& v : joined) {
            const auto e = glm::vec4(v) - glm::vec4(*(correct_iter++));
            rmse += glm::dot(e, e);
        }
        rmse = std::sqrt(rmse / (correct->width() * correct->height() * 3));
        CHECK(rmse < 1);
        const auto qimage = nucleus::tile::conversion::to_QImage(joined);
        qimage.save("merged.png");
    }
}

TEST_CASE("nucleus/tile/SchedulerDirector")
{
    SECTION("api")
    {
        std::shared_ptr<Scheduler> sch1 = default_scheduler();
        std::shared_ptr<Scheduler> sch2 = default_scheduler();
        SchedulerDirector d;
        d.check_in("sch1", sch1);
        d.check_in("sch2", sch2);
        CHECK(sch1->name() == "sch1");
        CHECK(sch2->name() == "sch2");
        REQUIRE(sch1->enabled() == true);
        REQUIRE(sch2->enabled() == true);
        d.visit([](Scheduler* sch) { sch->set_enabled(false); });
        CHECK(sch1->enabled() == false);
        CHECK(sch2->enabled() == false);
    }
    SECTION("no two entries with the same name")
    {
        SchedulerDirector reg;
        CHECK(reg.check_in("name", default_scheduler()));
        CHECK(!reg.check_in("name", default_scheduler()));
    }
}
