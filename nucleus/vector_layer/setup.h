/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2024 Adam Celarek
 * Copyright (C) 2024 Lucas Dworschak
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

#pragma once

#include "Scheduler.h"
#include <QThread>
#include <memory>
#include <nucleus/tile/QuadAssembler.h>
#include <nucleus/tile/RateLimiter.h>
#include <nucleus/tile/SlotLimiter.h>
#include <nucleus/tile/TileLoadService.h>
#include <nucleus/vector_layer/constants.h>

#include <QTimer>

#include <QDebug>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QtVersionChecks>

namespace nucleus::vector_layer::setup {

using TileLoadServicePtr = std::unique_ptr<nucleus::tile::TileLoadService>;

struct SchedulerHolder {
    std::shared_ptr<vector_layer::Scheduler> scheduler;
    TileLoadServicePtr tile_service;
};

SchedulerHolder scheduler(TileLoadServicePtr tile_service, const tile::utils::AabbDecoratorPtr& aabb_decorator, QThread* thread = nullptr)
{

    // Style style(":/vectorlayerstyles/basemap.json");
    Style style(":/vectorlayerstyles/openstreetmap.json");
    // Style style(":/vectorlayerstyles/qwant.json");
    // Style style(":/vectorlayerstyles/osm-bright.json");

    style.load();

    Scheduler::Settings settings;
    settings.max_zoom_level = constants::style_zoom_range.y;
    settings.tile_resolution = 256;
    settings.gpu_quad_limit = 512;
    auto scheduler = std::make_unique<nucleus::vector_layer::Scheduler>(settings, std::move(style));
    scheduler->set_aabb_decorator(aabb_decorator);

    {
        using nucleus::tile::QuadAssembler;
        using nucleus::tile::RateLimiter;
        using nucleus::tile::SlotLimiter;
        using nucleus::tile::TileLoadService;
        auto* sch = scheduler.get();
        auto* sl = new SlotLimiter(sch);
        auto* rl = new RateLimiter(sch);
        auto* qa = new QuadAssembler(sch);

        QObject::connect(sch, &Scheduler::quads_requested, sl, &SlotLimiter::request_quads);
        QObject::connect(sl, &SlotLimiter::quad_requested, rl, &RateLimiter::request_quad);
        QObject::connect(rl, &RateLimiter::quad_requested, qa, &QuadAssembler::load);
        QObject::connect(qa, &QuadAssembler::tile_requested, tile_service.get(), &TileLoadService::load);
        QObject::connect(tile_service.get(), &TileLoadService::load_finished, qa, &QuadAssembler::deliver_tile);

        QObject::connect(qa, &QuadAssembler::quad_loaded, sl, &SlotLimiter::deliver_quad);
        QObject::connect(sl, &SlotLimiter::quad_delivered, sch, &nucleus::vector_layer::Scheduler::receive_quad);
    }
    if (QNetworkInformation::loadDefaultBackend() && QNetworkInformation::instance()) {
        QNetworkInformation* n = QNetworkInformation::instance();
        scheduler->set_network_reachability(n->reachability());
        QObject::connect(n, &QNetworkInformation::reachabilityChanged, scheduler.get(), &Scheduler::set_network_reachability);
    }

    Q_UNUSED(thread);
#ifdef ALP_ENABLE_THREADING
#ifdef __EMSCRIPTEN__ // make request from main thread on webassembly due to QTBUG-109396
    tile_service->moveToThread(QCoreApplication::instance()->thread());
#else
    if (thread)
        tile_service->moveToThread(thread);
#endif
    if (thread)
        scheduler->moveToThread(thread);
#endif

    return { std::move(scheduler), std::move(tile_service) };
}
} // namespace nucleus::vector_layer::setup
