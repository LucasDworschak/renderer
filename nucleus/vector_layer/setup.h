/*****************************************************************************
 * AlpineMaps.org
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

#pragma once

#include "Scheduler.h"
#include <QThread>
#include <memory>
#include <nucleus/tile/QuadAssembler.h>
#include <nucleus/tile/RateLimiter.h>
#include <nucleus/tile/SlotLimiter.h>
#include <nucleus/tile/TileLoadService.h>

#include <QTimer>

#include <QDebug>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QtVersionChecks>

namespace nucleus::vector_layer::setup {

class TileLoadServiceMock : public nucleus::tile::TileLoadService {
public:
    TileLoadServiceMock(const QString& base_url, UrlPattern url_pattern, const QString& file_ending, const LoadBalancingTargets& load_balancing_targets = {})
        : nucleus::tile::TileLoadService(base_url, url_pattern, file_ending, load_balancing_targets)
    {
    }
public slots:
    void load_mock(const tile::Id& tile_id) const
    {

        QTimer* t = new QTimer();
        t->setSingleShot(true);

        // do some fancy stuff after timeout
        connect(t, &QTimer::timeout, [tile_id, t, this]() {
            const auto timestamp = utils::time_since_epoch();

            auto tile = std::make_shared<QByteArray>();
            emit load_finished({ tile_id, { tile::NetworkInfo::Status::Good, timestamp }, tile });
            t->deleteLater();
        });

        t->start(100);
    }
};
using TileLoadServicePtr = std::unique_ptr<nucleus::vector_layer::setup::TileLoadServiceMock>;

// using TileLoadServicePtr = std::unique_ptr<nucleus::tile::TileLoadService>;

struct SchedulerHolder {
    std::unique_ptr<vector_layer::Scheduler> scheduler;
    TileLoadServicePtr tile_service;
};

SchedulerHolder scheduler(std::string name, TileLoadServicePtr tile_service, const tile::utils::AabbDecoratorPtr& aabb_decorator, QThread* thread = nullptr)
{
    auto scheduler = std::make_unique<nucleus::vector_layer::Scheduler>(std::move(name));
    scheduler->read_disk_cache();
    scheduler->set_gpu_quad_limit(512);
    scheduler->set_ram_quad_limit(12000);
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
        QObject::connect(qa, &QuadAssembler::tile_requested, tile_service.get(), &TileLoadServiceMock::load_mock);
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