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

#include "RenderingContext.h"

#include "RenderThreadNotifier.h"
#include <QCoreApplication>
#include <QMutex>
#include <QOpenGLContext>
#include <QThread>
#include <gl_engine/Context.h>
#include <gl_engine/MapLabels.h>
#include <gl_engine/TextureLayer.h>
#include <gl_engine/TileGeometry.h>
#include <nucleus/DataQuerier.h>
#include <nucleus/camera/Controller.h>
#include <nucleus/camera/PositionStorage.h>
#include <nucleus/map_label/Filter.h>
#include <nucleus/map_label/setup.h>
#include <nucleus/picker/PickerManager.h>
#include <nucleus/tile_scheduler/OldScheduler.h>
#include <nucleus/tile_scheduler/TextureScheduler.h>
#include <nucleus/tile_scheduler/TileLoadService.h>
#include <nucleus/tile_scheduler/setup.h>
#include <nucleus/utils/thread.h>

using namespace nucleus::tile_scheduler;
using namespace nucleus::map_label;
using namespace nucleus::picker;
using nucleus::DataQuerier;

struct RenderingContext::Data {
    // WARNING: gl_engine::Context must be on the rendering thread!!
    mutable QMutex shared_ptr_mutex; // protects the shared_ptr
    std::shared_ptr<gl_engine::Context> engine_context;

    std::unique_ptr<QThread> scheduler_thread;
    // the ones below are on the scheduler thread.
    nucleus::tile_scheduler::setup::MonolithicScheduler scheduler;
    nucleus::tile_scheduler::setup::TextureSchedulerHolder ortho_texture;
    nucleus::map_label::setup::SchedulerHolder map_label;
    std::shared_ptr<nucleus::DataQuerier> data_querier;
    std::unique_ptr<nucleus::camera::Controller> camera_controller;
    std::shared_ptr<nucleus::map_label::Filter> label_filter;
    std::shared_ptr<nucleus::picker::PickerManager> picker_manager;
    std::shared_ptr<nucleus::tile_scheduler::utils::AabbDecorator> aabb_decorator;
};

RenderingContext::RenderingContext(QObject* parent)
    : QObject { parent }
    , m(std::make_unique<RenderingContext::Data>())
{
    qDebug("RenderingContext::RenderingContext()");
    using TilePattern = nucleus::tile_scheduler::TileLoadService::UrlPattern;
    assert(QThread::currentThread() == QCoreApplication::instance()->thread());

    auto terrain_service = std::make_unique<TileLoadService>("https://alpinemaps.cg.tuwien.ac.at/tiles/alpine_png/", TilePattern::ZXY, ".png");
    //    m->ortho_service.reset(new TileLoadService("https://tiles.bergfex.at/styles/bergfex-osm/", TileLoadService::UrlPattern::ZXY_yPointingSouth,
    //    ".jpeg")); m->ortho_service.reset(new TileLoadService("https://alpinemaps.cg.tuwien.ac.at/tiles/ortho/",
    //    TileLoadService::UrlPattern::ZYX_yPointingSouth, ".jpeg"));
    // m->ortho_service.reset(new TileLoadService("https://maps%1.wien.gv.at/basemap/bmaporthofoto30cm/normal/google3857/",
    //                                           TileLoadService::UrlPattern::ZYX_yPointingSouth,
    //                                           ".jpeg",
    //                                           {"", "1", "2", "3", "4"}));
    auto ortho_service
        = std::make_unique<TileLoadService>("https://gataki.cg.tuwien.ac.at/raw/basemap/tiles/", TileLoadService::UrlPattern::ZYX_yPointingSouth, ".jpeg");
    m->aabb_decorator = nucleus::tile_scheduler::setup::aabb_decorator();
    m->scheduler = nucleus::tile_scheduler::setup::monolithic(std::move(terrain_service), std::move(ortho_service), m->aabb_decorator);
    m->data_querier = std::make_shared<DataQuerier>(&m->scheduler.scheduler->ram_cache());

    m->ortho_texture = nucleus::tile_scheduler::setup::texture_scheduler(
        std::make_unique<TileLoadService>("https://gataki.cg.tuwien.ac.at/raw/basemap/tiles/", TileLoadService::UrlPattern::ZYX_yPointingSouth, ".jpeg"), m->aabb_decorator, m->scheduler.thread.get());

    m->map_label = nucleus::map_label::setup::scheduler(
        std::make_unique<TileLoadService>("https://osm.cg.tuwien.ac.at/vector_tiles/poi_v1/", TilePattern::ZXY_yPointingSouth, ""), m->aabb_decorator, m->data_querier, m->scheduler.thread.get());
    m->map_label.scheduler->set_geometry_ram_cache(&m->scheduler.scheduler->ram_cache());
    m->scheduler.scheduler->set_dataquerier(m->data_querier);

    m->picker_manager = std::make_shared<PickerManager>();
    m->label_filter = std::make_shared<Filter>();
    if (m->scheduler.thread) {
        m->picker_manager->moveToThread(m->scheduler.thread.get());
        m->label_filter->moveToThread(m->scheduler.thread.get());
    }
    // clang-format off
    connect(m->scheduler.scheduler.get(),     &nucleus::tile_scheduler::OldScheduler::gpu_quads_updated,     RenderThreadNotifier::instance(), &RenderThreadNotifier::notify);
    connect(m->ortho_texture.scheduler.get(), &nucleus::tile_scheduler::TextureScheduler::gpu_quads_updated, RenderThreadNotifier::instance(), &RenderThreadNotifier::notify);
    connect(m->map_label.scheduler.get(),     &nucleus::map_label::Scheduler::gpu_tiles_updated,             RenderThreadNotifier::instance(), &RenderThreadNotifier::notify);
    connect(m->map_label.scheduler.get(),     &nucleus::map_label::Scheduler::gpu_tiles_updated,             m->picker_manager.get(),          &PickerManager::update_quads);
    connect(m->map_label.scheduler.get(),     &nucleus::map_label::Scheduler::gpu_tiles_updated,             m->label_filter.get(),            &Filter::update_quads);
    // clang-format on

    if (QNetworkInformation::loadDefaultBackend() && QNetworkInformation::instance()) {
        QNetworkInformation* n = QNetworkInformation::instance();
        m->scheduler.scheduler->set_network_reachability(n->reachability());
        // clang-format off
        connect(n, &QNetworkInformation::reachabilityChanged, m->scheduler.scheduler.get(),                                  &OldScheduler::set_network_reachability);
        connect(n, &QNetworkInformation::reachabilityChanged, m->ortho_texture.scheduler.get(), &nucleus::tile_scheduler::TextureScheduler::set_network_reachability);
        connect(n, &QNetworkInformation::reachabilityChanged, m->map_label.scheduler.get(),                 &nucleus::map_label::Scheduler::set_network_reachability);
        // clang-format on
    }
}

RenderingContext::~RenderingContext()
{
    // being killed after app exited -> no event loop
}

RenderingContext* RenderingContext::instance()
{

    static RenderingContext s_instance;
    return &s_instance;
}

void RenderingContext::initialise()
{
    QMutexLocker locker(&m->shared_ptr_mutex);
    if (m->engine_context)
        return;

    m->engine_context = std::make_shared<gl_engine::Context>();
    // standard tiles
    m->engine_context->set_tile_geometry(std::make_shared<gl_engine::TileGeometry>());
    m->engine_context->set_ortho_layer(std::make_shared<gl_engine::TextureLayer>());
    m->engine_context->tile_geometry()->set_quad_limit(512);
    m->engine_context->tile_geometry()->set_aabb_decorator(m->aabb_decorator);
    m->engine_context->ortho_layer()->set_quad_limit(512);

    connect(m->scheduler.scheduler.get(), &nucleus::tile_scheduler::OldScheduler::gpu_quads_updated, m->engine_context->tile_geometry(), &gl_engine::TileGeometry::update_gpu_quads);
    connect(m->ortho_texture.scheduler.get(), &nucleus::tile_scheduler::TextureScheduler::gpu_quads_updated, m->engine_context->ortho_layer(), &gl_engine::TextureLayer::update_gpu_quads);
    nucleus::utils::thread::async_call(m->scheduler.scheduler.get(), [this]() { m->scheduler.scheduler->set_enabled(true); });
    const auto texture_compression = gl_engine::Texture::compression_algorithm();
    nucleus::utils::thread::async_call(m->ortho_texture.scheduler.get(), [this, texture_compression]() {
        m->ortho_texture.scheduler->set_texture_compression_algorithm(texture_compression);
        // m->ortho_texture.scheduler->set_enabled(true);
    });

    // labels
    m->engine_context->set_map_label_manager(std::make_unique<gl_engine::MapLabels>(m->aabb_decorator));
    connect(m->label_filter.get(), &Filter::filter_finished, m->engine_context->map_label_manager(), &gl_engine::MapLabels::update_labels);
    nucleus::utils::thread::async_call(m->map_label.scheduler.get(), [this]() { m->map_label.scheduler->set_enabled(true); });

    connect(QOpenGLContext::currentContext(), &QOpenGLContext::aboutToBeDestroyed, m->engine_context.get(), &nucleus::EngineContext::destroy);
    connect(QOpenGLContext::currentContext(), &QOpenGLContext::aboutToBeDestroyed, this, &RenderingContext::destroy);

    m->engine_context->initialise();
    nucleus::utils::thread::async_call(this, [this]() { emit this->initialised(); });
}

void RenderingContext::destroy()
{
    nucleus::utils::thread::sync_call(m->scheduler.scheduler.get(), [this]() {
        m->scheduler.scheduler.reset();
        m->camera_controller.reset();
        m->label_filter.reset();
        m->picker_manager.reset();
        m->map_label.scheduler.reset();
        m->ortho_texture.scheduler.reset();
    });
    if (m->scheduler.thread) {
        m->scheduler.thread->quit();
        m->scheduler.thread->wait(500); // msec
        m->scheduler.thread.reset();
    }
}

std::shared_ptr<gl_engine::Context> RenderingContext::engine_context() const
{
    QMutexLocker locker(&m->shared_ptr_mutex);
    return m->engine_context;
}

std::shared_ptr<nucleus::tile_scheduler::utils::AabbDecorator> RenderingContext::aabb_decorator() const
{
    QMutexLocker locker(&m->shared_ptr_mutex);
    return m->aabb_decorator;
}

std::shared_ptr<nucleus::DataQuerier> RenderingContext::data_querier() const
{
    QMutexLocker locker(&m->shared_ptr_mutex);
    return m->data_querier;
}

OldScheduler* RenderingContext::scheduler() const
{
    QMutexLocker locker(&m->shared_ptr_mutex);
    return m->scheduler.scheduler.get();
}

std::shared_ptr<nucleus::picker::PickerManager> RenderingContext::picker_manager() const
{
    QMutexLocker locker(&m->shared_ptr_mutex);
    return m->picker_manager;
}

std::shared_ptr<nucleus::map_label::Filter> RenderingContext::label_filter() const
{
    QMutexLocker locker(&m->shared_ptr_mutex);
    return m->label_filter;
}

nucleus::map_label::Scheduler* RenderingContext::map_label_scheduler() const
{
    QMutexLocker locker(&m->shared_ptr_mutex);
    return m->map_label.scheduler.get();
}
