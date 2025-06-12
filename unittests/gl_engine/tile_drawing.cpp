/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2025 Lucas Dworschak
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
#include "UnittestGLContext.h"

#include <vector>

#include <QImage>
#include <QSignalSpy>
#include <catch2/catch_test_macros.hpp>

#include <QNetworkInformation>

#include <gl_engine/Framebuffer.h>
#include <gl_engine/ShaderProgram.h>
#include <gl_engine/ShaderRegistry.h>
#include <gl_engine/Texture.h>
#include <gl_engine/TileGeometry.h>
#include <gl_engine/UniformBuffer.h>
#include <gl_engine/UniformBufferObjects.h>
#include <gl_engine/VectorLayer.h>
#include <gl_engine/helpers.h>
#include <nucleus/camera/PositionStorage.h>
#include <nucleus/srs.h>
#include <nucleus/tile/drawing.h>
#include <nucleus/tile/setup.h>
#include <nucleus/tile/utils.h>
#include <nucleus/utils/ColourTexture.h>
#include <nucleus/utils/terrain_mesh_index_generator.h>
#include <nucleus/vector_layer/Preprocessor.h>
#include <nucleus/vector_layer/constants.h>
#include <radix/TileHeights.h>
#include <radix/quad_tree.h>
#include <radix/tile.h>

#include <nucleus/vector_layer/Style.h>
#include <nucleus/vector_layer/setup.h>

namespace quad_tree = radix::quad_tree;
using gl_engine::Framebuffer;
using gl_engine::ShaderProgram;
using nucleus::srs::hash_uint16;
using nucleus::srs::pack;
using nucleus::srs::unpack;
using namespace nucleus::tile;
using nucleus::tile::utils::AabbDecorator;
using nucleus::tile::utils::refineFunctor;
using radix::TileHeights;

namespace {

void clear_buffer()
{
    QOpenGLExtraFunctions* f = QOpenGLContext::currentContext()->extraFunctions();

    f->glClearColor(0.0, 0.0, 0.0, 1.0);
    f->glClear(GL_COLOR_BUFFER_BIT);

    // Clear Albedo-Buffer
    const GLfloat clearAlbedoColor[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    f->glClearBufferfv(GL_COLOR, 0, clearAlbedoColor);
    // Clear Position-Buffer (IMPORTANT [4] to <0, such that i know by sign if fragment was processed)
    const GLfloat clearPositionColor[4] = { 0.0f, 0.0f, 0.0f, -1.0f };
    f->glClearBufferfv(GL_COLOR, 1, clearPositionColor);
    // Clear Normals-Buffer
    const GLuint clearNormalColor[2] = { 0u, 0u };
    f->glClearBufferuiv(GL_COLOR, 2, clearNormalColor);
    // Clear Encoded-Depth Buffer
    const GLfloat clearEncDepthColor[] = { 0.0f, 0.0f, 0.0f, 0.0f };
    f->glClearBufferfv(GL_COLOR, 3, clearEncDepthColor);
    // Clear Depth-Buffer
    f->glClearDepthf(0.0f); // reverse z
    f->glClear(GL_DEPTH_BUFFER_BIT);
}

// if you need/want other settings probably best to add an enum argument and add an if so you can freely choose which settings to apply
void apply_opengl_render_settings()
{
    QOpenGLExtraFunctions* f = QOpenGLContext::currentContext()->extraFunctions();

    f->glEnable(GL_CULL_FACE);
    f->glCullFace(GL_BACK);

    f->glEnable(GL_DEPTH_TEST);
    f->glDepthFunc(GL_GREATER); // reverse z
}

nucleus::camera::Definition create_camera_for_id(nucleus::tile::Id id)
{
    // auto neighbour_id = nucleus::tile::Id { id.zoom_level, { id.coords.x + 0, id.coords.y }, id.scheme };
    const auto centre = nucleus::srs::tile_bounds(id).centre();

    // 34816000 / pow(2, id.zoom_level) -> should guarantee that a 2 px black border is around a single tile of the given zoom level
    // nucleus::camera::Definition camera = { { centre.x, centre.y, 34816000 / pow(2, id.zoom_level) }, { centre.x, centre.y, 0 } };
    nucleus::camera::Definition camera = { { centre.x, centre.y, 34816000 / pow(2, id.zoom_level) }, { centre.x, centre.y, 0 } };

    // camera.set_viewport_size({ 256, 256 });
    camera.set_viewport_size({ 1024, 1024 });
    camera.set_field_of_view(60); // FOV 60 is the default fov used
    return camera;
}

std::vector<nucleus::tile::Id> calc_children_ids(const nucleus::tile::Id& id, unsigned target_zoom)
{
    if (id.zoom_level >= target_zoom)
        return { id };

    auto children = id.children();
    auto c0 = calc_children_ids(children[0], target_zoom);
    auto c1 = calc_children_ids(children[1], target_zoom);
    auto c2 = calc_children_ids(children[2], target_zoom);
    auto c3 = calc_children_ids(children[3], target_zoom);

    std::vector<nucleus::tile::Id> out;

    out.insert(out.end(), c0.begin(), c0.end());
    out.insert(out.end(), c1.begin(), c1.end());
    out.insert(out.end(), c2.begin(), c2.end());
    out.insert(out.end(), c3.begin(), c3.end());

    return out;
}

// IMPORTANT: we need to return the shared_ptr, because otherwise it will be destroyed and nothing will be drawn!!!
std::shared_ptr<gl_engine::UniformBuffer<gl_engine::uboCameraConfig>> bind_camera(
    gl_engine::ShaderRegistry* shader_registry, const nucleus::camera::Definition& camera)
{
    auto camera_config_ubo = std::make_shared<gl_engine::UniformBuffer<gl_engine::uboCameraConfig>>(3, "camera_config");
    camera_config_ubo->init();
    camera_config_ubo->bind_to_shader(shader_registry->all());

    gl_engine::uboCameraConfig* cc = &camera_config_ubo->data;
    cc->position = glm::vec4(camera.position(), 1.0);
    cc->view_matrix = camera.local_view_matrix();
    cc->proj_matrix = camera.projection_matrix();
    cc->view_proj_matrix = cc->proj_matrix * cc->view_matrix;
    cc->inv_view_proj_matrix = glm::inverse(cc->view_proj_matrix);
    cc->inv_view_matrix = glm::inverse(cc->view_matrix);
    cc->inv_proj_matrix = glm::inverse(cc->proj_matrix);
    cc->viewport_size = camera.viewport_size();
    cc->distance_scaling_factor = camera.distance_scale_factor();
    // cc->error_threshold_px = camera.error_threshold_px();
    camera_config_ubo->update_gpu_data();

    return camera_config_ubo;
}

enum ConfigMode { no_overlay, cell_size, uvs, cells };

// IMPORTANT: we need to return the shared_ptr, because otherwise it will be destroyed and nothing will be drawn!!!
std::shared_ptr<gl_engine::UniformBuffer<gl_engine::uboSharedConfig>> bind_shared_config(gl_engine::ShaderRegistry* shader_registry, ConfigMode mode)
{
    auto shared_config_ubo = std::make_shared<gl_engine::UniformBuffer<gl_engine::uboSharedConfig>>(2, "shared_config");
    shared_config_ubo->init();
    shared_config_ubo->bind_to_shader(shader_registry->all());

    gl_engine::uboSharedConfig* config = &shared_config_ubo->data;

    switch (mode) {
    case cell_size:
        config->m_overlay_mode = 202u;
        break;
    case uvs:
        config->m_overlay_mode = 200u;
        break;
    case cells:
        config->m_overlay_mode = 204u;
        break;
    case no_overlay:
    default:
        config->m_overlay_mode = 0u;
    }

    shared_config_ubo->update_gpu_data();

    return shared_config_ubo;
}

std::shared_ptr<gl_engine::TileGeometry> load_tile_geometry(
    const nucleus::tile::utils::AabbDecoratorPtr& aabb_decorator, const nucleus::camera::Definition& camera)
{
    auto tile_geometry = std::make_shared<gl_engine::TileGeometry>();
    { // setup TileGeometry
        tile_geometry->set_aabb_decorator(aabb_decorator);
        tile_geometry->set_tile_limit(1024);
        tile_geometry->init();
    }

    auto geometry_service = std::make_unique<TileLoadService>(
        "https://alpinemaps.cg.tuwien.ac.at/tiles/at_dtm_alpinemaps/", nucleus::tile::TileLoadService::UrlPattern::ZXY, ".png");

    auto geometry = nucleus::tile::setup::geometry_scheduler(std::move(geometry_service), aabb_decorator);
    QSignalSpy spy_stats(geometry.scheduler.get(), &nucleus::tile::GeometryScheduler::stats_ready);
    QSignalSpy spy(geometry.scheduler.get(), &nucleus::tile::GeometryScheduler::gpu_tiles_updated);

    geometry.scheduler->set_enabled(true);
    QNetworkInformation* n = QNetworkInformation::instance();
    geometry.scheduler->set_network_reachability(n->reachability());
    geometry.scheduler->set_ram_quad_limit(1024);
    geometry.scheduler->set_name("geometry");
    geometry.scheduler->read_disk_cache();

    geometry.scheduler->update_camera(camera);

    geometry.scheduler->send_quad_requests();

    QVariantMap stats;
    do {
        spy_stats.wait(300000);
        QList<QVariant> arguments = spy_stats.takeFirst();
        stats = arguments.at(1).value<QVariantMap>();
    } while (stats.contains("n_quads_requested") && stats["n_quads_requested"] != 0);

    if (spy.empty())
        spy.wait(300000);

    REQUIRE(spy.count() >= 1);
    for (int i = 0; i < spy.count(); i++) {
        QList<QVariant> arguments = spy.takeFirst();
        REQUIRE(arguments.size() == 2);
        std::vector<nucleus::tile::Id> deleted_tiles = arguments.at(0).value<std::vector<nucleus::tile::Id>>();
        std::vector<nucleus::tile::GpuGeometryTile> actual_new_tiles = arguments.at(1).value<std::vector<nucleus::tile::GpuGeometryTile>>();

        std::vector<nucleus::tile::GpuGeometryTile> new_tiles;
        for (const auto& tile : actual_new_tiles) {
            // replace tile with 0 height tiles
            new_tiles.push_back({ tile.id, tile.bounds, std::make_shared<const nucleus::Raster<uint16_t>>(glm::uvec2(65), uint16_t(0)) });
            // new_tiles.push_back(tile);
        }

        tile_geometry->update_gpu_tiles(deleted_tiles, new_tiles);
    }

    geometry.scheduler->persist_tiles(); // not sure yet if it is good to persist tiles in a test/testing tool
    geometry.scheduler.reset();

    return tile_geometry;
}

std::shared_ptr<gl_engine::VectorLayer> create_vectorlayer(gl_engine::ShaderRegistry* shader_registry, const nucleus::vector_layer::Style& style)
{
    auto vectorlayer = std::make_shared<gl_engine::VectorLayer>();
    vectorlayer->set_tile_limit(1024);
    vectorlayer->init(shader_registry);
    vectorlayer->update_style(style.styles());

    return vectorlayer;
}

void load_custom_vectortile(std::shared_ptr<gl_engine::VectorLayer> vectorlayer, nucleus::vector_layer::Style style, nucleus::tile::Id id, QFile file)
{

    file.open(QFile::ReadOnly);
    const auto bytes = file.readAll();

    nucleus::vector_layer::Preprocessor preprocessor(std::move(style));
    auto tile = preprocessor.preprocess(id, bytes);

    // add the current tile to deleted_tiles to ensure that we override any existing data
    const std::vector<radix::tile::Id> deleted_tiles { id.parent() }; // TODO currently only whole quads can be deleted...
    std::vector<nucleus::tile::GpuVectorLayerTile> new_tiles;
    new_tiles.push_back(tile);

    vectorlayer->update_gpu_tiles(deleted_tiles, new_tiles);
}

void load_vectortiles(std::shared_ptr<gl_engine::VectorLayer> vectorlayer,
    const nucleus::tile::utils::AabbDecoratorPtr& aabb_decorator,
    const nucleus::camera::Definition& camera)
{

    auto vector_layer_service
        = std::make_unique<TileLoadService>("http://localhost:3000/tiles/", nucleus::tile::TileLoadService::UrlPattern::ZXY_yPointingSouth, "");
    auto processor = nucleus::vector_layer::setup::scheduler(std::move(vector_layer_service), aabb_decorator);
    QSignalSpy spy_stats(processor.scheduler.get(), &nucleus::vector_layer::Scheduler::stats_ready);
    QSignalSpy spy(processor.scheduler.get(), &nucleus::vector_layer::Scheduler::gpu_tiles_updated);

    processor.scheduler->set_enabled(true);
    // vector_layer.scheduler->set_update_timeout(500);
    QNetworkInformation* n = QNetworkInformation::instance();
    processor.scheduler->set_network_reachability(n->reachability());
    processor.scheduler->set_ram_quad_limit(1024);
    processor.scheduler->set_name("vector");
    processor.scheduler->read_disk_cache();

    processor.scheduler->update_camera(camera);
    processor.scheduler->send_quad_requests();

    QVariantMap stats;
    do {
        spy_stats.wait(300000);
        QList<QVariant> arguments = spy_stats.takeFirst();
        stats = arguments.at(1).value<QVariantMap>();
    } while (stats.contains("n_quads_requested") && stats["n_quads_requested"] != 0);

    if (spy.empty())
        spy.wait(300000);

    REQUIRE(spy.count() >= 1);

    for (int i = 0; i < spy.count(); i++) {
        QList<QVariant> arguments = spy.takeFirst();
        REQUIRE(arguments.size() == 2);
        std::vector<nucleus::tile::Id> deleted_tiles = arguments.at(0).value<std::vector<nucleus::tile::Id>>();
        std::vector<nucleus::tile::GpuVectorLayerTile> new_tiles = arguments.at(1).value<std::vector<nucleus::tile::GpuVectorLayerTile>>();

        vectorlayer->update_gpu_tiles(deleted_tiles, new_tiles);
    }

    processor.scheduler->persist_tiles(); // not sure yet if it is good to persist tiles in a test/testing tool
    processor.scheduler.reset();
}

template <typename LoadTileFunc>
void draw_tile(const nucleus::tile::utils::AabbDecoratorPtr& aabb_decorator,
    const nucleus::camera::Definition& camera,
    std::vector<nucleus::tile::Id> ids_to_draw,
    LoadTileFunc load,
    std::vector<std::pair<QString, ConfigMode>> outputs)
{
    // define the camera

    // define what ids will ultimately be drawn
    auto draw_list = nucleus::tile::drawing::compute_bounds(nucleus::tile::drawing::limit({ ids_to_draw }, 1024u), aabb_decorator);
    // auto draw_list = drawing::compute_bounds(drawing::limit(drawing::generate_list(camera, aabb_decorator, 19), 1024u), aabb_decorator);
    draw_list = drawing::sort(drawing::cull(draw_list, camera), camera.position());

    // loading tiles and upload to gpu
    auto shader_registry = gl_engine::ShaderRegistry();
    auto tile_geometry = load_tile_geometry(aabb_decorator, camera);

    auto layer = load(&shader_registry);

    // bind framebuffer
    Framebuffer b(Framebuffer::DepthFormat::Float32, { Framebuffer::ColourFormat::RGB8 }, { 1024, 1024 });
    b.bind();

    // bind buffer and set opengl settings
    apply_opengl_render_settings();
    auto camera_conf = bind_camera(&shader_registry, camera);

    for (const auto& output : outputs) {
        auto shared_conf = bind_shared_config(&shader_registry, output.second);
        clear_buffer();
        layer->draw(*tile_geometry, camera, draw_list);
        const QImage render_result = b.read_colour_attachment(0);
        render_result.save(output.first);
    }

    Framebuffer::unbind();
}

} // namespace

TEST_CASE("gl_engine/tile_drawing", "[!mayfail]")
{
    UnittestGLContext::initialise();

    // aabbdecorator -> makes sure that everything form 0-4000 z is included
    TileHeights h;
    h.emplace({ 0, { 0, 0 } }, { 0, 4000 });
    const auto aabb_decorator = AabbDecorator::make(std::move(h));

    // unittests
    SECTION("tile drawing kals")
    {

        auto id_kals = nucleus::tile::Id { 18, { 140276, 92195 }, nucleus::tile::Scheme::SlippyMap }.to(nucleus::tile::Scheme::Tms);
        const auto camera = create_camera_for_id(id_kals);

        auto load = [&id_kals](gl_engine::ShaderRegistry* shader_registry) {
            nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
            // qDebug() << "otm";
            style.load();
            // nucleus::vector_layer::Style style_qwant(":/vectorlayerstyles/qwant.json");
            // qDebug() << "qwant";
            // style_qwant.load();
            // nucleus::vector_layer::Style style_bright(":/vectorlayerstyles/osm-bright.json");
            // qDebug() << "bright";
            // style_bright.load();

            auto layer = create_vectorlayer(shader_registry, style);

            load_custom_vectortile(
                layer, std::move(style), id_kals, QFile(QString("%1%2").arg(ALP_GL_TEST_DATA_DIR, "vectortile_openmaptile_18_140276_92195_kals.pbf")));

            return layer;
        };

        auto outputs = std::vector<std::pair<QString, ConfigMode>> { std::make_pair("vectortile_kals.png", ConfigMode::no_overlay),
            std::make_pair("vectortile_kals_cells.png", ConfigMode::cells) };

        draw_tile(aabb_decorator, camera, { id_kals }, load, outputs);
    }

    // SECTION("tile drawing airport")
    // {

    //     auto id_airport = nucleus::tile::Id { 18, { 143144, 90996 }, nucleus::tile::Scheme::SlippyMap }.to(nucleus::tile::Scheme::Tms);
    //     // auto id_airport = nucleus::tile::Id { 18, { 143144, 90996 }, nucleus::tile::Scheme::SlippyMap }.to(nucleus::tile::Scheme::Tms).parent().parent();
    //     // auto id_airport = nucleus::tile::Id { 14, { 8936, 5681 }, nucleus::tile::Scheme::SlippyMap }.to(nucleus::tile::Scheme::Tms);

    //     const auto camera = create_camera_for_id(id_airport);

    //     auto load = [&camera, &aabb_decorator](gl_engine::ShaderRegistry* shader_registry) {
    //         nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
    //         style.load();

    //         auto layer = create_vectorlayer(shader_registry, style);

    //         load_vectortiles(layer, aabb_decorator, camera);
    //         return layer;
    //     };

    //     auto outputs = std::vector<std::pair<QString, ConfigMode>> { std::make_pair("vectortile_airport.png", ConfigMode::no_overlay),
    //         std::make_pair("vectortile_airport_cells.png", ConfigMode::cells) };

    //     draw_tile(aabb_decorator, camera, { id_airport }, load, outputs);
    //     // draw_tile(aabb_decorator, camera, drawing::generate_list(camera, aabb_decorator, 19), load, outputs);
    // }

    SECTION("tile drawing skilift")
    {

        auto id_ski = nucleus::tile::Id { 18, { 140269, 92190 }, nucleus::tile::Scheme::SlippyMap }.to(nucleus::tile::Scheme::Tms);
        const auto camera = create_camera_for_id(id_ski.children()[2].children()[2].children()[0]);

        auto load = [&id_ski](gl_engine::ShaderRegistry* shader_registry) {
            nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
            style.load();

            auto layer = create_vectorlayer(shader_registry, style);

            load_custom_vectortile(layer, std::move(style), id_ski, QFile(QString("%1%2").arg(ALP_GL_TEST_DATA_DIR, "vectortile_skilift_18_140269_92190.pbf")));

            return layer;
        };

        auto outputs = std::vector<std::pair<QString, ConfigMode>> { std::make_pair("vectortile_ski.png", ConfigMode::no_overlay),
            std::make_pair("vectortile_ski_cells.png", ConfigMode::cells) };

        draw_tile(aabb_decorator, camera, { id_ski }, load, outputs);
    }

    SECTION("tile drawing mountain")
    {

        auto id_mountain = nucleus::tile::Id { 14, { 8781, 5760 }, nucleus::tile::Scheme::SlippyMap }.to(nucleus::tile::Scheme::Tms);
        const auto camera = create_camera_for_id(id_mountain);

        auto load = [&id_mountain](gl_engine::ShaderRegistry* shader_registry) {
            nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
            style.load();

            auto layer = create_vectorlayer(shader_registry, style);

            load_custom_vectortile(
                layer, std::move(style), id_mountain, QFile(QString("%1%2").arg(ALP_GL_TEST_DATA_DIR, "vectortile_mountain_14_8781_5760.pbf")));

            return layer;
        };

        auto outputs = std::vector<std::pair<QString, ConfigMode>> { std::make_pair("vectortile_mountain.png", ConfigMode::no_overlay) };

        draw_tile(aabb_decorator, camera, { id_mountain }, load, outputs);
    }

    SECTION("tile drawing vienna")
    {
        auto id_wien = nucleus::tile::Id { 14, { 8936, 5681 }, nucleus::tile::Scheme::SlippyMap }.to(nucleus::tile::Scheme::Tms);
        const auto camera = create_camera_for_id(id_wien);

        auto load = [&camera, &aabb_decorator](gl_engine::ShaderRegistry* shader_registry) {
            nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
            style.load();

            auto layer = create_vectorlayer(shader_registry, style);

            load_vectortiles(layer, aabb_decorator, camera);

            return layer;
        };

        auto outputs = std::vector<std::pair<QString, ConfigMode>> { std::make_pair("vectortile_vienna.png", ConfigMode::no_overlay),
            std::make_pair("vectortile_vienna_uvs.png", ConfigMode::uvs) };

        draw_tile(aabb_decorator, camera, calc_children_ids(id_wien, 16), load, outputs);
    }
}
