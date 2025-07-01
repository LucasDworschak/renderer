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
#include <gl_engine/TextureLayer.h>
#include <gl_engine/TileGeometry.h>
#include <gl_engine/UniformBuffer.h>
#include <gl_engine/UniformBufferObjects.h>
#include <gl_engine/VectorLayer.h>
#include <gl_engine/helpers.h>
#include <nucleus/camera/PositionStorage.h>
#include <nucleus/srs.h>
#include <nucleus/tile/conversion.h>
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

enum DrawMode { NoOverlay, CellSize, UVs, Cells, FloatZoom };

struct DrawConfig {
    QString name;
    DrawMode mode;
};

struct Config {
    bool real_heights = false;
    std::vector<DrawConfig> draw_config;
};

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


// IMPORTANT: we need to return the shared_ptr, because otherwise it will be destroyed and nothing will be drawn!!!
std::shared_ptr<gl_engine::UniformBuffer<gl_engine::uboSharedConfig>> bind_shared_config(gl_engine::ShaderRegistry* shader_registry, DrawMode mode)
{
    auto shared_config_ubo = std::make_shared<gl_engine::UniformBuffer<gl_engine::uboSharedConfig>>(2, "shared_config");
    shared_config_ubo->init();
    shared_config_ubo->bind_to_shader(shader_registry->all());

    gl_engine::uboSharedConfig* config = &shared_config_ubo->data;

    switch (mode) {
    case CellSize:
        config->m_overlay_mode = 202u;
        break;
    case UVs:
        config->m_overlay_mode = 200u;
        break;
    case Cells:
        config->m_overlay_mode = 204u;
        break;
    case FloatZoom:
        config->m_overlay_mode = 209u;
        break;
    case NoOverlay:
    default:
        config->m_overlay_mode = 0u;
    }

    shared_config_ubo->update_gpu_data();

    return shared_config_ubo;
}
template <typename LoadFunc>
void loading(LoadFunc load_func, QSignalSpy& spy_stats, QSignalSpy& spy_tile)
{
    spy_stats.wait(1000);
    QVariantMap stats;
    do {
        if (spy_stats.isEmpty())
            spy_stats.wait(10000);
        REQUIRE(spy_stats.count() >= 1);
        QList<QVariant> arguments = spy_stats.takeLast();
        stats = arguments.at(1).value<QVariantMap>();
        spy_stats.wait(1000);
    } while (!spy_stats.isEmpty() || (stats.contains("n_quads_requested") && stats["n_quads_requested"] != 0));

    spy_tile.wait(1000);
    do {
        for (int i = 0; i < spy_tile.count(); i++) {
            QList<QVariant> arguments = spy_tile.takeFirst();
            REQUIRE(arguments.size() == 2);

            load_func(arguments);
        }
        spy_tile.wait(1000);
    } while (!spy_tile.isEmpty());
}

std::shared_ptr<gl_engine::TileGeometry> load_tile_geometry(
    const nucleus::tile::utils::AabbDecoratorPtr& aabb_decorator, const nucleus::camera::Definition& camera, bool real_heights)
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

    auto load_func = [&tile_geometry, &real_heights](QList<QVariant> arguments) {
        std::vector<nucleus::tile::Id> deleted_tiles = arguments.at(0).value<std::vector<nucleus::tile::Id>>();
        std::vector<nucleus::tile::GpuGeometryTile> actual_new_tiles = arguments.at(1).value<std::vector<nucleus::tile::GpuGeometryTile>>();

        if (real_heights) {
            tile_geometry->update_gpu_tiles(deleted_tiles, actual_new_tiles);
        } else {
            std::vector<nucleus::tile::GpuGeometryTile> new_tiles;
            for (const auto& tile : actual_new_tiles) {
                // replace tile with 0 height tiles
                new_tiles.push_back({ tile.id, tile.bounds, std::make_shared<const nucleus::Raster<uint16_t>>(glm::uvec2(65), uint16_t(0)) });
            }

            tile_geometry->update_gpu_tiles(deleted_tiles, new_tiles);
        }
    };

    loading(load_func, spy_stats, spy);

    geometry.scheduler->persist_tiles(); // not sure yet if it is good to persist tiles in a test/testing tool
    geometry.scheduler.reset();

    return tile_geometry;
}

std::shared_ptr<gl_engine::TextureLayer> load_texture_layer(
    gl_engine::ShaderRegistry* shader_registry, const nucleus::tile::utils::AabbDecoratorPtr& aabb_decorator, const nucleus::camera::Definition& camera)
{
    auto texture_layer = std::make_shared<gl_engine::TextureLayer>(512);
    texture_layer->set_tile_limit(1024);
    texture_layer->init(shader_registry);

    auto ortho_service = std::make_unique<TileLoadService>(
        "https://mapsneu.wien.gv.at/basemap/bmapoberflaeche/grau/google3857/", nucleus::tile::TileLoadService::UrlPattern::ZYX_yPointingSouth, ".jpeg");

    auto ortho_texture = nucleus::tile::setup::texture_scheduler(std::move(ortho_service), aabb_decorator);

    QSignalSpy spy_stats(ortho_texture.scheduler.get(), &nucleus::tile::TextureScheduler::stats_ready);
    QSignalSpy spy(ortho_texture.scheduler.get(), &nucleus::tile::TextureScheduler::gpu_tiles_updated);

    ortho_texture.scheduler->set_enabled(true);
    QNetworkInformation* n = QNetworkInformation::instance();
    ortho_texture.scheduler->set_network_reachability(n->reachability());
    ortho_texture.scheduler->set_ram_quad_limit(1024);
    ortho_texture.scheduler->set_name("ortho");
    ortho_texture.scheduler->read_disk_cache();

    ortho_texture.scheduler->update_camera(camera);
    ortho_texture.scheduler->send_quad_requests();

    auto load_func = [&texture_layer](QList<QVariant> arguments) {
        std::vector<nucleus::tile::Id> deleted_tiles = arguments.at(0).value<std::vector<nucleus::tile::Id>>();
        std::vector<GpuTextureTile> actual_new_tiles = arguments.at(1).value<std::vector<GpuTextureTile>>();

        // if (real_ortho) {
        texture_layer->update_gpu_tiles(deleted_tiles, actual_new_tiles);
        // } else {
        //     std::vector<GpuTextureTile> new_tiles;
        //     for (const auto& tile : actual_new_tiles) {
        //         // replace tile with 0 height tiles
        //         new_tiles.push_back({ tile.id, std::make_shared<const nucleus::Raster<uint16_t>>(glm::uvec2(256), uint16_t(0)) });
        //     }

        //     texture_layer->update_gpu_tiles(deleted_tiles, new_tiles);
        // }
    };

    loading(load_func, spy_stats, spy);

    ortho_texture.scheduler->persist_tiles(); // not sure yet if it is good to persist tiles in a test/testing tool
    ortho_texture.scheduler.reset();

    return texture_layer;
}

std::shared_ptr<gl_engine::VectorLayer> create_vectorlayer(gl_engine::ShaderRegistry* shader_registry, const nucleus::vector_layer::Style& style)
{
    auto vectorlayer = std::make_shared<gl_engine::VectorLayer>();

    auto defines_map = gl_engine::VectorLayer::default_defines();
    defines_map[QString("display_mode")] = QString::number(1);

    vectorlayer->set_defines(defines_map);

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

    auto vector_layer_service = std::make_unique<TileLoadService>(
        "https://osm.cg.tuwien.ac.at/vector_tiles/vector_layer_v1/", nucleus::tile::TileLoadService::UrlPattern::ZXY_yPointingSouth, "");
    auto processor = nucleus::vector_layer::setup::scheduler(std::move(vector_layer_service), aabb_decorator);

    QSignalSpy spy_stats(processor.scheduler.get(), &nucleus::vector_layer::Scheduler::stats_ready);
    QSignalSpy spy(processor.scheduler.get(), &nucleus::vector_layer::Scheduler::gpu_tiles_updated);

    // vector_layer.scheduler->set_update_timeout(500);
    QNetworkInformation* n = QNetworkInformation::instance();
    processor.scheduler->set_network_reachability(n->reachability());
    processor.scheduler->set_ram_quad_limit(1024);
    processor.scheduler->set_name("vector");
    processor.scheduler->read_disk_cache();
    processor.scheduler->set_enabled(true);

    processor.scheduler->update_camera(camera);
    processor.scheduler->send_quad_requests();

    auto load_func = [&vectorlayer](QList<QVariant> arguments) {
        std::vector<nucleus::tile::Id> deleted_tiles = {};
        std::vector<nucleus::tile::GpuVectorLayerTile> new_tiles = arguments.at(1).value<std::vector<nucleus::tile::GpuVectorLayerTile>>();

        vectorlayer->update_gpu_tiles(deleted_tiles, new_tiles);
    };

    loading(load_func, spy_stats, spy);

    processor.scheduler->persist_tiles(); // not sure yet if it is good to persist tiles in a test/testing tool
    processor.scheduler.reset();
}

void compare_images(QImage render_result, QString name)
{
    auto file = QFile(QString("%1%2").arg(ALP_GL_TEST_DATA_DIR, name + ".png"));

    file.open(QFile::ReadOnly);
    const auto bytes = file.readAll();
    auto comparison_image = QImage::fromData(bytes);

    if (comparison_image.isNull()) {
        qDebug() << "no comparison image found for" << name << "-> only writing out result";
        render_result.save(name + ".png");

        CHECK(false);
        return;
    }

    // convert to same format
    comparison_image = comparison_image.convertedTo(QImage::Format_RGBA8888);
    render_result = render_result.convertedTo(QImage::Format_RGBA8888);

    // determine if render result is equal to the comparison image
    CHECK(comparison_image == render_result);

    if (comparison_image != render_result) {

        // make difference comparison_image
        // use the comparison image as the input and output of the diff raster

        // output the current reference image -> easier to compare between render and reference
        comparison_image.save(name + "_reference.png");

        auto diff_raster = nucleus::tile::conversion::to_rgba8raster(comparison_image);
        auto render_raster = nucleus::tile::conversion::to_rgba8raster(render_result);

        const auto w0 = 0;
        const auto h0 = 0;
        const auto w1 = diff_raster.width();
        const auto h1 = diff_raster.height();
        for (unsigned x = w0; x < w1; x++) {
            for (unsigned y = h0; y < h1; y++) {
                diff_raster.pixel({ x, y }) = abs(diff_raster.pixel({ x, y }) - render_raster.pixel({ x, y }));
                diff_raster.pixel({ x, y }).w = 255; // ensure that we have an opaque image -> easier to see differences
            }
        }

        auto diff_image = nucleus::tile::conversion::to_QImage(diff_raster);

        // save the current render result
        diff_image.save(name + "_diff.png");

        render_result.save(name + "_render.png");
    }
}

template <typename LoadTileFunc>
void draw_tile(const nucleus::tile::utils::AabbDecoratorPtr& aabb_decorator,
    const nucleus::camera::Definition& camera,
    std::vector<nucleus::tile::Id> ids_to_draw,
    LoadTileFunc load,
    Config config)
{
    // define the camera

    // define what ids will ultimately be drawn
    auto draw_list = nucleus::tile::drawing::compute_bounds(nucleus::tile::drawing::limit(ids_to_draw, 1024u), aabb_decorator);
    // auto draw_list = drawing::compute_bounds(drawing::limit(drawing::generate_list(camera, aabb_decorator, 19), 1024u), aabb_decorator);
    draw_list = drawing::sort(drawing::cull(draw_list, camera), camera.position());

    // loading tiles and upload to gpu
    auto shader_registry = gl_engine::ShaderRegistry();
    auto tile_geometry = load_tile_geometry(aabb_decorator, camera, config.real_heights);
    auto texture_layer = load_texture_layer(&shader_registry, aabb_decorator, camera);

    auto layer = load(&shader_registry);

    // bind framebuffer
    // Framebuffer b(Framebuffer::DepthFormat::Float32, { Framebuffer::ColourFormat::RGB8 }, { 1024, 1024 });
    Framebuffer b(Framebuffer::DepthFormat::Float32,
        std::vector {
            Framebuffer::ColourFormat::RGB8, // Albedo
            // Framebuffer::ColourFormat::RGBA32F, // Position WCS and distance (distance is optional, but i use it directly for a little speed improvement)
            // Framebuffer::ColourFormat::RG16UI, // Octahedron Normals
            // Framebuffer::ColourFormat::RGBA8, // Discretized Encoded Depth for readback IMPORTANT: IF YOU MOVE THIS YOU HAVE TO ADAPT THE GET DEPTH FUNCTION
            // TextureDefinition { Framebuffer::ColourFormat::R32UI }, // VertexID
            // Framebuffer::ColourFormat::RGB8, // vector map colour
            // Framebuffer::ColourFormat::RGBA8, // vector map colour
        },
        { camera.viewport_size() });

    b.bind();

    // bind buffer and set opengl settings
    apply_opengl_render_settings();
    auto camera_conf = bind_camera(&shader_registry, camera);

    for (const auto& config : config.draw_config) {
        auto shared_conf = bind_shared_config(&shader_registry, config.mode);
        clear_buffer();
        layer->draw(*tile_geometry, *texture_layer, camera, draw_list);
        const QImage render_result = b.read_colour_attachment(0);

        compare_images(render_result, config.name);
    }

    Framebuffer::unbind();
}

} // namespace

TEST_CASE("gl_engine/tile_drawing", "[!mayfail]")
{
    SECTION("refine ids")
    {
        auto id_wien = nucleus::tile::Id { 14, { 8936, 5681 }, nucleus::tile::Scheme::SlippyMap }.to(nucleus::tile::Scheme::Tms);

        auto children = calc_children_ids(id_wien, 16);

        CHECK(children.size() == 16);
    }

    UnittestGLContext::initialise();

    // aabbdecorator -> makes sure that everything form 0-4000 z is included
    TileHeights h;
    h.emplace({ 0, { 0, 0 } }, { 0, 4000 });
    const auto aabb_decorator = AabbDecorator::make(std::move(h));

    // SECTION("draw custom cell")
    // {
    //     auto id_ski = nucleus::tile::Id { 18, { 140269, 92190 }, nucleus::tile::Scheme::SlippyMap }.to(nucleus::tile::Scheme::Tms);
    //     // const auto camera = create_camera_for_id(id_ski.children()[2].children()[2].children()[2]);
    //     const auto camera = create_camera_for_id(id_ski.children()[0].children()[2].children()[0]);
    //     // const auto camera = create_camera_for_id(id_ski);

    //     auto load = [&id_ski](gl_engine::ShaderRegistry* shader_registry) {
    //         nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
    //         style.load();

    //         auto layer = create_vectorlayer(shader_registry, style);

    //         nucleus::vector_layer::Preprocessor preprocessor(std::move(style));
    //         // auto tile = preprocessor.preprocess(id_ski, bytes);

    //         const nucleus::vector_layer::ClipperPaths polygon = { {
    //             { 156 * nucleus::vector_layer::constants::scale_polygons, 912 * nucleus::vector_layer::constants::scale_polygons },
    //             { 68 * nucleus::vector_layer::constants::scale_polygons, 616 * nucleus::vector_layer::constants::scale_polygons },
    //             { 150 * nucleus::vector_layer::constants::scale_polygons, 616 * nucleus::vector_layer::constants::scale_polygons },
    //         } };
    //         // const nucleus::vector_layer::ClipperPaths polygon = { {
    //         //     { 72, 631 },
    //         //     { 81, 662 },
    //         //     { 82, 631 },
    //         // } };
    //         // const nucleus::vector_layer::ClipperPaths polygon = { { { 300, 300 }, { 500, 500 }, { 500, 800 } } };
    //         // const nucleus::vector_layer::ClipperPaths polygon = { { { 48 + 16 + 4, 48 - 4 }, { 48 - 4, 48 - 4 }, { 48 + 8, 48 + 16 + 4 } } };
    //         // const nucleus::vector_layer::ClipperPaths polygon = { { { 48 + 16, 48 }, { 48, 48 }, { 48 + 8, 48 + 16 } } };

    //         // auto b = int16_t(nucleus::vector_layer::constants::aa_border);
    //         // const nucleus::vector_layer::ClipperPaths polygon = { { { 48 + 16 - b, 48 + b }, { 48 + b, 48 + b }, { 48 + 8, 48 + 16 - b } } };

    //         nucleus::vector_layer::VectorLayers tile_data;
    //         std::vector<nucleus::vector_layer::ClipperRect> bounds;
    //         radix::geometry::Aabb2i aabb({});
    //         constexpr float scale = nucleus::vector_layer::constants::scale_polygons;
    //         constexpr auto cell_scale = float(nucleus::vector_layer::constants::grid_size) / (float(nucleus::vector_layer::constants::tile_extent) * scale);

    //         bounds.reserve(polygon.size());
    //         for (size_t j = 0; j < polygon.size(); j++) {
    //             const auto bound = Clipper2Lib::GetBounds(polygon[j]);

    //             aabb.expand_by({ std::floor(float(bound.left) * cell_scale), std::floor(float(bound.top) * cell_scale) });
    //             aabb.expand_by({ std::ceil(float(bound.right) * cell_scale), std::ceil(float(bound.bottom) * cell_scale) });
    //             bounds.push_back(bound);
    //         }

    //         std::pair<uint32_t, uint32_t> style_layer = { 148u << 1, 0u };
    //         // std::pair<uint32_t, uint32_t> style_layer = { 0, 0 };
    //         tile_data[style_layer.second].emplace_back(polygon, bounds, aabb, style_layer, true);

    //         preprocessor.preprocess_geometry(tile_data);
    //         auto tile = preprocessor.create_gpu_tile();
    //         tile.id = id_ski;

    //         // add the current tile to deleted_tiles to ensure that we override any existing data
    //         const std::vector<radix::tile::Id> deleted_tiles {};
    //         std::vector<nucleus::tile::GpuVectorLayerTile> new_tiles;
    //         new_tiles.push_back(tile);

    //         layer->update_gpu_tiles(deleted_tiles, new_tiles);

    //         return layer;
    //     };

    //     auto config = Config { false, std::vector<DrawConfig>{{"vectortile_custom_cell", DrawMode::NoOverlay }}};

    //     draw_tile(aabb_decorator, camera, { id_ski }, load, config);
    // }

    // unittests
    SECTION("tile drawing straight thin line")
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

        auto config = Config { false, std::vector<DrawConfig> { { "vectortile_straight_thin_line", DrawMode::NoOverlay } } };

        draw_tile(aabb_decorator, camera, { id_ski }, load, config);
    }

    SECTION("tile drawing straight polygon")
    {
        auto id_ski = nucleus::tile::Id { 18, { 140269, 92190 }, nucleus::tile::Scheme::SlippyMap }.to(nucleus::tile::Scheme::Tms);
        const auto camera = create_camera_for_id(id_ski.children()[0].children()[2].children()[0]);

        auto load = [&id_ski](gl_engine::ShaderRegistry* shader_registry) {
            nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
            style.load();

            auto layer = create_vectorlayer(shader_registry, style);
            load_custom_vectortile(layer, std::move(style), id_ski, QFile(QString("%1%2").arg(ALP_GL_TEST_DATA_DIR, "vectortile_skilift_18_140269_92190.pbf")));

            return layer;
        };

        auto config = Config { false, std::vector<DrawConfig> { { "vectortile_straight_polygon", DrawMode::NoOverlay } } };

        draw_tile(aabb_decorator, camera, { id_ski }, load, config);
    }

    SECTION("tile drawing vienna uniform tiles")
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

        auto config = Config { false, std::vector<DrawConfig> { { "vectortile_vienna_uniform", DrawMode::NoOverlay } } };

        draw_tile(aabb_decorator, camera, calc_children_ids(id_wien, 16), load, config);
    }

    SECTION("tile drawing vienna view")
    {
        auto camera = nucleus::camera::stored_positions::wien();
        camera.set_viewport_size({ 1920, 1080 });
        camera.set_field_of_view(60); // FOV 60 is the default fov used

        auto load = [&camera, &aabb_decorator](gl_engine::ShaderRegistry* shader_registry) {
            nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
            style.load();

            auto layer = create_vectorlayer(shader_registry, style);
            load_vectortiles(layer, aabb_decorator, camera);

            return layer;
        };

        auto config = Config { true, std::vector<DrawConfig> { { "vectortile_vienna_view", DrawMode::NoOverlay } } };

        const auto draw_list = drawing::generate_list(camera, aabb_decorator, 19);

        draw_tile(aabb_decorator, camera, draw_list, load, config);
    }

    SECTION("tile drawing grossglockner view")
    {
        auto camera = nucleus::camera::stored_positions::grossglockner();
        camera.set_viewport_size({ 1920, 1080 });
        camera.set_field_of_view(60); // FOV 60 is the default fov used

        auto load = [&camera, &aabb_decorator](gl_engine::ShaderRegistry* shader_registry) {
            nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
            style.load();

            auto layer = create_vectorlayer(shader_registry, style);
            load_vectortiles(layer, aabb_decorator, camera);

            return layer;
        };

        auto config = Config { true, std::vector<DrawConfig> { { "vectortile_grossglockner_view", DrawMode::NoOverlay } } };

        const auto draw_list = drawing::generate_list(camera, aabb_decorator, 19);

        draw_tile(aabb_decorator, camera, draw_list, load, config);
    }

    SECTION("tile drawing weichtalhaus view")
    {
        auto camera = nucleus::camera::stored_positions::weichtalhaus();
        camera.set_viewport_size({ 1920, 1080 });
        camera.set_field_of_view(60); // FOV 60 is the default fov used

        auto load = [&camera, &aabb_decorator](gl_engine::ShaderRegistry* shader_registry) {
            nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
            style.load();

            auto layer = create_vectorlayer(shader_registry, style);
            load_vectortiles(layer, aabb_decorator, camera);

            return layer;
        };

        auto config = Config { true, std::vector<DrawConfig> { { "vectortile_weichtalhaus_view", DrawMode::NoOverlay } } };

        const auto draw_list = drawing::generate_list(camera, aabb_decorator, 19);

        draw_tile(aabb_decorator, camera, draw_list, load, config);
    }

    SECTION("tile drawing landstrasse view")
    {
        const auto coords_lookat = nucleus::srs::lat_long_alt_to_world({ 48.20618344255035, 16.385128312902415, 0 });
        const auto coords_position = nucleus::srs::lat_long_alt_to_world({ 48.20618344255035, 16.385128312902415, 200 });
        auto camera = nucleus::camera::Definition { coords_position, coords_lookat };

        camera.set_viewport_size({ 1920, 1080 });
        camera.set_field_of_view(60); // FOV 60 is the default fov used

        auto load = [&camera, &aabb_decorator](gl_engine::ShaderRegistry* shader_registry) {
            nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
            style.load();

            auto layer = create_vectorlayer(shader_registry, style);
            load_vectortiles(layer, aabb_decorator, camera);

            return layer;
        };

        auto config = Config { true, std::vector<DrawConfig> { { "vectortile_landstrasse_view", DrawMode::NoOverlay } } };

        const auto draw_list = drawing::generate_list(camera, aabb_decorator, 19);

        draw_tile(aabb_decorator, camera, draw_list, load, config);
    }

    SECTION("tile drawing highway view")
    {
        nucleus::camera::Definition camera = { { 1.81929e+06, 6.11431e+06, 1671.63 }, { 1.81929e+06, 6.11431e+06 + 100, 1671.63 - 100 } };

        camera.set_viewport_size({ 1920, 1080 });
        camera.set_field_of_view(60); // FOV 60 is the default fov used

        auto load = [&camera, &aabb_decorator](gl_engine::ShaderRegistry* shader_registry) {
            nucleus::vector_layer::Style style(":/vectorlayerstyles/openstreetmap.json");
            style.load();

            auto layer = create_vectorlayer(shader_registry, style);
            load_vectortiles(layer, aabb_decorator, camera);

            return layer;
        };

        auto config = Config { false, std::vector<DrawConfig> { { "vectortile_highway", DrawMode::NoOverlay } } };

        const auto draw_list = drawing::generate_list(camera, aabb_decorator, 19);

        draw_tile(aabb_decorator, camera, draw_list, load, config);
    }
}
