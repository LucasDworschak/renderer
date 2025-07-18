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
#include "UnittestGLContext.h"

#include <unordered_set>
#include <vector>

#include <QImage>
#include <catch2/catch_test_macros.hpp>

#include <gl_engine/Framebuffer.h>
#include <gl_engine/ShaderProgram.h>
#include <gl_engine/Texture.h>
#include <gl_engine/VectorLayer.h>
#include <gl_engine/helpers.h>
#include <nucleus/camera/PositionStorage.h>
#include <nucleus/srs.h>
#include <nucleus/tile/utils.h>
#include <nucleus/utils/ColourTexture.h>
#include <nucleus/vector_layer/Preprocessor.h>
#include <nucleus/vector_layer/constants.h>
#include <radix/TileHeights.h>
#include <radix/quad_tree.h>
#include <radix/tile.h>

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
const QString default_vertex_shader = R"(
out highp vec2 texcoords;
void main() {
    vec2 vertices[3]=vec2[3](vec2(-1.0, -1.0), vec2(3.0, -1.0), vec2(-1.0, 3.0));
    gl_Position = vec4(vertices[gl_VertexID], 0.0, 1.0);
    texcoords = 0.5 * gl_Position.xy + vec2(0.5);
})";

ShaderProgram create_debug_shader(
    const QString& fragment_shader, const QString& vertex_shader = default_vertex_shader, const std::vector<QString>& defines = {})
{
    ShaderProgram tmp(vertex_shader, fragment_shader, gl_engine::ShaderCodeSource::PLAINTEXT, defines);
    return tmp;
}

void hashing_cpp_same_as_glsl(const Id& id)
{
    {
        Framebuffer b(Framebuffer::DepthFormat::None, { Framebuffer::ColourFormat::RGBA8 }, { 1, 1 });
        b.bind();

        ShaderProgram shader = create_debug_shader(QString(R"(
            #include "tile_id.glsl"

            out lowp vec4 out_color;
            void main() {
                mediump uint hash_ref = %1u;
                mediump uint hash_glsl = hash_tile_id(uvec3(%2u, %3u, %4u));
                out_color = vec4((hash_ref == hash_glsl) ? 121.0 / 255.0 : 9.0 / 255.0, 0, 0, 1);
            }
        )")
                .arg(hash_uint16(id))
                .arg(id.coords.x)
                .arg(id.coords.y)
                .arg(id.zoom_level));
        shader.bind();
        gl_engine::helpers::create_screen_quad_geometry().draw();

        const QImage render_result = b.read_colour_attachment(0);
        Framebuffer::unbind();
        CHECK(qRed(render_result.pixel(0, 0)) == 121);
    }
    {
        Framebuffer b(Framebuffer::DepthFormat::None, { Framebuffer::ColourFormat::RGBA8 }, { 1, 1 });
        b.bind();

        ShaderProgram shader = create_debug_shader(QString(R"(
            out lowp vec4 out_color;
            flat in mediump uint hash;
            void main() {
                mediump uint hash_ref = %1u;
                out_color = vec4((hash_ref == hash) ? 121.0 / 255.0 : 9.0 / 255.0, 0, 0, 1);
            }
            )")
                                                       .arg(hash_uint16(id)),
            QString(R"(
            #include "tile_id.glsl"

            out highp vec2 texcoords;
            flat out mediump uint hash;
            void main() {
                hash = hash_tile_id(uvec3(%1u, %2u, %3u));
                vec2 vertices[3]=vec2[3](vec2(-1.0, -1.0), vec2(3.0, -1.0), vec2(-1.0, 3.0));
                gl_Position = vec4(vertices[gl_VertexID], 0.0, 1.0);
                texcoords = 0.5 * gl_Position.xy + vec2(0.5);
            })")
                .arg(id.coords.x)
                .arg(id.coords.y)
                .arg(id.zoom_level));
        shader.bind();
        gl_engine::helpers::create_screen_quad_geometry().draw();

        const QImage render_result = b.read_colour_attachment(0);
        Framebuffer::unbind();
        CHECK(qRed(render_result.pixel(0, 0)) == 121);
    }
}

void packing_cpp_same_as_glsl(const Id& id)
{
    {
        Framebuffer b(Framebuffer::DepthFormat::None, { Framebuffer::ColourFormat::RGBA8 }, { 1, 1 });
        b.bind();

        ShaderProgram shader = create_debug_shader(QString(R"(
            #include "tile_id.glsl"

            out lowp vec4 out_color;
            void main() {
                highp uvec2 cpp_packed_id = uvec2(%1u, %2u);
                highp uvec3 id = uvec3(%3u, %4u, %5u);
                highp uvec3 unpacked_id = unpack_tile_id(cpp_packed_id);
                bool unpack_ok = id == unpacked_id;
                highp uvec2 packed_id = pack_tile_id(id);
                bool pack_ok = packed_id == cpp_packed_id;
                out_color = vec4(unpack_ok ? 121.0 / 255.0 : 9.0 / 255.0, pack_ok ? 122.0 / 255.0 : 9.0 / 255.0, 0, 1);
            }
        )")
                                                       .arg(pack(id).x)
                                                       .arg(pack(id).y)
                                                       .arg(id.coords.x)
                                                       .arg(id.coords.y)
                                                       .arg(id.zoom_level));
        shader.bind();
        gl_engine::helpers::create_screen_quad_geometry().draw();

        const QImage render_result = b.read_colour_attachment(0);
        Framebuffer::unbind();
        CHECK(qRed(render_result.pixel(0, 0)) == 121);
        CHECK(qGreen(render_result.pixel(0, 0)) == 122);
    }
    {
        Framebuffer b(Framebuffer::DepthFormat::None, { Framebuffer::ColourFormat::RGBA8 }, { 1, 1 });
        b.bind();

        ShaderProgram shader = create_debug_shader(QString(R"(
            out lowp vec4 out_color;
            flat in lowp vec2 ok;
            void main() {
                out_color = vec4(ok, 0, 1);
            }
            )"),
            QString(R"(
            #include "tile_id.glsl"

            out highp vec2 texcoords;
            flat out lowp vec2 ok;
            void main() {
                highp uvec2 cpp_packed_id = uvec2(%1u, %2u);
                highp uvec3 id = uvec3(%3u, %4u, %5u);
                highp uvec3 unpacked_id = unpack_tile_id(cpp_packed_id);
                bool unpack_ok = id == unpacked_id;
                highp uvec2 packed_id = pack_tile_id(id);
                bool pack_ok = packed_id == cpp_packed_id;
                ok = vec2(unpack_ok ? 121.0 / 255.0 : 9.0 / 255.0, pack_ok ? 122.0 / 255.0 : 9.0 / 255.0);

                vec2 vertices[3]=vec2[3](vec2(-1.0, -1.0), vec2(3.0, -1.0), vec2(-1.0, 3.0));
                gl_Position = vec4(vertices[gl_VertexID], 0.0, 1.0);
                texcoords = 0.5 * gl_Position.xy + vec2(0.5);
            })")
                .arg(pack(id).x)
                .arg(pack(id).y)
                .arg(id.coords.x)
                .arg(id.coords.y)
                .arg(id.zoom_level));
        shader.bind();
        gl_engine::helpers::create_screen_quad_geometry().draw();

        const QImage render_result = b.read_colour_attachment(0);
        Framebuffer::unbind();
        CHECK(qRed(render_result.pixel(0, 0)) == 121);
        CHECK(qGreen(render_result.pixel(0, 0)) == 122);
    }
}

void vectorlayer_packing_cpp_same_as_glsl(const nucleus::vector_layer::VectorLayerData& data)
{
    auto defines_map = gl_engine::VectorLayer::default_defines();
    defines_map[QString("display_mode")] = QString::number(1);
    std::vector<QString> defines;
    for (const auto& define : defines_map) {
        defines.push_back("#define " + define.first + " " + define.second);
    }

    {
        Framebuffer b(Framebuffer::DepthFormat::None, { Framebuffer::ColourFormat::RGBA8 }, { 1, 1 });
        b.bind();

        glm::uvec2 packed_data;
        if (data.is_polygon)
            packed_data = nucleus::vector_layer::Preprocessor::pack_triangle_data(data);
        else
            packed_data = nucleus::vector_layer::Preprocessor::pack_line_data(data.a, data.b, data.style_index);

        qDebug() << packed_data.x << packed_data.y;

        ShaderProgram shader = create_debug_shader(QString(R"(
            #include "vector_layer.glsl"

            out lowp vec4 out_color;
            void main() {
                highp uvec2 cpp_packed_data = uvec2(%1u, %2u);
                VectorLayerData data = VectorLayerData(ivec2(%3, %4),ivec2(%5, %6),ivec2(%7, %8), %9u, bool(%10));
                VectorLayerData unpacked_data = unpack_data(cpp_packed_data);

                bool unpack_ok = data == unpacked_data;
                highp uvec2 packed_data = pack_vectorlayer_data(data);
                bool pack_ok = packed_data == cpp_packed_data;
                out_color = vec4(unpack_ok ? 121.0 / 255.0 : 9.0 / 255.0, pack_ok ? 122.0 / 255.0 : 9.0 / 255.0, 0, 1);
            }
        )")
                                                       .arg(packed_data.x)
                                                       .arg(packed_data.y)
                                                       .arg(data.a.x)
                                                       .arg(data.a.y)
                                                       .arg(data.b.x)
                                                       .arg(data.b.y)
                                                       .arg(data.c.x)
                                                       .arg(data.c.y)
                                                       .arg(data.style_index)
                                                       .arg(data.is_polygon),
            default_vertex_shader,
            defines);
        shader.bind();
        gl_engine::helpers::create_screen_quad_geometry().draw();

        const QImage render_result = b.read_colour_attachment(0);
        Framebuffer::unbind();
        CHECK(qRed(render_result.pixel(0, 0)) == 121);
        CHECK(qGreen(render_result.pixel(0, 0)) == 122);
    }
    {
        Framebuffer b(Framebuffer::DepthFormat::None, { Framebuffer::ColourFormat::RGBA8 }, { 1, 1 });
        b.bind();

        glm::uvec2 packed_data;
        if (data.is_polygon)
            packed_data = nucleus::vector_layer::Preprocessor::pack_triangle_data(data);
        else
            packed_data = nucleus::vector_layer::Preprocessor::pack_line_data(data.a, data.b, data.style_index);

        ShaderProgram shader = create_debug_shader(QString(R"(
            out lowp vec4 out_color;
            flat in lowp vec2 ok;
            void main() {
                out_color = vec4(ok, 0, 1);
            }
            )"),
            QString(R"(
            #include "vector_layer.glsl"

            out highp vec2 texcoords;
            flat out lowp vec2 ok;
            void main() {
                highp uvec2 cpp_packed_data = uvec2(%1u, %2u);
                VectorLayerData data = VectorLayerData(ivec2(%3, %4),ivec2(%5, %6),ivec2(%7, %8), %9u, bool(%10));
                VectorLayerData unpacked_data = unpack_data(cpp_packed_data);

                bool unpack_ok = data == unpacked_data;
                highp uvec2 packed_data = pack_vectorlayer_data(data);
                bool pack_ok = packed_data == cpp_packed_data;
                ok = vec2(unpack_ok ? 121.0 / 255.0 : 9.0 / 255.0, pack_ok ? 122.0 / 255.0 : 9.0 / 255.0);

                vec2 vertices[3]=vec2[3](vec2(-1.0, -1.0), vec2(3.0, -1.0), vec2(-1.0, 3.0));
                gl_Position = vec4(vertices[gl_VertexID], 0.0, 1.0);
                texcoords = 0.5 * gl_Position.xy + vec2(0.5);
            })")
                .arg(packed_data.x)
                .arg(packed_data.y)
                .arg(data.a.x)
                .arg(data.a.y)
                .arg(data.b.x)
                .arg(data.b.y)
                .arg(data.c.x)
                .arg(data.c.y)
                .arg(data.style_index)
                .arg(data.is_polygon),
            defines);
        shader.bind();
        gl_engine::helpers::create_screen_quad_geometry().draw();

        const QImage render_result = b.read_colour_attachment(0);
        Framebuffer::unbind();
        CHECK(qRed(render_result.pixel(0, 0)) == 121);
        CHECK(qGreen(render_result.pixel(0, 0)) == 122);
    }
}

} // namespace

TEST_CASE("glsl tile functions")
{
    UnittestGLContext::initialise();
    std::unordered_set<Id, Id::Hasher> ids;
    {
        TileHeights h;
        h.emplace({ 0, { 0, 0 } }, { 100, 4000 });
        auto aabb_decorator = AabbDecorator::make(std::move(h));

        const auto add_tiles = [&](auto camera) {
            camera.set_viewport_size({ 1920, 1080 });
            camera.set_pixel_error_threshold(3);
            const auto all_leaves = quad_tree::onTheFlyTraverse(Id { 0, { 0, 0 } }, refineFunctor(camera, aabb_decorator, 64, 18), [&ids](const Id& v) {
                ids.insert(v);
                return v.children();
            });
        };
        add_tiles(nucleus::camera::stored_positions::grossglockner());
    }

    SECTION("hashing c++ same as glsl")
    {
        hashing_cpp_same_as_glsl({ 0, { 0, 0 } });
        hashing_cpp_same_as_glsl({ 1, { 0, 0 } });
        hashing_cpp_same_as_glsl({ 1, { 1, 1 } });
        hashing_cpp_same_as_glsl({ 14, { 2673, 12038 } });
        hashing_cpp_same_as_glsl({ 20, { 430489, 100204 } });
        for (const auto id : ids) {
            hashing_cpp_same_as_glsl(id);
        }
    }

    SECTION("packing c++ same as glsl")
    {
        packing_cpp_same_as_glsl({ 0, { 0, 0 } });
        packing_cpp_same_as_glsl({ 1, { 0, 0 } });
        packing_cpp_same_as_glsl({ 1, { 1, 1 } });
        packing_cpp_same_as_glsl({ 14, { 2673, 12038 } });
        packing_cpp_same_as_glsl({ 20, { 430489, 100204 } });
        for (const auto id : ids) {
            packing_cpp_same_as_glsl(id);
        }
    }

    SECTION("vectorlayer packing c++ same as glsl")
    {
        const auto data = std::vector<nucleus::vector_layer::VectorLayerData> {
            { glm::ivec2(8, 3), glm::ivec2(40, 36), glm::ivec2(28, 44), 289u, true },
            { glm::ivec2(34, 60), glm::ivec2(50, 52), glm::ivec2(1, 24), 565u, true },
            { glm::ivec2(16, 32), glm::ivec2(40, 5), glm::ivec2(36, 46), 630u, true },
            { glm::ivec2(46, 45), glm::ivec2(31, 4), glm::ivec2(21, 61), 546u, true },
            { glm::ivec2(27, 16), glm::ivec2(46, 7), glm::ivec2(44, 24), 610u, true },
            { glm::ivec2(48, 49), glm::ivec2(4, 11), glm::ivec2(25, 45), 103u, true },
            { glm::ivec2(44, 63), glm::ivec2(2, 10), glm::ivec2(30, 24), 852u, true },
            { glm::ivec2(55, 17), glm::ivec2(2, 47), glm::ivec2(21, 3), 912u, true },
            { glm::ivec2(32, 56), glm::ivec2(52, 56), glm::ivec2(48, 61), 124u, true },
            { glm::ivec2(45, 35), glm::ivec2(19, 46), glm::ivec2(55, 54), 482u, true },
            { glm::ivec2(5, 33), glm::ivec2(14, 32), glm::ivec2(51, 4), 478u, true },
            { glm::ivec2(43, 56), glm::ivec2(12, 49), glm::ivec2(16, 64), 389u, true },
            { glm::ivec2(63, 38), glm::ivec2(26, 60), glm::ivec2(37, 9), 192u, true },
            { glm::ivec2(27, 27), glm::ivec2(39, 33), glm::ivec2(63, 16), 179u, true },
            { glm::ivec2(0, 18), glm::ivec2(55, 14), glm::ivec2(30, 57), 34u, true },
            { glm::ivec2(17, 63), glm::ivec2(40, 33), glm::ivec2(12, 46), 15u, true },
            { glm::ivec2(5, 19), glm::ivec2(64, 34), glm::ivec2(53, 58), 303u, true },
            { glm::ivec2(40, 44), glm::ivec2(32, 29), glm::ivec2(38, 46), 151u, true },
            { glm::ivec2(16, 0), glm::ivec2(44, 55), glm::ivec2(48, 46), 819u, true },

            // testing lines (note b and c coordinates need to be the same)
            { glm::ivec2(5, 33), glm::ivec2(14, 32), glm::ivec2(0, 0), 478u, false },
            { glm::ivec2(43, 56), glm::ivec2(12, 49), glm::ivec2(0, 0), 389u, false },
            { glm::ivec2(63, 38), glm::ivec2(26, 60), glm::ivec2(0, 0), 192u, false },
            { glm::ivec2(27, 27), glm::ivec2(39, 33), glm::ivec2(0, 0), 179u, false },
            { glm::ivec2(0, 18), glm::ivec2(55, 14), glm::ivec2(0, 0), 34u, false },
            { glm::ivec2(17, 63), glm::ivec2(40, 33), glm::ivec2(0, 0), 15u, false },
            { glm::ivec2(5, 19), glm::ivec2(64, 34), glm::ivec2(0, 0), 303u, false },
            { glm::ivec2(40, 44), glm::ivec2(32, 29), glm::ivec2(0, 0), 151u, false },
            { glm::ivec2(16, 0), glm::ivec2(44, 55), glm::ivec2(0, 0), 819u, false },

        };

        for (const auto& d : data) {
            vectorlayer_packing_cpp_same_as_glsl(d);
        }
    }
}
