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
#include <gl_engine/helpers.h>
#include <nucleus/camera/PositionStorage.h>
#include <nucleus/srs.h>
#include <nucleus/tile/utils.h>
#include <nucleus/utils/ColourTexture.h>
#include <nucleus/vector_layer/Preprocessor.h>
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
ShaderProgram create_debug_shader(const QString& fragment_shader, const QString& vertex_shader = R"(
out highp vec2 texcoords;
void main() {
    vec2 vertices[3]=vec2[3](vec2(-1.0, -1.0), vec2(3.0, -1.0), vec2(-1.0, 3.0));
    gl_Position = vec4(vertices[gl_VertexID], 0.0, 1.0);
    texcoords = 0.5 * gl_Position.xy + vec2(0.5);
})")
{
    ShaderProgram tmp(vertex_shader, fragment_shader, gl_engine::ShaderCodeSource::PLAINTEXT);
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

void vectorlayer_packing_cpp_same_as_glsl(const std::tuple<glm::vec2, glm::vec2, glm::vec2, uint32_t>& data)
{
    {
        Framebuffer b(Framebuffer::DepthFormat::None, { Framebuffer::ColourFormat::RGBA8 }, { 1, 1 });
        b.bind();

        ShaderProgram shader = create_debug_shader(QString(R"(
            #include "vector_layer.glsl"

            out lowp vec4 out_color;
            void main() {
                highp uvec3 cpp_packed_data = uvec3(%1u, %2u, %3u);
                VectorLayerData data = VectorLayerData(ivec2(%4u, %5u),ivec2(%6u, %7u),ivec2(%8u, %9u), %10u);
                VectorLayerData unpacked_data = unpack_vectorlayer_data(cpp_packed_data);

                bool unpack_ok = data == unpacked_data;
                highp uvec3 packed_data = pack_vectorlayer_data(data);
                bool pack_ok = packed_data == cpp_packed_data;
                out_color = vec4(unpack_ok ? 121.0 / 255.0 : 9.0 / 255.0, pack_ok ? 122.0 / 255.0 : 9.0 / 255.0, 0, 1);
            }
        )")
                .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).x)
                .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).y)
                .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).z)
                .arg(std::get<0>(data).x)
                .arg(std::get<0>(data).y)
                .arg(std::get<1>(data).x)
                .arg(std::get<1>(data).y)
                .arg(std::get<2>(data).x)
                .arg(std::get<2>(data).y)
                .arg(std::get<3>(data)));
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
            #include "vector_layer.glsl"

            out highp vec2 texcoords;
            flat out lowp vec2 ok;
            void main() {
                highp uvec3 cpp_packed_data = uvec3(%1u, %2u, %3u);
                VectorLayerData data = VectorLayerData(ivec2(%4u, %5u),ivec2(%6u, %7u),ivec2(%8u, %9u), %10u);
                VectorLayerData unpacked_data = unpack_vectorlayer_data(cpp_packed_data);

                bool unpack_ok = data == unpacked_data;
                highp uvec3 packed_data = pack_vectorlayer_data(data);
                bool pack_ok = packed_data == cpp_packed_data;
                ok = vec2(unpack_ok ? 121.0 / 255.0 : 9.0 / 255.0, pack_ok ? 122.0 / 255.0 : 9.0 / 255.0);

                vec2 vertices[3]=vec2[3](vec2(-1.0, -1.0), vec2(3.0, -1.0), vec2(-1.0, 3.0));
                gl_Position = vec4(vertices[gl_VertexID], 0.0, 1.0);
                texcoords = 0.5 * gl_Position.xy + vec2(0.5);
            })")
                .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).x)
                .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).y)
                .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).z)
                .arg(std::get<0>(data).x)
                .arg(std::get<0>(data).y)
                .arg(std::get<1>(data).x)
                .arg(std::get<1>(data).y)
                .arg(std::get<2>(data).x)
                .arg(std::get<2>(data).y)
                .arg(std::get<3>(data)));
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
            const auto all_leaves = quad_tree::onTheFlyTraverse(Id { 0, { 0, 0 } }, refineFunctor(camera, aabb_decorator, 3, 64), [&ids](const Id& v) {
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
        const auto data = std::vector<std::tuple<glm::vec2, glm::vec2, glm::vec2, uint32_t>> { { glm::vec2(3637, 2602), glm::vec2(1816, 1568), glm::vec2(1048, 15), 1497 },
            { glm::vec2(2423, 2882), glm::vec2(2432, -3051), glm::vec2(-1249, 496), 171 },
            { glm::vec2(-1863, 1487), glm::vec2(3195, 2862), glm::vec2(751, -1107), 1733 },
            { glm::vec2(1976, -3478), glm::vec2(-3281, 2244), glm::vec2(906, -2299), 3564 },
            { glm::vec2(-3279, 772), glm::vec2(3118, -1859), glm::vec2(-140, 1085), 2519 },
            { glm::vec2(844, 2525), glm::vec2(4008, 616), glm::vec2(1371, 2862), 3525 },
            { glm::vec2(-770, 3363), glm::vec2(2393, -3696), glm::vec2(3884, -3747), 2627 },
            { glm::vec2(1787, -3877), glm::vec2(-4002, 357), glm::vec2(-1805, 1013), 2277 },
            { glm::vec2(3565, 4005), glm::vec2(-3499, 158), glm::vec2(3524, 2306), 1203 },
            { glm::vec2(1849, 3503), glm::vec2(1141, 321), glm::vec2(3914, -617), 1291 },
            { glm::vec2(-796, -3154), glm::vec2(4005, 2381), glm::vec2(-2881, 1556), 1104 },
            { glm::vec2(3248, 3527), glm::vec2(-1668, -2124), glm::vec2(413, 3317), 1538 },
            { glm::vec2(2044, -3993), glm::vec2(216, 363), glm::vec2(3365, -1947), 828 },
            { glm::vec2(-315, 3953), glm::vec2(162, 2457), glm::vec2(-523, 2958), 3046 },
            { glm::vec2(2101, 1822), glm::vec2(520, -348), glm::vec2(3279, -1007), 853 },
            { glm::vec2(507, 3648), glm::vec2(-143, 3612), glm::vec2(-3560, 853), 2070 },
            { glm::vec2(965, -1242), glm::vec2(860, 3314), glm::vec2(1488, 1195), 684 },
            { glm::vec2(-331, 2837), glm::vec2(1067, 2650), glm::vec2(255, -1344), 2669 },
            { glm::vec2(3788, 1900), glm::vec2(764, -1566), glm::vec2(1482, 129), 1863 },
            { glm::vec2(1545, -1604), glm::vec2(-2837, 4094), glm::vec2(-4060, 3153), 1554 },
            { glm::vec2(-3839, 754), glm::vec2(2806, 254), glm::vec2(709, -3693), 2456 },
            { glm::vec2(2545, 817), glm::vec2(597, -3773), glm::vec2(456, 219), 312 },
            { glm::vec2(3504, 1471), glm::vec2(-3802, 3467), glm::vec2(-2926, 1818), 1352 },
            { glm::vec2(-1158, -3941), glm::vec2(2641, -608), glm::vec2(2567, 1767), 1402 },
            { glm::vec2(-3772, -3133), glm::vec2(-898, -3267), glm::vec2(-955, -3670), 3917 },
            { glm::vec2(2778, 2997), glm::vec2(1602, -3743), glm::vec2(3543, 4086), 3988 } };

        for (const auto d : data) {
            vectorlayer_packing_cpp_same_as_glsl(d);
        }
    }
}
