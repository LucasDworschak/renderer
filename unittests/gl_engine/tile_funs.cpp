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

// void vectorlayer_packing_cpp_same_as_glsl(const std::tuple<glm::ivec2, glm::ivec2, glm::ivec2, uint32_t>& data)
// {
//     {
//         Framebuffer b(Framebuffer::DepthFormat::None, { Framebuffer::ColourFormat::RGBA8 }, { 1, 1 });
//         b.bind();

//         ShaderProgram shader = create_debug_shader(QString(R"(
//             #include "vector_layer.glsl"

//             out lowp vec4 out_color;
//             void main() {
//                 highp uvec3 cpp_packed_data = uvec3(%1u, %2u, %3u);
//                 VectorLayerData data = VectorLayerData(ivec2(%4, %5),ivec2(%6, %7),ivec2(%8, %9), %10u);
//                 VectorLayerData unpacked_data = unpack_vectorlayer_data(cpp_packed_data);

//                 bool unpack_ok = data == unpacked_data;
//                 highp uvec3 packed_data = pack_vectorlayer_data(data);
//                 bool pack_ok = packed_data == cpp_packed_data;
//                 out_color = vec4(unpack_ok ? 121.0 / 255.0 : 9.0 / 255.0, pack_ok ? 122.0 / 255.0 : 9.0 / 255.0, 0, 1);
//             }
//         )")
//                 .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).x)
//                 .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).y)
//                 .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).z)
//                 .arg(std::get<0>(data).x)
//                 .arg(std::get<0>(data).y)
//                 .arg(std::get<1>(data).x)
//                 .arg(std::get<1>(data).y)
//                 .arg(std::get<2>(data).x)
//                 .arg(std::get<2>(data).y)
//                 .arg(std::get<3>(data)));
//         shader.bind();
//         gl_engine::helpers::create_screen_quad_geometry().draw();

//         const QImage render_result = b.read_colour_attachment(0);
//         Framebuffer::unbind();
//         CHECK(qRed(render_result.pixel(0, 0)) == 121);
//         CHECK(qGreen(render_result.pixel(0, 0)) == 122);
//     }
//     {
//         Framebuffer b(Framebuffer::DepthFormat::None, { Framebuffer::ColourFormat::RGBA8 }, { 1, 1 });
//         b.bind();

//         ShaderProgram shader = create_debug_shader(QString(R"(
//             out lowp vec4 out_color;
//             flat in lowp vec2 ok;
//             void main() {
//                 out_color = vec4(ok, 0, 1);
//             }
//             )"),
//             QString(R"(
//             #include "vector_layer.glsl"

//             out highp vec2 texcoords;
//             flat out lowp vec2 ok;
//             void main() {
//                 highp uvec3 cpp_packed_data = uvec3(%1u, %2u, %3u);
//                 VectorLayerData data = VectorLayerData(ivec2(%4, %5),ivec2(%6, %7),ivec2(%8, %9), %10u);
//                 VectorLayerData unpacked_data = unpack_vectorlayer_data(cpp_packed_data);

//                 bool unpack_ok = data == unpacked_data;
//                 highp uvec3 packed_data = pack_vectorlayer_data(data);
//                 bool pack_ok = packed_data == cpp_packed_data;
//                 ok = vec2(unpack_ok ? 121.0 / 255.0 : 9.0 / 255.0, pack_ok ? 122.0 / 255.0 : 9.0 / 255.0);

//                 vec2 vertices[3]=vec2[3](vec2(-1.0, -1.0), vec2(3.0, -1.0), vec2(-1.0, 3.0));
//                 gl_Position = vec4(vertices[gl_VertexID], 0.0, 1.0);
//                 texcoords = 0.5 * gl_Position.xy + vec2(0.5);
//             })")
//                 .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).x)
//                 .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).y)
//                 .arg(nucleus::vector_layer::details::pack_triangle_data(std::get<0>(data), std::get<1>(data), std::get<2>(data), std::get<3>(data)).z)
//                 .arg(std::get<0>(data).x)
//                 .arg(std::get<0>(data).y)
//                 .arg(std::get<1>(data).x)
//                 .arg(std::get<1>(data).y)
//                 .arg(std::get<2>(data).x)
//                 .arg(std::get<2>(data).y)
//                 .arg(std::get<3>(data)));
//         shader.bind();
//         gl_engine::helpers::create_screen_quad_geometry().draw();

//         const QImage render_result = b.read_colour_attachment(0);
//         Framebuffer::unbind();
//         CHECK(qRed(render_result.pixel(0, 0)) == 121);
//         CHECK(qGreen(render_result.pixel(0, 0)) == 122);
//     }
// }

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

    // SECTION("vectorlayer packing c++ same as glsl") // TODO uncomment again
    // {
    //     const auto data = std::vector<std::tuple<glm::ivec2, glm::ivec2, glm::ivec2, uint32_t>> {
    //         { glm::ivec2(2568, -731), glm::ivec2(-428, 3500), glm::ivec2(4535, 298), 4540 },
    //         { glm::ivec2(4116, 2789), glm::ivec2(3248, 2225), glm::ivec2(1464, 4824), 4721 },
    //         { glm::ivec2(-344, 2433), glm::ivec2(71, 3082), glm::ivec2(2572, 1766), 810 },
    //         { glm::ivec2(2707, 2513), glm::ivec2(1179, 4218), glm::ivec2(1744, 3508), 4507 },
    //         { glm::ivec2(3697, 1196), glm::ivec2(1888, 3972), glm::ivec2(4821, 247), 4881 },
    //         { glm::ivec2(811, -881), glm::ivec2(2687, 2825), glm::ivec2(854, 2548), 562 },
    //         { glm::ivec2(3388, 1904), glm::ivec2(2635, 767), glm::ivec2(4655, 719), 5140 },
    //         { glm::ivec2(-293, 4040), glm::ivec2(4769, 3429), glm::ivec2(1977, 3670), 4525 },
    //         { glm::ivec2(1034, 707), glm::ivec2(1429, 3563), glm::ivec2(71, 1919), 117 },
    //         { glm::ivec2(38, 4924), glm::ivec2(3496, 4930), glm::ivec2(2838, 2492), 177 },
    //         { glm::ivec2(1608, -356), glm::ivec2(3086, 4179), glm::ivec2(44, 751), 3165 },
    //         { glm::ivec2(1135, -1), glm::ivec2(1326, -571), glm::ivec2(3842, 4239), 1930 },
    //         { glm::ivec2(2198, 10), glm::ivec2(107, 226), glm::ivec2(2811, 2905), 4730 },
    //         { glm::ivec2(3992, -282), glm::ivec2(3357, 957), glm::ivec2(2681, 2104), 1769 },
    //         { glm::ivec2(-629, 819), glm::ivec2(-562, 4504), glm::ivec2(1408, 3934), 37 },
    //         { glm::ivec2(1029, 866), glm::ivec2(402, 1028), glm::ivec2(2452, -702), 3540 },
    //         { glm::ivec2(3772, 1429), glm::ivec2(4331, 3193), glm::ivec2(3084, 4924), 2149 },
    //         { glm::ivec2(-14, -852), glm::ivec2(159, 1886), glm::ivec2(-777, 2938), 4673 },
    //         { glm::ivec2(4790, 1825), glm::ivec2(84, 3121), glm::ivec2(4828, 3936), 1745 },
    //         { glm::ivec2(-79, 130), glm::ivec2(4261, -959), glm::ivec2(2291, 3268), 3762 },
    //         { glm::ivec2(4608, -306), glm::ivec2(4122, 1717), glm::ivec2(1955, -316), 436 },
    //         { glm::ivec2(-578, 468), glm::ivec2(-79, -569), glm::ivec2(975, 1022), 2952 },
    //         { glm::ivec2(410, 2462), glm::ivec2(-642, 2471), glm::ivec2(2871, 74), 2302 },
    //         { glm::ivec2(2180, 4306), glm::ivec2(-338, 78), glm::ivec2(3464, 2530), 3651 },
    //         { glm::ivec2(3768, 1701), glm::ivec2(1574, 1065), glm::ivec2(1603, -19), 3327 },
    //         { glm::ivec2(3525, 589), glm::ivec2(2822, 2682), glm::ivec2(1646, 547), 657 },
    //         { glm::ivec2(3803, 1514), glm::ivec2(-291, 2755), glm::ivec2(-126, 510), 3260 },
    //         { glm::ivec2(-91, -295), glm::ivec2(-974, 3803), glm::ivec2(4989, 4987), 230 },
    //         { glm::ivec2(3177, 1344), glm::ivec2(3243, 3118), glm::ivec2(2261, 2246), 5770 },
    //         { glm::ivec2(4861, 182), glm::ivec2(1017, 2182), glm::ivec2(556, 3455), 3232 },
    //         { glm::ivec2(-224, 1241), glm::ivec2(4800, 4865), glm::ivec2(3207, 621), 1205 },
    //         { glm::ivec2(2199, 1398), glm::ivec2(-436, 1025), glm::ivec2(-981, 2726), 3924 },
    //         { glm::ivec2(3889, 5), glm::ivec2(4635, 4479), glm::ivec2(91, 481), 3735 },
    //         { glm::ivec2(3193, -210), glm::ivec2(4159, 1941), glm::ivec2(3662, 3806), 2182 },
    //     };

    //     // const auto x = nucleus::vector_layer::details::pack_triangle_data(glm::ivec2(2568, -731), glm::ivec2(-428, 3500), glm::ivec2(4535, 298), 4540);
    //     // qDebug() << x.x << x.y << x.z;

    //     for (const auto& d : data) {
    //         vectorlayer_packing_cpp_same_as_glsl(d);
    //     }
    // }
}
