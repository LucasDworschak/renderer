/*****************************************************************************
* Alpine Renderer
* Copyright (C) 2022 Adam Celarek
* Copyright (C) 2023 Gerald Kimmersdorfer
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

#include "shared_config.glsl"
#include "hashing.glsl"
#include "camera_config.glsl"
#include "tile.glsl"

uniform highp usampler2DArray grid_sampler;
uniform highp usampler2D vector_map_meta_sampler;
uniform highp usampler2D vector_map_tile_id_sampler;

out highp vec2 var_uv;
out highp vec3 var_pos_cws;
out highp vec3 var_normal;
flat out highp uvec3 var_tile_id;
#if CURTAIN_DEBUG_MODE > 0
out lowp float is_curtain;
#endif
flat out lowp vec3 vertex_color;
flat out highp uvec3 grid_data;

lowp ivec2 to_dict_pixel(mediump uint hash) {
    return ivec2(int(hash & 255u), int(hash >> 8u));
}

bool find_tile(inout highp uvec3 tile_id, out lowp ivec2 dict_px, inout highp vec2 uv) {
    uvec2 missing_packed_tile_id = uvec2((-1u) & 65535u, (-1u) & 65535u);
    uint iter = 0u;
    do {
        mediump uint hash = hash_tile_id(tile_id);
        highp uvec2 wanted_packed_tile_id = pack_tile_id(tile_id);
        highp uvec2 found_packed_tile_id = texelFetch(vector_map_tile_id_sampler, to_dict_pixel(hash), 0).xy;
        while(found_packed_tile_id != wanted_packed_tile_id && found_packed_tile_id != missing_packed_tile_id) {
            hash++;
            found_packed_tile_id = texelFetch(vector_map_tile_id_sampler, to_dict_pixel(hash), 0).xy;
            if (iter++ > 50u) {
                break;
            }
        }
        if (found_packed_tile_id == wanted_packed_tile_id) {
            dict_px = to_dict_pixel(hash);
            tile_id = unpack_tile_id(wanted_packed_tile_id);
            return true;
        }
    }
    while (decrease_zoom_level_by_one(tile_id, uv));
    return false;
}

highp uvec3 u32_2_to_u16_u24_u24(highp uvec2 data){
    highp uvec3 res;

    res.x = data.x >> 16;

    res.y = (data.x & ((1u << 16) - 1u)) << 8;
    res.y = res.y | data.y >> 24;

    res.z = data.y & ((1u << 24) - 1u);

    return res;
}

void main() {
    float n_quads_per_direction;
    float quad_width;
    float quad_height;
    float altitude_correction_factor;
    var_pos_cws = camera_world_space_position(var_uv, n_quads_per_direction, quad_width, quad_height, altitude_correction_factor);

    if (conf.normal_mode == 1u) {
        var_normal = normal_by_finite_difference_method(var_uv, n_quads_per_direction, quad_width, quad_height, altitude_correction_factor);
    }

    var_tile_id = unpack_tile_id(packed_tile_id);
    gl_Position = camera.view_proj_matrix * vec4(var_pos_cws, 1);

    vertex_color = vec3(0.0);
    // switch(conf.overlay_mode) {
    //     case 2u: vertex_color = color_from_id_hash(uint(packed_tile_id.x ^ packed_tile_id.y)); break;
    //     case 3u: vertex_color = color_from_id_hash(uint(var_tile_id.z)); break;
    //     case 4u: vertex_color = color_from_id_hash(uint(gl_VertexID)); break;
    //     case 5u: vertex_color = vec3(texture(height_tex_sampler, vec3(var_uv, height_texture_layer)).rrr) / 65535.0; break;
    // }

    highp uvec3 tile_id = var_tile_id;
    highp vec2 uv = var_uv;

    lowp ivec2 dict_px;
    if (find_tile(tile_id, dict_px, uv)) {
        highp uvec3 vector_meta = u32_2_to_u16_u24_u24(texelFetch(vector_map_meta_sampler, dict_px, 0).rg);

        if(highp uint(-1) == vector_meta.x) // check for invalid data
            grid_data = uvec3(0u);
        else
        {
            grid_data = uvec3(texture(grid_sampler, vec3(uv, vector_meta.x)).r, vector_meta.yz);
        }
    }
    else {
        vertex_color = vec3(1.0, 0.0, 0.5);
    }

    // vertex_color = vec3(0.0);
}
