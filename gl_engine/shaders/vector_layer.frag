/*****************************************************************************
* AlpineMaps.org
* Copyright (C) 2022 Adam Celarek
* Copyright (C) 2023 Gerald Kimmersdorfer
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

#include "shared_config.glsl"
#include "camera_config.glsl"
#include "encoder.glsl"
#include "tile_id.glsl"

uniform highp usampler2DArray grid_sampler;
uniform highp usampler2D vector_map_meta_sampler;
uniform highp usampler2D vector_map_tile_id_sampler;

uniform highp usampler2D triangle_index_sampler;
uniform highp usampler2D triangle_data_sampler;

layout (location = 0) out lowp vec3 texout_albedo;
layout (location = 1) out highp vec4 texout_position;
layout (location = 2) out highp uvec2 texout_normal;
layout (location = 3) out lowp vec4 texout_depth;

flat in highp uvec3 var_tile_id;
in highp vec2 var_uv;
in highp vec3 var_pos_cws;
in highp vec3 var_normal;
#if CURTAIN_DEBUG_MODE > 0
in lowp float is_curtain;
#endif
flat in lowp vec3 vertex_color;

highp float calculate_falloff(highp float dist, highp float from, highp float to) {
    return clamp(1.0 - (dist - from) / (to - from), 0.0, 1.0);
}

highp vec3 normal_by_fragment_position_interpolation() {
    highp vec3 dFdxPos = dFdx(var_pos_cws);
    highp vec3 dFdyPos = dFdy(var_pos_cws);
    return normalize(cross(dFdxPos, dFdyPos));
}

// https://iquilezles.org/articles/distfunctions2d/
float circle(vec2 uv, vec2 pos, float r)
{
    return distance(uv, pos) - r;
}

float sdTriangle( in vec2 p, in vec2 p0, in vec2 p1, in vec2 p2 )
{
    vec2 e0 = p1-p0, e1 = p2-p1, e2 = p0-p2;
    vec2 v0 = p -p0, v1 = p -p1, v2 = p -p2;
    vec2 pq0 = v0 - e0*clamp( dot(v0,e0)/dot(e0,e0), 0.0, 1.0 );
    vec2 pq1 = v1 - e1*clamp( dot(v1,e1)/dot(e1,e1), 0.0, 1.0 );
    vec2 pq2 = v2 - e2*clamp( dot(v2,e2)/dot(e2,e2), 0.0, 1.0 );
    float s = sign( e0.x*e2.y - e0.y*e2.x );
    vec2 d = min(min(vec2(dot(pq0,pq0), s*(v0.x*e0.y-v0.y*e0.x)),
                     vec2(dot(pq1,pq1), s*(v1.x*e1.y-v1.y*e1.x))),
                     vec2(dot(pq2,pq2), s*(v2.x*e2.y-v2.y*e2.x)));
    return -sqrt(d.x)*sign(d.y);
}

highp uvec2 to_offset_size(highp uint combined) {
    // note: offset (x coord) is 24 bit -> we have to use highp
    return uvec2(uint(combined >> 8), uint(combined & 255u));
}

mediump ivec2 to_dict_pixel_2048(mediump uint hash) {
    return ivec2(int(hash & 2047u), int(hash >> 11u));
}

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
#if CURTAIN_DEBUG_MODE == 2
    if (is_curtain == 0.0) {
        discard;
    }
#endif







    // get grid acceleration structure data
    highp uvec3 tile_id = var_tile_id;
    highp vec2 uv = var_uv;

    lowp ivec2 dict_px;
    highp uvec3 grid_data = uvec3(0u);
    if (find_tile(tile_id, dict_px, uv)) {
        highp uvec3 vector_meta = u32_2_to_u16_u24_u24(texelFetch(vector_map_meta_sampler, dict_px, 0).rg);
        if(vector_meta.x != highp uint(-1)) // check for valid data
        {
            grid_data = uvec3(texture(grid_sampler, vec3(uv, vector_meta.x)).r, vector_meta.yz);
        }
    }

    texout_albedo = vec3(0.0, 0.0, 0.0);

    // using the grid data we now want to traverse all triangles referenced in grid cell and draw them.
    if(grid_data.x != uint(-1)) // only if we have data here
    {
        // calculate the offset and size of the triangle index list
        // the offset in grid_data.x is local to the current tile, the value in grid_data.y is global to all tiles
        highp uvec2 offset_size = to_offset_size(grid_data.x) + uvec2(grid_data.y, 0);

        lowp vec3 raw_grid = vec3(float(offset_size.y),0,0);// DEBUG

        float alpha = 0.0;

        for(highp uint i = offset_size.x; i < offset_size.x + offset_size.y; i++)
        {
            highp uint triangle_index = texelFetch(triangle_index_sampler, to_dict_pixel_2048(i), 0).r;

            highp vec2 v0;
            v0.x = uintBitsToFloat(texelFetch(triangle_data_sampler, to_dict_pixel_2048(triangle_index+0u), 0).r) / 64.0f;
            v0.y = uintBitsToFloat(texelFetch(triangle_data_sampler, to_dict_pixel_2048(triangle_index+1u), 0).r) / 64.0f;

            highp vec2 v1;
            v1.x = uintBitsToFloat(texelFetch(triangle_data_sampler, to_dict_pixel_2048(triangle_index+2u), 0).r) / 64.0f;
            v1.y = uintBitsToFloat(texelFetch(triangle_data_sampler, to_dict_pixel_2048(triangle_index+3u), 0).r) / 64.0f;

            highp vec2 v2;
            v2.x = uintBitsToFloat(texelFetch(triangle_data_sampler, to_dict_pixel_2048(triangle_index+4u), 0).r) / 64.0f;
            v2.y = uintBitsToFloat(texelFetch(triangle_data_sampler, to_dict_pixel_2048(triangle_index+5u), 0).r) / 64.0f;

            highp uint style_index = texelFetch(triangle_data_sampler, to_dict_pixel_2048(triangle_index+6u), 0).r;

            highp float c1 = 1.0 - step(0.0, sdTriangle(var_uv, v0, v1, v2) - 0.078125f);
            vec3 triangle_out = vec3(0.0f, c1, 0.0f);
            alpha += c1;

            texout_albedo = mix(raw_grid, triangle_out, 0.5);

            if(alpha > 1.0)
                break; // early exit if alpha is 1;
        }
    }

    // Write Position (and distance) in gbuffer
    highp float dist = length(var_pos_cws);
    texout_position = vec4(var_pos_cws, dist);

    // Write and encode normal in gbuffer
    highp vec3 normal = vec3(0.0);
    if (conf.normal_mode == 0u) normal = normal_by_fragment_position_interpolation();
    else normal = var_normal;
    texout_normal = octNormalEncode2u16(normal);

    // Write and encode distance for readback
    texout_depth = vec4(depthWSEncode2n8(dist), 0.0, 0.0);

    // HANDLE OVERLAYS (and mix it with the albedo color) THAT CAN JUST BE DONE IN THIS STAGE
    // (because of DATA thats not forwarded)
    // NOTE: Performancewise its generally better to handle overlays in the compose step! (screenspace effect)
    if (conf.overlay_mode > 0u && conf.overlay_mode < 100u) {
        lowp vec3 overlay_color = vec3(0.0);
        switch(conf.overlay_mode) {
            case 1u: overlay_color = normal * 0.5 + 0.5; break;
            default: overlay_color = vertex_color;
        }
        texout_albedo = mix(texout_albedo, overlay_color, conf.overlay_strength);
    }

#if CURTAIN_DEBUG_MODE == 1
    if (is_curtain > 0.0) {
        texout_albedo = vec3(1.0, 0.0, 0.0);
        return;
    }
#endif

}
