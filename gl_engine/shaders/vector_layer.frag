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
#include "vector_layer.glsl"

#include "hashing.glsl" // DEBUG

uniform highp usampler2DArray acceleration_grid_sampler;
uniform highp usampler2D array_index_sampler;
uniform highp usampler2D vector_map_tile_id_sampler;

uniform highp usampler2D styles_sampler;

uniform highp usampler2DArray index_buffer_sampler_0;
uniform highp usampler2DArray index_buffer_sampler_1;
uniform highp usampler2DArray index_buffer_sampler_2;
uniform highp usampler2DArray index_buffer_sampler_3;
uniform highp usampler2DArray vertex_buffer_sampler_0;
uniform highp usampler2DArray vertex_buffer_sampler_1;
uniform highp usampler2DArray vertex_buffer_sampler_2;
uniform highp usampler2DArray vertex_buffer_sampler_3;

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

struct Style_Data
{
    lowp vec4 fill_color;
    lowp vec4 outline_color;
    lowp float outline_width;
    lowp vec2 outline_dash;
};

struct Layer_Style
{
    highp uint last_style;
    lowp float layer_alpha;
    Style_Data current_layer_style;
};


///////////////////////////////////////////////
// CONSTANTS

const lowp int grid_size = 64;
const lowp uint data_size = 1u; // how many texels are needed to store the data for one triangle

const lowp uint sampler_offset = 16u - 2u; // used to calculate how many bits are used to determine the sampler index, and how many are used for the layer
const highp uint layer_mask = ((1u << sampler_offset) - 1u);

const lowp float style_precision = 100;

///////////////////////////////////////////////


highp float calculate_falloff(highp float dist, highp float from, highp float to) {
    return clamp(1.0 - (dist - from) / (to - from), 0.0, 1.0);
}

highp vec3 normal_by_fragment_position_interpolation() {
    highp vec3 dFdxPos = dFdx(var_pos_cws);
    highp vec3 dFdyPos = dFdy(var_pos_cws);
    return normalize(cross(dFdxPos, dFdyPos));
}

// https://iquilezles.org/articles/distfunctions2d/
float circle(vec2 uv, vec2 pos, float r) // DEBUG
{
    return distance(uv, pos) - r;
}

float sdLine( in vec2 p, in vec2 a, in vec2 b )
{
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
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

mediump ivec2 to_dict_pixel_64(mediump uint hash) {
    return ivec2(int(hash & 63u), int(hash >> 6u));
}
mediump ivec2 to_dict_pixel_128(mediump uint hash) {
    return ivec2(int(hash & 127u), int(hash >> 7u));
}
mediump ivec2 to_dict_pixel_256(mediump uint hash) {
    return ivec2(int(hash & 255u), int(hash >> 8u));
}
mediump ivec2 to_dict_pixel_512(mediump uint hash) {
    return ivec2(int(hash & 511u), int(hash >> 9u));
}
mediump ivec2 to_dict_pixel_1024(mediump uint hash) {
    return ivec2(int(hash & 1023u), int(hash >> 10u));
}
mediump ivec2 to_dict_pixel_2048(mediump uint hash) {
    return ivec2(int(hash & 2047u), int(hash >> 11u));
}

lowp ivec2 to_dict_pixel(mediump uint hash) {
    return ivec2(int(hash & 255u), int(hash >> 8u));
}

highp uint index_sample(lowp uint sampler_index, highp uint pixel_index, highp uint texture_layer)
{
    if(sampler_index == 0u)
    {
        return texelFetch(index_buffer_sampler_0, ivec3(to_dict_pixel_64(pixel_index), texture_layer), 0).r * data_size;
    }
    else if(sampler_index == 1u)
    {
        return texelFetch(index_buffer_sampler_1, ivec3(to_dict_pixel_128(pixel_index), texture_layer), 0).r * data_size;
    }
    else if(sampler_index == 2u)
    {
        return texelFetch(index_buffer_sampler_2, ivec3(to_dict_pixel_256(pixel_index), texture_layer), 0).r * data_size;
    }
    else
    {
        return texelFetch(index_buffer_sampler_3, ivec3(to_dict_pixel_512(pixel_index), texture_layer), 0).r * data_size;
    }
}

VectorLayerData vertex_sample(lowp uint sampler_index, highp uint index, highp uint texture_layer)
{

    if(sampler_index == 0u)
    {
        return unpack_vectorlayer_data(texelFetch(vertex_buffer_sampler_0, ivec3(to_dict_pixel_64(index), texture_layer), 0).rgb);
    }
    else if(sampler_index == 1u)
    {
        return unpack_vectorlayer_data(texelFetch(vertex_buffer_sampler_1, ivec3(to_dict_pixel_128(index), texture_layer), 0).rgb);
    }
    else if(sampler_index == 2u)
    {
        return unpack_vectorlayer_data(texelFetch(vertex_buffer_sampler_2, ivec3(to_dict_pixel_256(index), texture_layer), 0).rgb);
    }
    else
    {
        return unpack_vectorlayer_data(texelFetch(vertex_buffer_sampler_3, ivec3(to_dict_pixel_512(index), texture_layer), 0).rgb);
    }
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

Style_Data parse_style(highp uint style_index) {
    Style_Data style;

    highp uvec4 style_data = texelFetch(styles_sampler, ivec2(to_dict_pixel_64(style_index)), 0);
    style.fill_color = vec4((style_data.r & 4278190080u) >> 24, (style_data.r & 16711680u) >> 16, (style_data.r & 65280u) >> 8, style_data.r & 255u) / vec4(255.0f);
    style.outline_color = vec4((style_data.g & 4278190080u) >> 24, (style_data.g & 16711680u) >> 16, (style_data.g & 65280u) >> 8, style_data.g & 255u) / vec4(255.0f);

    style.outline_width = float(style_data.b) / style_precision;
    // lowp vec2 outline_dash; // TODO

    return style;
}


bool prepare_layer_style(highp uint style_index, inout Layer_Style layer_style, inout lowp vec3 pixel_color)
{
    // we need to make sure that a layer with < 1 opacity does only fill the correct amount of opacity
    // we therefore fill colors per layerstyle
    if(layer_style.last_style == style_index)
    {
        if(layer_style.layer_alpha >= 1.0)
        {
            // we already fully filled the current layer -> go to the next layer
            return true;
        }
    }
    else
    {
        // we encountered a new layer

        // mix the previous layer color information with output
        pixel_color = mix(pixel_color, layer_style.current_layer_style.fill_color.rgb, layer_style.layer_alpha * layer_style.current_layer_style.fill_color.a);

        // get and store new style info
        layer_style.last_style = style_index;
        layer_style.current_layer_style = parse_style(style_index);
        // how much alpha per layer we accumulate
        layer_style.layer_alpha = 0.0;
    }

    return false;

}



void main() {
#if CURTAIN_DEBUG_MODE == 2
    if (is_curtain == 0.0) {
        discard;
    }
#endif


    // highp float depth = depthWSEncode2n8(length(var_pos_cws));
    highp float depth = length(var_pos_cws);







    // get grid acceleration structure data
    highp uvec3 tile_id = var_tile_id;
    highp vec2 uv = var_uv;

    lowp ivec2 dict_px;
    highp uvec2 offset_size = uvec2(0u);
    texout_albedo = vec3(242.0f/255.0f, 239.0f/255.0f, 233.0f/255.0f);

    if (find_tile(tile_id, dict_px, uv)) {


        highp uvec2 texture_layer = texelFetch(array_index_sampler, dict_px, 0).xy;
        if(texture_layer.x != highp uint(-1) && texture_layer.y != highp uint(-1)) // check for valid data
        {
            // acceleration_grid_sampler contains the offset and the number of triangles of the current grid cell
            offset_size = to_offset_size(texelFetch(acceleration_grid_sampler, ivec3(uv*vec2(grid_size,grid_size), texture_layer.x & layer_mask),0).r);

            // using the grid data we now want to traverse all triangles referenced in grid cell and draw them.
            if(offset_size.y != uint(0)) // only if we have data here
            {
                // lowp vec3 raw_grid = vec3(float(offset_size.y),0,0);// DEBUG
                lowp vec3 raw_grid = vec3(1,0,0);// DEBUG
                ivec2 grid_cell = ivec2(uv*vec2(grid_size,grid_size)); // DEBUG
                lowp vec3 cells = color_from_id_hash(uint(grid_cell.x ^ grid_cell.y)); // DEBUG
                vec3 triangle_lines_out = vec3(0.0f, 0.0, 0.0f); // DEBUG

                lowp vec3 pixel_color = vec3(0.0f, 0.0, 0.0f);
                float pixel_alpha = 0.0;

                // get the buffer index and extract the correct texture_layer.y
                lowp uint sampler_buffer_index = (texture_layer.y & ((highp uint(-1u) << sampler_offset))) >> sampler_offset;
                texture_layer.y = texture_layer.y & layer_mask;

                lowp vec3 layer_debug = vec3(0,0,0);// DEBUG -> which cascade is being used
                {
                    if(sampler_buffer_index == 0u)
                        layer_debug = vec3(0,1,0); // green
                    else if(sampler_buffer_index == 1u)
                        layer_debug = vec3(1,1,0); // yellow
                    else if(sampler_buffer_index == 2u)
                        layer_debug = vec3(1,0.5,0); // orange
                    else if(sampler_buffer_index == 3u)
                        layer_debug = vec3(1,0,0); // red
                    else
                        layer_debug = vec3(1,0,1); // purple -> should never happen -> unrecognized index
                }
                lowp vec3 cell_debug = vec3(0,0,0);// DEBUG -> how many triangles per cell
                {
                    if(offset_size.y < 32u)
                        cell_debug = vec3(0,1,0); // green
                    else if(offset_size.y  < 64u)
                        cell_debug = vec3(1,1,0); // yellow
                    else if(offset_size.y  < 128u)
                        cell_debug = vec3(1,0.5,0); // orange
                    else if(offset_size.y  < 255u)
                        cell_debug = vec3(1,0,0); // red
                    else
                        cell_debug = vec3(1,0,1); // purple -> should never happen -> unrecognized index
                }

                Layer_Style layer_style;

                layer_style.last_style = -1u;
                layer_style.current_layer_style = Style_Data(vec4(0.0), vec4(0.0), 0.0, vec2(0.0));
                // layer_style.current_layer_style.fill_color = vec4(0.0);
                layer_style.layer_alpha = 0.0;

                highp float geometry_influence = 0.0;

                for(highp uint i = offset_size.x; i < offset_size.x + offset_size.y; i++)
                // for(highp uint i = offset_size.x+ offset_size.y; i --> offset_size.x ; ) // reverse traversal
                {                    
                    highp uint index = index_sample(sampler_buffer_index, i, texture_layer.y);
                    bool is_polygon = (index & 1u) == 1u;
                    index = index >> 1;

                    float d = 0.0;
                    highp uint style_index = -1u;

                    if(is_polygon)
                    {
                        VectorLayerData triangle_data = vertex_sample(sampler_buffer_index, index, texture_layer.y);

                        highp vec2 v0 = triangle_data.a / vec2(tile_extent);
                        highp vec2 v1 = triangle_data.b / vec2(tile_extent);
                        highp vec2 v2 = triangle_data.c / vec2(tile_extent);

                        float thickness = 0.0;
                        d = sdTriangle(uv, v0, v1, v2) - thickness;
                        // d = 100.0;

                        geometry_influence = 1.0 - step(0.0, d);

                        // polygon does not influence pixel at all -> we do not draw it
                        if(geometry_influence <= 0.0)
                            continue;

                        // calling it here prevents getting the layerstyle if we do not need it yet
                        bool check_next_geometry = prepare_layer_style(triangle_data.style_index, layer_style, pixel_color);
                        if(check_next_geometry)
                            continue;

                    }
                    else
                    {
                        VectorLayerData line_data = vertex_sample(sampler_buffer_index, index, texture_layer.y);

                        highp vec2 v0 = line_data.a / vec2(tile_extent);
                        highp vec2 v1 = line_data.b / vec2(tile_extent);

                        // needs to be applied here to get the thickness of the line
                        bool check_next_geometry = prepare_layer_style(line_data.style_index, layer_style, pixel_color);
                        if(check_next_geometry)
                            continue;


                        float thickness = layer_style.current_layer_style.outline_width / tile_extent;
                        d = sdLine(uv, v0, v1) - thickness;

                        // d = 100.0;
                        // style_index = 36u;


                        geometry_influence = 1.0 - step(0.0, d);

                        // polygon does not influence pixel at all -> we do not draw it
                        if(geometry_influence <= 0.0)
                            continue;


                    }





                    { // DEBUG -> triangle lines
                        highp float t_line = 0.0;
                        if(d <= 0.0001 && d >= -0.0001)
                            t_line = 0.2;
                        triangle_lines_out += vec3(0.0f, t_line, 0.0f);
                    }







                    pixel_alpha += geometry_influence * layer_style.current_layer_style.fill_color.a;
                    if(pixel_alpha > 1.0)
                        geometry_influence = pixel_alpha - 1.0;

                    layer_style.layer_alpha += geometry_influence;

                    if(pixel_alpha >= 1.0)
                        break; // early exit if alpha is >= 1;

                }


                // mix the last layer we parsed
                pixel_color = mix(pixel_color, layer_style.current_layer_style.fill_color.rgb, layer_style.layer_alpha * layer_style.current_layer_style.fill_color.a);

                // mix polygon color with background
                texout_albedo = mix(texout_albedo, pixel_color, pixel_alpha);

                // texout_albedo = pixel_color;

                // texout_albedo = mix(cells, texout_albedo, 0.9);// DEBUG
                // texout_albedo = mix(texout_albedo, triangle_lines_out, 0.5);// DEBUG

                // vec3 grid_start = color_from_id_hash(offset_size.x);
                // vec3 grid_start = color_from_id_hash(offset_size.y);

                // texout_albedo = grid_start;

                // texout_albedo = mix(grid_start, pixel_color, 0.5);// DEBUG
                // texout_albedo = mix(raw_grid, pixel_color, 0.5);// DEBUG
                // texout_albedo = raw_grid;// DEBUG
                // texout_albedo = vec3(pixel_alpha);// DEBUG
                // texout_albedo = triangle_lines_out;// DEBUG
                // texout_albedo = layer_debug;// DEBUG
                // texout_albedo = cell_debug;// DEBUG
                 // texout_albedo = mix(texout_albedo, color_from_id_hash(texture_layer.y), 0.2);

            }
            else
            {
                // we aren't processing anything
                // texout_albedo = vec3(0.0, 0.0, 0.0);// DEBUG
            }

            // float min_size = 50;
            // float max_size = 160;
            // texout_albedo = vec3(((clamp(offset_size.y, min_size, max_size) - min_size) / (max_size-min_size)) / 2.0 + (0.5 * ((offset_size.y > min_size)? 1.0 : 0.0)), 0.0 ,0.0);

        }
    }


    // DEBUG tile bounds
    // if(uv.x <= 0.01 || uv.y <= 0.01)
    // {
    //     texout_albedo = vec3(1.0, 0.0, 0.0);
    // }






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
