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

#line 20031

uniform lowp sampler2DArray texture_sampler;
uniform highp usampler2D instanced_texture_array_index_sampler;
uniform highp usampler2D instanced_texture_zoom_sampler;

uniform highp usampler2DArray acceleration_grid_sampler;
uniform highp usampler2D instanced_texture_array_index_sampler_vector;
uniform highp usampler2D instanced_texture_zoom_sampler_vector;

uniform highp usampler2D styles_sampler;

uniform highp usampler2DArray geometry_buffer_sampler_0;
uniform highp usampler2DArray geometry_buffer_sampler_1;
uniform highp usampler2DArray geometry_buffer_sampler_2;
uniform highp usampler2DArray geometry_buffer_sampler_3;

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
flat in highp uint instance_id;

struct LayerStyle {
    highp uint index;
    lowp vec4 color;
    lowp float line_width;
    lowp vec2 dash_info;
    lowp int line_caps;
};

const lowp int n_aa_samples_row_cols = 4;
const lowp int n_aa_samples = n_aa_samples_row_cols*n_aa_samples_row_cols;
const lowp float aa_sample_dist = 1.0; // 0.5 goes from -0.25 to +0.25 of the current uv coordinate

// SDF_MODE 0: naive multisample antialiasing (golden -> but worse performnace)
// SDF_MODE 1: derivative multisample antialiasing
#define SDF_MODE 1

const mediump float division_by_n_samples = 1.0 / float(n_aa_samples);
const lowp float style_precision_mult = 1.0 / float(style_precision);
const mediump float inv_tile_extent = 1.0 / tile_extent;




///////////////////////////////////////////////
// CPP CONFIG CONSTANTS

// defined in c++ code
// const highp float tile_extent
// const highp vec2 grid_size
// const lowp float style_precision
// const lowp int style_bits
// const lowp float max_zoom
// const lowp int zoom_blend_steps

///////////////////////////////////////////////

const highp vec2 cell_size_uv = vec2(1.0) / vec2(grid_size);
const highp vec2 aa_cell_overlap_uv = vec2(aa_border) * vec2(cell_size_uv);

const highp uint layer_mask = ((1u << sampler_offset) - 1u);
const highp uint bit_mask_ones = -1u;

const highp float full_threshold = 0.99;



highp float calculate_falloff(highp float dist, highp float from, highp float to) {
    return clamp(1.0 - (dist - from) / (to - from), 0.0, 1.0);
}

highp vec3 normal_by_fragment_position_interpolation() {
    highp vec3 dFdxPos = dFdx(var_pos_cws);
    highp vec3 dFdyPos = dFdy(var_pos_cws);
    return normalize(cross(dFdxPos, dFdyPos));
}

#if SDF_MODE == 0
void calculate_sample_positions(out highp vec2 aa_sample_positions[n_aa_samples], highp vec2 uv, highp ivec2 grid_lookup)
{
    if(n_aa_samples == 1)
    {
        aa_sample_positions[0] = uv;
        return;
    }

    highp vec2 min_cell = vec2(grid_lookup) * cell_size_uv;
    highp vec2 max_cell = min_cell + cell_size_uv;
    min_cell -= vec2(aa_cell_overlap_uv);
    max_cell += vec2(aa_cell_overlap_uv);


    highp vec2 grad_u = vec2(dFdx(uv.x), dFdy(uv.x));
    highp vec2 grad_v = vec2(dFdx(uv.y), dFdy(uv.y));

    highp float aa_sample_dist_increments = aa_sample_dist / float(n_aa_samples_row_cols-1.0);

    for (int x = 0; x < n_aa_samples_row_cols; ++x) {
        for (int y = 0; y < n_aa_samples_row_cols; ++y) {
            lowp int index = x*n_aa_samples_row_cols + y;

            highp vec2 aa_sample_multiplier = vec2(-aa_sample_dist / 2.0) + vec2(x,y) * vec2(aa_sample_dist_increments);

            aa_sample_positions[index].x = uv.x + dot(grad_u, aa_sample_multiplier);
            aa_sample_positions[index].y = uv.y + dot(grad_v, aa_sample_multiplier);

            // clip to current cell
            aa_sample_positions[index] = min(max_cell, max(min_cell, aa_sample_positions[index]));

        }
    }

}
#endif
#if SDF_MODE == 1
void calculate_sample_multipliers(out highp vec2 aa_sample_multipliers[n_aa_samples])
{
    if(n_aa_samples == 1)
    {
        aa_sample_multipliers[0] = vec2(0.0);
        return;
    }

    highp float aa_sample_dist_increments = aa_sample_dist / float(n_aa_samples_row_cols-1);

    for (int x = 0; x < n_aa_samples_row_cols; ++x) {
        for (int y = 0; y < n_aa_samples_row_cols; ++y) {
            lowp int index = x*n_aa_samples_row_cols + y;
            aa_sample_multipliers[index] = vec2(-aa_sample_dist / 2.0) + vec2(x,y) * vec2(aa_sample_dist_increments);
        }
    }

}
#endif


#if SDF_MODE == 0
struct SDFData
{
    VectorLayerData vertices;

    highp vec2 e0;
    highp vec2 e1;
    highp vec2 e2;

    highp float dot_e0;
    highp float dot_e1;
    highp float dot_e2;
};

SDFData prepare_sd_Line_Triangle(VectorLayerData geom_data)
{
    SDFData data;

    data.vertices = geom_data;

    data.e0 = data.vertices.b-data.vertices.a;
    data.dot_e0 = dot(data.e0,data.e0);

    if(geom_data.is_polygon)
    {
        data.e1 = data.vertices.c-data.vertices.b;
        data.e2 = data.vertices.a-data.vertices.c;
        data.dot_e1 = dot(data.e1,data.e1);
        data.dot_e2 = dot(data.e2,data.e2);
    }

    return data;
}

// https://iquilezles.org/articles/distfunctions2d/
highp float sd_Line_Triangle( in highp vec2 uv, SDFData data, bool triangle )
{
    highp vec2 v0 = uv-data.vertices.a;
    highp vec2 pq0 = v0 - data.e0*clamp( dot(v0,data.e0)/data.dot_e0, 0.0, 1.0 );

    highp float poly_sign = 1.0;
    highp float result = 1.0;

    if(triangle)
    {
        highp vec2 v1 = uv  - data.vertices.b;
        highp vec2 v2 = uv  - data.vertices.c;
        highp vec2 pq1 = v1 - data.e1*clamp( dot(v1,data.e1)/data.dot_e1, 0.0, 1.0 );
        highp vec2 pq2 = v2 - data.e2*clamp( dot(v2,data.e2)/data.dot_e2, 0.0, 1.0 );
        highp float s = sign( data.e0.x*data.e2.y - data.e0.y*data.e2.x );
        highp vec2 d0 = vec2(dot(pq0,pq0), s*(v0.x*data.e0.y-v0.y*data.e0.x));
        highp vec2 d1 = vec2(dot(pq1,pq1), s*(v1.x*data.e1.y-v1.y*data.e1.x));
        highp vec2 d = min(d0,d1);
        highp vec2 d2 = vec2(dot(pq2,pq2), s*(v2.x*data.e2.y-v2.y*data.e2.x));
        d = min(d,d2);

        poly_sign = -sign(d.y);
        result = d.x;

        // return -sqrt(d.x)*sign(d.y); // TODO sqrt and length is basically the same instruction -> get them out of the if
    }
    else{
        poly_sign = 1.0;
        result = dot(pq0,pq0);

        // return sqrt(dot(pq0,pq0));
    }

    return sqrt(result)*poly_sign;

}
#endif
#if SDF_MODE == 1

highp float grad_clamp_fun(highp float v, highp float min, highp float max, highp float incoming_grad)
{
    // return incoming_grad * // TODO convert to * instead of if
    if (v > min && v < max)
        return incoming_grad;
    return 0.0;
}


highp vec3 sdf_with_grad(VectorLayerData data, highp vec2 uv, highp float incoming_grad)
{
    highp vec2 e0 = data.b - data.a;
    highp vec2 v0 = uv - data.a;
    highp float dot0 = dot(v0, e0);
    highp float one_over_dot0 = 1.0 / dot(e0, e0);
    highp float div0 = dot0 * one_over_dot0;
    highp vec2 pq0 = v0 - e0 * clamp(div0, 0.0, 1.0);
    highp float dot_pq0_pq0 = dot(pq0, pq0);

    highp float poly_sign = 1.0;
    highp float distance_sq = 1.0;

    highp vec2 grad_uv = vec2(0.0);
    highp vec2 grad_pq0 = vec2(0.0);
    if (data.is_polygon) {
        highp vec2 e1 = data.c - data.b;
        highp vec2 e2 = data.a - data.c;
        highp vec2 v1 = uv - data.b;
        highp vec2 v2 = uv - data.c;
        highp float dot1 = dot(v1, e1);
        highp float dot2 = dot(v2, e2);
        highp float one_over_dot1 = 1.0 / dot(e1, e1);
        highp float one_over_dot2 = 1.0 / dot(e2, e2);
        highp float div1 = dot1 * one_over_dot1;
        highp float div2 = dot2 * one_over_dot2;
        highp float clamp1 = clamp(div1, 0.0, 1.0);
        highp float clamp2 = clamp(div2, 0.0, 1.0);
        highp vec2 pq1 = v1 - e1 * clamp1;
        highp vec2 pq2 = v2 - e2 * clamp2;
        highp float s = sign(e0.x * e2.y - e0.y * e2.x);
        highp vec2 d0 = vec2(dot_pq0_pq0, s * (v0.x * e0.y - v0.y * e0.x));
        highp vec2 d1 = vec2(dot(pq1, pq1), s * (v1.x * e1.y - v1.y * e1.x));
        highp vec2 d2 = vec2(dot(pq2, pq2), s * (v2.x * e2.y - v2.y * e2.x));
        highp vec2 d = min(min(d0, d1), d2);

        poly_sign = -sign(d.y);
        distance_sq = d.x;

        // gradient computation
        highp float grad_d0_x = 0.0;
        highp float grad_d1_x = 0.0;
        highp float grad_d2_x = 0.0;

        if (d0.x <= d1.x && d0.x <= d2.x) {
            grad_d0_x = incoming_grad;
        } else if (d1.x < d0.x && d1.x <= d2.x) {
            grad_d1_x = incoming_grad;
        } else {
            grad_d2_x = incoming_grad;
        }
        grad_pq0 = 2.0 * pq0 * grad_d0_x;
        highp vec2 grad_pq1 = 2.0 * pq1 * grad_d1_x;
        highp vec2 grad_pq2 = 2.0 * pq2 * grad_d2_x;

        highp vec2 grad_v1 = grad_pq1;
        highp vec2 grad_e1 = -grad_pq1 * clamp1;
        highp float grad_clamp1 = -dot(e1, grad_pq1);

        highp vec2 grad_v2 = grad_pq2;
        highp vec2 grad_e2 = -grad_pq2 * clamp2;
        highp float grad_clamp2 = -dot(e2, grad_pq2);

        highp float grad_div1 = grad_clamp_fun(div1, 0.0, 1.0, grad_clamp1);

        highp float grad_div2 = grad_clamp_fun(div2, 0.0, 1.0, grad_clamp2);

        highp float grad_dot1 = grad_div1 * one_over_dot1;

        highp float grad_dot2 = grad_div2 * one_over_dot2;

        grad_v1 += e1 * grad_dot1;
        grad_e1 += v1 * grad_dot1;

        grad_v2 += e2 * grad_dot2;
        grad_e2 += v2 * grad_dot2;

        grad_uv += grad_v1 + grad_v2;

    } else {
        grad_pq0 += 2.0 * pq0 * incoming_grad;
        distance_sq = dot_pq0_pq0;
    }

    highp vec2 grad_v0 = grad_pq0;
    highp float grad_clamp = -dot(grad_pq0, e0);
    highp float grad_div0 = grad_clamp_fun(div0, 0.0, 1.0, grad_clamp);
    highp float grad_dot0 = grad_div0 * one_over_dot0;
    grad_v0 += e0 * grad_dot0;
    grad_uv += grad_v0;
    highp float sdf_val = sqrt(distance_sq) * poly_sign;

    return vec3(sdf_val, grad_uv / (2.0 * sdf_val));
}
#endif


highp uvec2 to_offset_size(highp uint combined) {
    // note: offset (x coord) is 24 bit -> we have to use highp
    return uvec2(uint(combined >> 8), uint(combined & 255u));
}

highp uvec2 fetch_raw_geometry_data(lowp uint sampler_index, highp uint index, highp uint texture_layer)
{

    // for constants::data_size: 128u, 256u, 512u, 1024u
    // mediump ivec3 dict_px = ivec3(int(index & ((128u<<sampler_index)-1u)), int(index >> (7u+sampler_index)), texture_layer);

     // for constants::data_size: 64u, 128u, 256u, 512u
    mediump ivec3 dict_px = ivec3(int(index & ((64u<<sampler_index)-1u)), int(index >> (6u+sampler_index)), texture_layer);

    switch (sampler_index) {
        case 0u:
            return texelFetch(geometry_buffer_sampler_0, dict_px, 0).rg;
            break;
        case 1u:
            return texelFetch(geometry_buffer_sampler_1, dict_px, 0).rg;
            break;
        case 2u:
            return texelFetch(geometry_buffer_sampler_2, dict_px, 0).rg;
            break;
        default:
            return texelFetch(geometry_buffer_sampler_3, dict_px, 0).rg;
    }

    return uvec2(0u);
}

mediump float float_zoom_interpolation()
{
    // 3d
    highp float dist_camera = length(var_pos_cws.xyz); // TODO divide by cosinus to view direction -> goal: stretched if viewed at steep angle
    // cosinus should not be 0 -> abflachen
    // 2d
    // highp float dist_camera = length(var_pos_cws.xy);

    // TODO move error_threshold_px to camera_config
    // const highp float error_threshold_px = 1.0 / 0.1;
    const highp float error_threshold_px = 1.0 / 0.5;
    // const highp float error_threshold_px = 1.0 / 2.0;

    const highp float sqrt2 = 1.414213562373095;
    const highp float cEarthCircumference = 40075016.685578486;
    const highp float tile_size = 256.0;

    highp float camera_factors = camera.viewport_size.y * 0.5 * camera.distance_scaling_factor;
    highp float static_factors = camera_factors * sqrt2 * cEarthCircumference / tile_size;
    highp float z = log2(static_factors / dist_camera / error_threshold_px)+1.0;

    return clamp(z, 0.0, float(max_zoom));
}

highp uvec3 u32_2_to_u16_u24_u24(highp uvec2 data){
    highp uvec3 res;

    res.x = data.x >> 16;

    res.y = (data.x & ((1u << 16) - 1u)) << 8;
    res.y = res.y | data.y >> 24;

    res.z = data.y & ((1u << 24) - 1u);

    return res;
}

mediump ivec2 to_dict_pixel_64(mediump uint hash) {
    return ivec2(int(hash & 63u), int(hash >> 6u));
}

void parse_style(out LayerStyle style, highp uint style_index, mediump float zoom_offset, mediump float zoom_blend, lowp vec4 ortho_color, bool is_polygon)
{
    // calculate an integer zoom offset for lower and higher style indices and clamp
    lowp int zoom_offset_lower = max(int(floor(zoom_offset-1.0)), -mipmap_levels+1);
    lowp int zoom_offset_higher = max(int(floor(zoom_offset-0.0)), -mipmap_levels+1);

    highp uint style_index_lower = uint(int(style_index) + zoom_offset_lower);
    highp uint style_index_higher = uint(int(style_index) + zoom_offset_higher);

    // get the actual data
    highp uvec4 style_data_lower = texelFetch(styles_sampler, ivec2(to_dict_pixel_64(style_index_lower)), 0);
    highp uvec4 style_data_higher = texelFetch(styles_sampler, ivec2(to_dict_pixel_64(style_index_higher)), 0);

    ///////////////////////////////////////
    // colors
    lowp vec4 color_lower = vec4((style_data_lower.r & 4278190080u) >> 24, (style_data_lower.r & 16711680u) >> 16, (style_data_lower.r & 65280u) >> 8, style_data_lower.r & 255u) / vec4(255.0f);
    lowp vec4 color_higher = vec4((style_data_higher.r & 4278190080u) >> 24, (style_data_higher.r & 16711680u) >> 16, (style_data_higher.r & 65280u) >> 8, style_data_higher.r & 255u) / vec4(255.0f);

    ///////////////////////////////////////
    // line_width
    // saved as tile_extent dependent
    // by dividing by tile_extent we get the width we want to draw
    // by further dividing the tile_extent by 2^zoom_offset, we reduce the tile_extent and increase the line width.
    // we have to increase the zoom_offset by one in order to use the same size as in the preprocessor
    // by using inv_tile_extent and 0.5 as base for pow, we essentially convert a division to a multiplication operation
    mediump float zoomed_tile_extent_lower = inv_tile_extent * pow(0.5, float(zoom_offset_lower+1));
    mediump float zoomed_tile_extent_higher = inv_tile_extent * pow(0.5, float(zoom_offset_higher+1));

    lowp float outline_width_lower = float(style_data_lower.b) * style_precision_mult * zoomed_tile_extent_lower;
    lowp float outline_width_higher = float(style_data_higher.b) * style_precision_mult * zoomed_tile_extent_higher;

    ///////////////////////////////////////
    // actual mix lower and higher style and store info in layerstyle
    style.index = style_index; // setting index to index from geometry -> needed to only parse new styles
    style.color = mix(color_lower, color_higher, zoom_blend);
    if(is_polygon)
        style.color *= ortho_color;
    style.line_width = mix(outline_width_lower, outline_width_higher, zoom_blend);
    // style.dash_info; // TODO
    // style.line_caps; // TODO
}

void alpha_blend(inout lowp vec4 pixel_color, LayerStyle style, lowp float intersection_percentage)
{
    pixel_color = pixel_color + ((1.0-pixel_color.a) * style.color * intersection_percentage);
}


void merge_intersection(inout mediump float accumulated, mediump float current_sample)
{
    // TODO try max(accumulated, current_sample) instead of min(1, sum(accumulated + current_sample))
    accumulated = min(1.0, accumulated + current_sample);
    // accumulated = max(accumulated, current_sample); // -> worsens performance -> not sure if more accurated though
}



highp float calculate_aa_half_radius(highp vec2 uv)
{
  //////////////////////
  // jacobian determinant
  //////////////////////
  highp vec2 uv_x = dFdx(uv);
  highp vec2 uv_y = dFdy(uv);

  highp mat2 jacobian = mat2(uv_x.x,uv_x.y,uv_y.x,uv_y.y);

  // return 0.0;
  return sqrt(abs(determinant(jacobian))) / 2.0;
}

void main() {
#if CURTAIN_DEBUG_MODE == 2
    if (is_curtain == 0.0) {
        discard;
    }
#endif


    // highp float depth = depthWSEncode2n8(length(var_pos_cws));
    highp float depth = length(var_pos_cws);



    lowp vec3 debug_cacscade_layer = vec3(0,0,0);// DEBUG -> which cascade is being used
    lowp vec3 debug_cell_size = vec3(0,0,0);// DEBUG -> how many triangles per cell
    lowp vec3 debug_index_buffer_start = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp vec3 debug_index_buffer_size = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp vec3 debug_texture_layer = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp int debug_draw_calls = 0; // DEBUG



    // get grid acceleration structure data
    highp uvec3 tile_id = var_tile_id;
    highp vec2 uv = var_uv;


    lowp ivec2 dict_px;
    highp uvec2 offset_size = uvec2(0u);

    lowp vec4 pixel_color = vec4(0.0f, 0.0, 0.0f, 0.0f);

    // openmaptile
    lowp vec3 background_color = vec3(242.0f/255.0f, 239.0f/255.0f, 233.0f/255.0f);
    // qwant
    // lowp vec3 background_color = vec3(248.0f/255.0f, 248.0f/255.0f, 248.0f/255.0f);
    // osm-bright
    // lowp vec3 background_color = vec3(248.0f/255.0f, 244.0f/255.0f, 240.0f/255.0f);

    // ORTHO color -> note we wrap tile_id in uvec3 and use ortho_uv, to not change those variables for the vector color
    highp vec2 ortho_uv = uv;
    highp uvec3 temp_tile_id = uvec3(tile_id);
    decrease_zoom_level_until(temp_tile_id, ortho_uv, texelFetch(instanced_texture_zoom_sampler, ivec2(instance_id, 0), 0).x);
    highp float texture_layer_f = float(texelFetch(instanced_texture_array_index_sampler, ivec2(instance_id, 0), 0).x);
    lowp vec4 ortho_color = vec4(texture(texture_sampler, vec3(ortho_uv, texture_layer_f)).rgb, 1.0);
    ortho_color = mix(ortho_color, conf.material_color, conf.material_color.a);

    // VECTOR color
    decrease_zoom_level_until(tile_id, uv, texelFetch(instanced_texture_zoom_sampler_vector, ivec2(instance_id, 0), 0).x);


    mediump float float_zoom = float_zoom_interpolation();
    // float_zoom = tile_id.z;     // DEBUG
    lowp uint mipmap_level = uint(clamp(int(ceil(float(tile_id.z)-float_zoom)), 0,mipmap_levels-1));

    highp uvec2 texture_layer = texelFetch(instanced_texture_array_index_sampler_vector, ivec2(instance_id, mipmap_level), 0).xy;

    decrease_zoom_level_until(tile_id, uv, tile_id.z - mipmap_level);
    mediump float zoom_offset = float_zoom-float(tile_id.z);



    highp float aa_half_radius = calculate_aa_half_radius(uv);

    if(texture_layer.x != bit_mask_ones && texture_layer.y != bit_mask_ones) // check for valid data
    {

        // acceleration_grid_sampler contains the offset and the number of triangles of the current grid cell
        highp ivec2 grid_lookup = ivec2(grid_size*uv);
        offset_size = to_offset_size(texelFetch(acceleration_grid_sampler, ivec3(grid_lookup.x, grid_lookup.y, texture_layer.x & layer_mask),0).r);

#if SDF_MODE == 0
        highp vec2 aa_sample_positions[n_aa_samples]; // uv space
        calculate_sample_positions(aa_sample_positions, uv, grid_lookup);
#endif
#if SDF_MODE == 1
        highp vec2 aa_sample_multipliers[n_aa_samples]; // uv space
        calculate_sample_multipliers(aa_sample_multipliers);
        highp vec2 duvdx = dFdx(uv);
        highp vec2 duvdy = dFdy(uv);
#endif




        // using the grid data we now want to traverse all triangles referenced in grid cell and draw them.
        if(offset_size.y != uint(0)) // only if we have data here
        {
            // lowp vec3 raw_grid = vec3(float(offset_size.y),0,0);// DEBUG
            lowp vec3 raw_grid = vec3(1,0,0);// DEBUG
            lowp ivec2 grid_cell = ivec2(grid_lookup.x, grid_lookup.y);
            lowp vec3 cells = color_from_id_hash(uint(grid_cell.x ^ grid_cell.y)); // DEBUG
            lowp vec2 grid_cell_float = vec2(grid_cell);



            // get the buffer index and extract the correct texture_layer.y
            lowp uint sampler_buffer_index = (texture_layer.y & ((bit_mask_ones << sampler_offset))) >> sampler_offset;
            texture_layer.y = texture_layer.y & layer_mask;

            { // DEBUG
                debug_texture_layer = color_from_id_hash(texture_layer.y);

                {
                    if(sampler_buffer_index == 0u)
                        debug_cacscade_layer = vec3(0,1,0); // green
                    else if(sampler_buffer_index == 1u)
                        debug_cacscade_layer = vec3(1,1,0); // yellow
                    else if(sampler_buffer_index == 2u)
                        debug_cacscade_layer = vec3(1,0.5,0); // orange
                    else if(sampler_buffer_index == 3u)
                        debug_cacscade_layer = vec3(1,0,0); // red
                    else
                        debug_cacscade_layer = vec3(1,0,1); // purple -> should never happen -> unrecognized index
                }

                {
                    if(offset_size.y < 32u)
                        debug_cell_size = vec3(0,1,0); // green
                    else if(offset_size.y < 64u)
                        debug_cell_size = vec3(1,1,0); // yellow
                    else if(offset_size.y  < 128u)
                        debug_cell_size = vec3(1,0.5,0); // orange
                    else if(offset_size.y < 255u)
                        debug_cell_size = vec3(1,0,0); // red
                    else
                        debug_cell_size = vec3(1,0,1); // purple -> should never happen -> unrecognized index
                }

                debug_index_buffer_start = color_from_id_hash(offset_size.x);
                debug_index_buffer_size = color_from_id_hash(offset_size.y);
            } // DEBUG END

            LayerStyle style;
            style.index = -1u;
            style.color = vec4(0.0);
            style.line_width = 0.0;
            style.dash_info = vec2(1.0, 0.0);
            style.line_caps = 0;

            mediump float intersection_percentage = 0.0;
            mediump float zoom_blend = fract(zoom_offset);


            // for(highp uint i = offset_size.x; i < offset_size.x + min(6u,offset_size.y); i++) // show only x layers
            // for(highp uint i = offset_size.x+ offset_size.y; i --> offset_size.x ; ) // reverse traversal
            for(highp uint i = offset_size.x; i < offset_size.x + offset_size.y; i++)
            {
                debug_draw_calls = debug_draw_calls + 1;

                highp uvec2 raw_geom_data = fetch_raw_geometry_data(sampler_buffer_index,i, texture_layer.y);
                highp uint style_index = unpack_style_index(raw_geom_data);

                if (style_index != style.index) {
                    // we changed style -> blend previous style, reset layer infos and parse the new style
                    alpha_blend(pixel_color, style, intersection_percentage);
                    intersection_percentage = 0.0;
                    if (pixel_color.a > full_threshold) {
                        break;
                    }

                    parse_style(style, style_index, zoom_offset, zoom_blend, ortho_color, is_polygon(raw_geom_data)); // mix floating zoom levels; for polygons, mul poly color with surface shading texture color
                }


#if SDF_MODE == 1
                { // derivative aa
                    VectorLayerData geom_data = unpack_data(raw_geom_data, grid_cell_float); // TODO move below code inside unpack_data


                    highp vec3 dist_and_grad = sdf_with_grad(geom_data, uv, 1.0);


                    highp float accumulated = 0.0;


                    for (lowp int j = 0; j < n_aa_samples; ++j)
                    {
                        highp float dDist_dx = dot(dist_and_grad.yz, duvdx);
                        highp float dDist_dy = dot(dist_and_grad.yz, duvdy);

                        highp float d = dist_and_grad.x + aa_sample_multipliers[j].x * dDist_dx + aa_sample_multipliers[j].y * dDist_dy;
                        if(!geom_data.is_polygon)
                            d = abs(d);
                        d = d - style.line_width;




                        // mediump float percentage = 1.0 - smoothstep(-aa_half_radius, aa_half_radius, d);
                        mediump float percentage = 1.0 - step(0.0,d);
                        accumulated += percentage;
                    }
                    merge_intersection(intersection_percentage, accumulated * division_by_n_samples);

                }
#endif
#if SDF_MODE == 0
                { // naive aa
                    VectorLayerData geom_data = unpack_data(raw_geom_data, grid_cell_float);
                    SDFData prepared_sdf_data = prepare_sd_Line_Triangle(geom_data);

                    highp float accumulated = 0.0;
                    for (lowp int j = 0; j < n_aa_samples; ++j)
                    {
                        highp float d = sd_Line_Triangle(aa_sample_positions[j], prepared_sdf_data, geom_data.is_polygon) - style.line_width;
                        mediump float percentage = 1.0 - smoothstep(-aa_half_radius, aa_half_radius, d);
                        accumulated += percentage;
                    }
                    merge_intersection(intersection_percentage, accumulated * division_by_n_samples); // try current min(1, sum(geometries)) and alternative max(geometries)
                }
#endif




            }

            // blend the last layer (this only does something if we did not break out of the loop early)
            alpha_blend(pixel_color, style, intersection_percentage);
        }
    }


    // blend with background color and write pixel color to output
    texout_albedo = vec3(pixel_color.rgb + ((1.0-pixel_color.a)*background_color * ortho_color.rgb));




    if (conf.overlay_mode > 199u && conf.overlay_mode < 300u) {
        lowp vec3 zoom_debug_color =  color_from_id_hash(uint(float_zoom));

        lowp vec3 overlay_color = vec3(0.0);
        switch(conf.overlay_mode) {
            case 200u: overlay_color = vec3(uv, 0.0);break;
            case 201u: overlay_color = debug_cacscade_layer;break;
            case 202u: overlay_color = debug_cell_size;break;
            case 204u: overlay_color = debug_index_buffer_start;break;
            case 205u: overlay_color = debug_index_buffer_size;break;
            case 206u: overlay_color = debug_texture_layer;break;
            case 207u: overlay_color = vec3(pixel_color.a, 0.0, 0.0);break;

            // case 208u: handled below;
            // case 209u: overlay_color = vec3(float_zoom, 0.0, 0.0);break;
            case 209u: overlay_color = mix(vec3(1,1,1), zoom_debug_color, 1.0-fract(float_zoom));break;
            default: overlay_color = vertex_color;
        }
        texout_albedo.rgb = mix(texout_albedo.rgb, overlay_color, conf.overlay_strength);


        if(conf.overlay_mode == 208u)
        {
            lowp float upper_limit = conf.overlay_strength * 255.0;
            if(debug_draw_calls > int(upper_limit))
                texout_albedo.rgb = vec3(1.0, 0.0, 0.0);
            else if (debug_draw_calls == 0)
                texout_albedo.rgb = vec3(0.0, 0.0, 0.0);
            else
                texout_albedo.rgb = vec3(0.0, float(debug_draw_calls) / upper_limit, 0.0);
        }
        else if(conf.overlay_mode == 210u)
        {
            // comparison between float zoom and tile zoom
            texout_albedo.rgb = mix(zoom_debug_color, color_from_id_hash(uint(var_tile_id.z)), conf.overlay_strength);
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
        texout_albedo.rgb = mix(texout_albedo.rgb, overlay_color, conf.overlay_strength);
    }

#if CURTAIN_DEBUG_MODE == 1
    if (is_curtain > 0.0) {
        texout_albedo = vec3(1.0, 0.0, 0.0, 1.0);
        return;
    }
#endif

}
