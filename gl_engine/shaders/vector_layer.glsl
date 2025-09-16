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

#line 10020

#define PI 3.1415926535

// TODO only here for testing purposes -> remove the define again
#define SDF_MODE 0



// SDF_MODE 0: naive multisample antialiasing (golden -> but worse performnace)
// SDF_MODE 1: derivative multisample antialiasing
#ifndef SDF_MODE
#define SDF_MODE 1
#endif

// 0 ortho only
// 1 mixed
// 2 vector only
#ifndef VIEW_MODE
#define VIEW_MODE 1
#endif


uniform highp usampler2D styles_sampler;

uniform highp usampler2DArray geometry_buffer_sampler_0;
uniform highp usampler2DArray geometry_buffer_sampler_1;
uniform highp usampler2DArray geometry_buffer_sampler_2;


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

#ifndef n_multisamples
const lowp int n_aa_samples_row_cols = 4;
#else
const lowp int n_aa_samples_row_cols = n_multisamples;
#endif
const lowp int n_aa_samples = n_aa_samples_row_cols*n_aa_samples_row_cols;
const lowp float aa_sample_dist = 1.0; // 0.5 goes from -0.25 to +0.25 of the current uv coordinate

const mediump float division_by_n_samples = 1.0 / float(n_aa_samples);
const lowp float style_precision_mult = 1.0 / float(style_precision);
const mediump float inv_tile_extent = 1.0 / tile_extent;


const highp vec2 cell_size_uv = vec2(1.0) / vec2(grid_size);
const highp vec2 aa_cell_overlap_uv = vec2(aa_border) * vec2(cell_size_uv);

const highp uint layer_mask = ((1u << sampler_offset) - 1u);
const highp uint bit_mask_ones = -1u;

const highp float full_threshold = 0.99;

// openmaptile
const lowp vec3 background_color = vec3(242.0f/255.0f, 239.0f/255.0f, 233.0f/255.0f);
// qwant
// const lowp vec3 background_color = vec3(248.0f/255.0f, 248.0f/255.0f, 248.0f/255.0f);
// osm-bright
// const lowp vec3 background_color = vec3(248.0f/255.0f, 244.0f/255.0f, 240.0f/255.0f);


struct VectorLayerData{
    highp vec2 a;
    highp vec2 b;
    highp vec2 c;

    highp uint style_index;
    bool is_polygon;

    bool line_cap0;
    bool line_cap1;
};

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

    bool line_cap0;
    bool line_cap1;
};
#endif

struct LayerStyle {
    highp uint index;
    lowp vec4 color;
    lowp float line_width;
    lowp vec2 dash_info;
    lowp bool round_line_caps;
};



/////////////////////////////////////////////
// constants for data packing/unpacking

// defined in c++ code
// const lowp int all_bits = 32; // per output channel
// const lowp int coordinate_bits = 8;

const highp vec2 cell_size = vec2(tile_extent) / vec2(grid_size);


const highp uint coordinate_bitmask = (1u << coordinate_bits_polygons) - 1u;
const highp uint coordinate_bitmask_lines = (1u << coordinate_bits_lines) - 1u;
const lowp int remaining_coordinate_bits_lines = coordinate_bits_lines - coordinate_bits_polygons;
const highp uint remaining_coordinate_bitmask_lines = (1u << remaining_coordinate_bits_lines) - 1u;
const highp uint is_polygon_bitmask = (1u << style_bits);

const lowp int coordinate_shift1 = all_bits - coordinate_bits_polygons;
const lowp int coordinate_shift2 = all_bits - (2 * coordinate_bits_polygons);
const lowp int coordinate_shift3 = all_bits - (3 * coordinate_bits_polygons);
const lowp int coordinate_shift4 = all_bits - (4 * coordinate_bits_polygons);

const highp uint coordinate_bitmask_shift1 = coordinate_bitmask << coordinate_shift1;
const highp uint coordinate_bitmask_shift2 = coordinate_bitmask << coordinate_shift2;
const highp uint coordinate_bitmask_shift3 = coordinate_bitmask << coordinate_shift3;
const highp uint coordinate_bitmask_shift4 = coordinate_bitmask << coordinate_shift4;
const highp uint remaining_coordinates_bitmask_shift = remaining_coordinate_bitmask_lines << remaining_coordinate_bits_lines;

// for packed.y -> we need to divide coordinate_bits_polygons by 2
const lowp int coordinate_shift1_lines = all_bits - int(0.5 * float(coordinate_bits_polygons));
const lowp int coordinate_shift2_lines = all_bits - int(1.0 * float(coordinate_bits_polygons));
const lowp int coordinate_shift3_lines = all_bits - int(1.5 * float(coordinate_bits_polygons));
const lowp int coordinate_shift4_lines = all_bits - int(2.0 * float(coordinate_bits_polygons));

const highp int cell_width_polygons = int((tile_extent * scale_polygons) / grid_size.x);
const highp int cell_width_lines = int((tile_extent * scale_lines) / grid_size.x);

const highp int max_cell_width_polygons = (1 << (coordinate_bits_polygons));
const highp int geometry_offset_polygons = (max_cell_width_polygons - cell_width_polygons) / 2;
const highp int max_cell_width_line = (1 << (coordinate_bits_lines));
const highp int geometry_offset_line = (max_cell_width_line - cell_width_lines) / 2;

const highp uint style_bitmask = ((1u << style_bits) - 1u);

const highp uint line_cap0_mask = 1u << (style_bits + 1);
const highp uint line_cap1_mask = 1u << (style_bits + 2);


// Style unpacking
const lowp uint style_width_offset = 17u;

const lowp uint style_dash_ratio_offset = 9u;
const lowp uint style_dash_ratio_mask = (1u << 17u) - 1u;

const lowp uint style_dash_sum_offset = 1u;
const lowp uint style_dash_sum_mask = (1u << 9u) - 1u;

const lowp uint style_cap_mask = 1u;


// end constants for data packing/unpacking
/////////////////////////////////////////////

highp uvec2 pack_vectorlayer_data(VectorLayerData data)
{
    highp uvec2 packed_data;

    ivec2 a = ivec2(data.a);
    ivec2 b = ivec2(data.b);
    ivec2 c = ivec2(data.c);

    if (data.is_polygon) {
        a += geometry_offset_polygons;
        b += geometry_offset_polygons;
        c += geometry_offset_polygons;
    } else {
        a += geometry_offset_line;
        b += geometry_offset_line;
    }

    packed_data.x = uint(a.x) << coordinate_shift1;
    packed_data.x = packed_data.x | ((uint(a.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.x = packed_data.x | ((uint(b.x) & coordinate_bitmask) << coordinate_shift3);
    packed_data.x = packed_data.x | ((uint(b.y) & coordinate_bitmask) << coordinate_shift4);

    if (data.is_polygon) {
        packed_data.y = uint(c.x) << coordinate_shift1;
        packed_data.y = packed_data.y | ((uint(c.y) & coordinate_bitmask) << coordinate_shift2);
    } else {
        packed_data.y = ((uint(a.x) >> coordinate_bits_polygons) << coordinate_shift1_lines);
        packed_data.y = packed_data.y | ((uint(a.y) >> coordinate_bits_polygons) << coordinate_shift2_lines);
        packed_data.y = packed_data.y | ((uint(b.x) >> coordinate_bits_polygons) << coordinate_shift3_lines);
        packed_data.y = packed_data.y | ((uint(b.y) >> coordinate_bits_polygons) << coordinate_shift4_lines);
    }

    highp uint is_polygon = (data.is_polygon ? 1u : 0u) << style_bits;

    packed_data.y = packed_data.y | is_polygon | data.style_index;

    return packed_data;
}

highp uint unpack_style_index(highp uvec2 packed_data)
{
    return (packed_data.y & style_bitmask) * uint(max_zoom+1);
}

bool is_polygon(highp uvec2 packed_data)
{
    return (packed_data.y & is_polygon_bitmask) != 0u;
}

VectorLayerData unpack_data(highp uvec2 packed_data, lowp vec2 grid_cell)
{
    VectorLayerData unpacked_data;

    ivec2 a;
    ivec2 b;
    ivec2 c = ivec2(0,0);

    a.x = int((packed_data.x & (coordinate_bitmask_shift1)) >> coordinate_shift1);
    a.y = int((packed_data.x & (coordinate_bitmask_shift2)) >> coordinate_shift2);
    b.x = int((packed_data.x & (coordinate_bitmask_shift3)) >> coordinate_shift3);
    b.y = int((packed_data.x & (coordinate_bitmask_shift4)) >> coordinate_shift4);

    highp uvec2 c_u;
    c_u.x = (packed_data.y & (coordinate_bitmask_shift1)) >> coordinate_shift1;
    c_u.y = (packed_data.y & (coordinate_bitmask_shift2)) >> coordinate_shift2;

    unpacked_data.style_index = unpack_style_index(packed_data);

    unpacked_data.is_polygon = is_polygon(packed_data);

    if (unpacked_data.is_polygon) {
        a -= geometry_offset_polygons;
        b -= geometry_offset_polygons;
        c = ivec2(c_u) - geometry_offset_polygons;
    } else {
        // unpack most significant coordinates of the line and add them to the unpacked lines
        a.x = a.x | int(((c_u.x & (remaining_coordinates_bitmask_shift)) << remaining_coordinate_bits_lines));
        a.y = a.y | int(((c_u.x & remaining_coordinate_bitmask_lines) << coordinate_bits_polygons));
        b.x = b.x | int(((c_u.y & (remaining_coordinates_bitmask_shift)) << remaining_coordinate_bits_lines));
        b.y = b.y | int(((c_u.y & remaining_coordinate_bitmask_lines) << coordinate_bits_polygons));

        a -= geometry_offset_line;
        b -= geometry_offset_line;

        unpacked_data.line_cap0 = (packed_data.y & line_cap0_mask) != 0u;
        unpacked_data.line_cap1 = (packed_data.y & line_cap1_mask) != 0u;
    }

    highp float tile_scale = float(scale_lines) * (1.0-float(unpacked_data.is_polygon)) + float(scale_polygons) * float(unpacked_data.is_polygon);

    highp vec2 cell_offset = grid_cell * cell_size * tile_scale;
    highp float division_extent = 1.0 / (tile_extent * tile_scale);

    unpacked_data.a = (vec2(a) + cell_offset) * division_extent;
    unpacked_data.b = (vec2(b) + cell_offset) * division_extent;
    unpacked_data.c = (vec2(c) + cell_offset) * division_extent;

    return unpacked_data;
}


// for the base unittest we are not interested in testing if the polygons scaling and uv scaling works
// we only are interested if the data we send to the gpu is correctly encoded and decoded without anly loss of precision
VectorLayerData normalize_unpack_for_unittest(VectorLayerData unpacked_data, lowp vec2 grid_cell)
{
    highp float tile_scale = float(scale_lines) * (1.0-float(unpacked_data.is_polygon)) + float(scale_polygons) * float(unpacked_data.is_polygon);

    highp vec2 cell_offset = grid_cell * cell_size * tile_scale;
    highp float division_extent = 1.0 / (tile_extent * tile_scale);

    unpacked_data.a = (vec2(unpacked_data.a) / division_extent) - cell_offset;
    unpacked_data.b = (vec2(unpacked_data.b) / division_extent) - cell_offset;
    unpacked_data.c = (vec2(unpacked_data.c) / division_extent) - cell_offset;

    unpacked_data.style_index /= uint(max_zoom+1);

    return unpacked_data;
}





#if SDF_MODE == 0
void calculate_sample_positions(out highp vec2 aa_sample_positions[n_aa_samples], highp vec2 uv, highp ivec2 grid_cell)
{
    if(n_aa_samples <= 1)
    {
        aa_sample_positions[0] = uv;
        return;
    }

    highp vec2 min_cell = vec2(grid_cell) * cell_size_uv;
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

SDFData prepare_sd_Line_Triangle(VectorLayerData geom_data)
{
    SDFData data;

    data.vertices = geom_data;

    data.e0 = data.vertices.b-data.vertices.a;
    data.dot_e0 = dot(data.e0,data.e0);

    data.line_cap0 = geom_data.line_cap0;
    data.line_cap1 = geom_data.line_cap1;

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
highp float sd_Line_Triangle( in highp vec2 uv, SDFData data, bool triangle, highp float line_width, lowp vec2 dash_info, bool round_line_caps)
{
    highp vec2 v0 = uv - data.vertices.a;
    highp vec2 v1 = uv - data.vertices.b;
    highp float h = clamp( dot(v0,data.e0)/data.dot_e0, 0.0, 1.0 );
    highp vec2 pq0 = v0 - data.e0*h;

    highp float poly_sign = 1.0;
    highp float mask = 1.0;
    highp float result = 1.0;

    if(triangle)
    {
        highp vec2 v2 = uv - data.vertices.c;
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

    }
    else{

        highp float line_length = length(data.e0);
        highp float round_factor = 1.0 + 450.0*float(!round_line_caps);
        const highp float round_factor_dashes = 5000.0;

        highp float amount_dash_gap_pairs = ceil(line_length/dash_info.y);
        // + 0.01 -> small delta to remove artifacts if there shouldn't be any dashes
        highp float dash_period = cos(PI*h*amount_dash_gap_pairs*2.0)+cos((1.0-dash_info.x)*PI)+0.01;
        highp float dashes = tanh((dash_period)*round_factor_dashes);


        float line_endings = 1.0;
        if(!round_line_caps)
        {
            if(data.line_cap0)
                line_endings *= dot(normalize(data.e0), v0);
            if(data.line_cap1)
                line_endings *= dot(normalize(-data.e0), v1);
        }
        line_endings = (sign(line_endings)+1.0) / 2.0;
        // line_endings = 1.0;


        mask = line_endings*dashes;
        poly_sign = 1.0;
        result = dot(pq0,pq0);
    }

    return sqrt(result)*poly_sign - (line_width * mask);

}
#endif
#if SDF_MODE == 1
void calculate_sample_multipliers(out highp vec2 aa_sample_multipliers[n_aa_samples])
{
    if(n_aa_samples <= 1)
    {
        aa_sample_multipliers[0] = vec2(0.0);
        return;
    }


    highp float aa_sample_dist_increments = aa_sample_dist / float(n_aa_samples_row_cols);

    highp vec2 start = vec2(-aa_sample_dist / 2.0 + aa_sample_dist_increments / 2.0);
    for (int x = 0; x < n_aa_samples_row_cols; ++x) {
        for (int y = 0; y < n_aa_samples_row_cols; ++y) {
            lowp int index = x*n_aa_samples_row_cols + y;
            aa_sample_multipliers[index] = start + vec2(x,y) * vec2(aa_sample_dist_increments);
        }
    }

}

highp float grad_clamp_fun(highp float v, highp float min, highp float max, highp float incoming_grad)
{
    // return incoming_grad * // TODO convert to * instead of if

    // return float(v > min && v < max) * incoming_grad;
    // return incoming_grad * step(v,min) * step(max,v);

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

    // for constants::data_size: 128u, 256u, 512u
    mediump ivec3 dict_px = ivec3(int(index & ((128u<<sampler_index)-1u)), int(index >> (7u+sampler_index)), texture_layer);

     // for constants::data_size: 64u, 128u, 256u
    // highp ivec3 dict_px = ivec3(int(index & ((64u<<sampler_index)-1u)), int(index >> (6u+sampler_index)), texture_layer);

    switch (sampler_index) {
        case 0u:
            return texelFetch(geometry_buffer_sampler_0, dict_px, 0).rg;
            break;
        case 1u:
            return texelFetch(geometry_buffer_sampler_1, dict_px, 0).rg;
            break;
        default:
            return texelFetch(geometry_buffer_sampler_2, dict_px, 0).rg;
    }

    return uvec2(0u);
}


mediump ivec2 to_dict_pixel_128(mediump uint hash) {
    return ivec2(int(hash & 127u), int(hash >> 7u));
}

void parse_style(out LayerStyle style, highp uint style_index, mediump float zoom_offset, mediump float zoom_blend, lowp vec4 ortho_color, mediump float cos_smoothing_factor, bool is_polygon)
{
    // calculate an integer zoom offset for lower and higher style indices and clamp
    lowp int zoom_offset_lower = max(int(floor(zoom_offset-1.0)), -mipmap_levels+1);
    lowp int zoom_offset_higher = max(int(floor(zoom_offset-0.0)), -mipmap_levels+1);

    highp uint style_index_lower = uint(int(style_index) + zoom_offset_lower);
    highp uint style_index_higher = uint(int(style_index) + zoom_offset_higher);

    // get the actual data
    highp uvec4 style_data_lower = texelFetch(styles_sampler, ivec2(to_dict_pixel_128(style_index_lower)), 0);
    highp uvec4 style_data_higher = texelFetch(styles_sampler, ivec2(to_dict_pixel_128(style_index_higher)), 0);

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

    lowp float outline_width_lower = float(style_data_lower.g >> style_width_offset) * style_precision_mult * zoomed_tile_extent_lower;
    lowp float outline_width_higher = float(style_data_higher.g >> style_width_offset) * style_precision_mult * zoomed_tile_extent_higher;

    ///////////////////////////////////////
    // actual mix lower and higher style and store info in layerstyle
    style.index = style_index; // setting index to index from geometry -> needed to only parse new styles
    // calculate color by blending lower/higher
    // and multiply orthocolor or a line smoothing factor from viewing angle depending if line or polygon
    style.color = mix(color_lower, color_higher, zoom_blend) * mix(vec4(cos_smoothing_factor), ortho_color, float(is_polygon));

    style.line_width = mix(outline_width_lower, outline_width_higher, zoom_blend);

    lowp float dash_ratio_lower = float((style_data_lower.g & style_dash_ratio_mask) >> style_dash_ratio_offset) * style_precision_mult;
    lowp float dash_ratio_higher = float((style_data_higher.g & style_dash_ratio_mask) >> style_dash_ratio_offset) * style_precision_mult;
    lowp float dash_sum_lower = float((style_data_lower.g & style_dash_sum_mask) >> style_dash_sum_offset) * style_precision_mult;
    lowp float dash_sum_higher = float((style_data_higher.g & style_dash_sum_mask) >> style_dash_sum_offset) * style_precision_mult;

    // TODO maybe one mix is faster?
    style.dash_info = vec2(mix(dash_ratio_lower, dash_ratio_higher, zoom_blend), mix(dash_sum_lower, dash_sum_higher, zoom_blend));

    // for line caps we only really need one style since we assume that they do not change between zoom levels
    style.round_line_caps = (style_data_higher.g & style_cap_mask) == 1u;
}

lowp float hit_percentage(highp uint intersections)
{
    if(intersections == 0u)
        return 0.0; // TODO check if early exit here is faster

    lowp uint bits_hit = 0u;

    // NOTE: this is essentially what bitCount method does, but bitCount does not exist for webassembly builds.
    for (lowp int i = 0; i < n_aa_samples; ++i)
    {
        bits_hit += (intersections >> i) & 1u;
    }

    return float(bits_hit) / float(n_aa_samples);
}


void alpha_blend(inout lowp vec4 pixel_color, LayerStyle style, highp uint intersections)
{
    // we store which sample has hit the geometry -> if two geometries hit the same sample we only store one hit
    // this should make multisample anti-aliasing a bit better
    lowp float intersection_percentage = hit_percentage(intersections);

    pixel_color = pixel_color + ((1.0-pixel_color.a) * style.color * intersection_percentage);
}


struct DrawMeta
{
    lowp uint sampler_buffer_index;
    highp uint texture_layer;
    highp uint tile_zoom;
    lowp vec4 ortho_color;
    lowp vec2 grid_cell_float;
    highp vec2 aa_sample_multipliers[n_aa_samples];
    highp vec2 aa_sample_positions[n_aa_samples];
    mediump float cos_smoothing_factor;
    highp vec2 duvdx;
    highp vec2 duvdy;
    mediump float zoom_offset;
    mediump float zoom_blend;
};


bool draw_layer(inout lowp vec4 pixel_color, inout highp uint intersections, inout LayerStyle style, highp vec2 uv, highp uint i, DrawMeta meta)
{
    highp uvec2 raw_geom_data = fetch_raw_geometry_data(meta.sampler_buffer_index, i, meta.texture_layer);
    highp uint style_index = unpack_style_index(raw_geom_data) + meta.tile_zoom;

    if (style_index != style.index) {
        // we changed style -> blend previous style, reset layer infos and parse the new style

        alpha_blend(pixel_color, style, intersections);
        intersections = 0u;
        if (pixel_color.a > full_threshold) {
            return true;
        }

        parse_style(style, style_index, meta.zoom_offset, meta.zoom_blend, meta.ortho_color, meta.cos_smoothing_factor, is_polygon(raw_geom_data)); // mix floating zoom levels; for polygons, mul poly color with surface shading texture color
    }


#if SDF_MODE == 1
    { // derivative aa
        VectorLayerData geom_data = unpack_data(raw_geom_data, meta.grid_cell_float);

        highp vec3 dist_and_grad = sdf_with_grad(geom_data, uv, 1.0);

        highp float dDist_dx = dot(dist_and_grad.yz, meta.duvdx);
        highp float dDist_dy = dot(dist_and_grad.yz, meta.duvdy);

        for (lowp int j = 0; j < n_aa_samples; ++j)
        {
            highp float d = dist_and_grad.x + meta.aa_sample_multipliers[j].x * dDist_dx + meta.aa_sample_multipliers[j].y * dDist_dy;

            // for lines use the abs(d) -> lines should not be able to be negative -> but it is possible with derivative -> we have to correct this
            d = mix(abs(d), d, float(geom_data.is_polygon)) - style.line_width;

            // highp uint geometry_hit = uint(1.0 - step(0.0,d));
            highp uint geometry_hit = uint(d<0.0);// TODO check if faster
            intersections |= geometry_hit << j;
        }
    }
#endif
#if SDF_MODE == 0
    { // naive aa
        VectorLayerData geom_data = unpack_data(raw_geom_data, meta.grid_cell_float);
        SDFData prepared_sdf_data = prepare_sd_Line_Triangle(geom_data);

        for (lowp int j = 0; j < n_aa_samples; ++j)
        {
            highp float d = sd_Line_Triangle(meta.aa_sample_positions[j], prepared_sdf_data, geom_data.is_polygon, style.line_width, style.dash_info, style.round_line_caps);

            highp uint geometry_hit = uint(1.0 - step(0.0,d));
            intersections |= geometry_hit << j;
        }
    }
#endif

     return false;

}


