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


// SDF_MODE 0: naive multisample antialiasing (golden -> but worse performnace)
// SDF_MODE 1: screenspace transformation -> smoothstep
#ifndef SDF_MODE
#define SDF_MODE 1
#endif

// 0 ortho only
// 1 mixed
// 2 vector only
#ifndef VIEW_MODE
#define VIEW_MODE 1
#endif

// SAMPLE_DISTRIBUTION 0: UNIFORM SAMPLES
// SAMPLE_DISTRIBUTION 1: RANDOM SAMPLES
#ifndef SAMPLE_DISTRIBUTION
#define SAMPLE_DISTRIBUTION 1
#endif

// COVERAGE_SIMPLE
// 0 -> we are using subtraction and multiplication of other halfspaces (correct but costly)
// 1 -> we are only computing the coverage of the biggest halfspace distance (expect aliasing at thin lines) -> less noticeable if we reduce opacity of far away geometry at shallow angles (cos_smoothing_factor)
#ifndef COVERAGE_SIMPLE
#define COVERAGE_SIMPLE 1
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

// NOTE: n_aa_samples_row_cols defines how many samples per uniform row -> they will be squared
// This is also currently married to the random samples -> e.g. setting it to 4 will result in 16 random samples
// the random samples are also currently limited to 16 -> if we want to increase the samples we would also need to increase the random_samples variable
#ifndef n_multisamples
const lowp int n_aa_samples_row_cols = 4;
#else
const lowp int n_aa_samples_row_cols = n_multisamples;
#endif
const lowp int n_aa_samples = n_aa_samples_row_cols*n_aa_samples_row_cols;
// determined at compilation of the shader from cpp side
// const lowp float aa_sample_dist = 1.0; // 0.5 goes from -0.25 to +0.25 of the current uv coordinate

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

const highp vec2[16] random_samples = vec2[16](vec2(-0.12781827,-0.07746551),vec2(-0.56825564,-0.50103748),vec2(0.21763546,0.22721745),vec2(0.06014171,0.00090832),vec2(0.38402289,0.78174722),vec2(0.40367972,-0.20728750),vec2(0.72827647,-0.75817493),vec2(-0.55837659,-0.06544222),vec2(-0.51244102,0.36216960),vec2(-0.19719190,0.25865584),vec2(-0.06226554,-0.12613159),vec2(-0.10024053,0.68217920),vec2(0.39563396,0.31162355),vec2(0.40185253,-0.29009449),vec2(0.32805709,-0.38292970),vec2(-0.92934095,0.41766199));


struct VectorLayerData{
    highp vec2 a;
    highp vec2 b;
    highp vec2 c;

    bvec3 additional_info;

    highp uint style_index;
    bool is_polygon;
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
};
#endif

struct LayerStyle {
    highp uint index;
    lowp vec4 color;
    lowp float line_width;
    lowp vec2 dash_info;
    bool round_line_caps;
};


struct DrawMeta
{
    lowp uint sampler_buffer_index;
    highp uint texture_layer;
    highp uint tile_zoom;
    lowp vec4 ortho_color;
    lowp vec2 grid_cell_float;
#if SDF_MODE == 0
    highp vec2 aa_sample_multipliers[n_aa_samples];
    highp vec2 aa_sample_positions[n_aa_samples];
    highp vec2 duvdx;
    highp vec2 duvdy;
#else
    highp mat3x3 uv2fragspace_normal_matrix;
#endif
    mediump float cos_smoothing_factor;
    mediump float zoom_offset;
    mediump float zoom_blend;
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

const highp uint additonal_info0_mask = 1u << (style_bits + 3);
const highp uint additonal_info1_mask = 1u << (style_bits + 2);
const highp uint additonal_info2_mask = 1u << (style_bits + 1);


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

// TODO we are converting every geometry to uv space and also all the line widths are converted to uv space
// is it possible to stay in tile space and only convert the uv coordinate to tile space? -> this can be done before the loop
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

    unpacked_data.additional_info.x = (packed_data.y & additonal_info0_mask) != 0u;
    unpacked_data.additional_info.y = (packed_data.y & additonal_info1_mask) != 0u;
    unpacked_data.additional_info.z = (packed_data.y & additonal_info2_mask) != 0u;

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

// https://thebookofshaders.com/10/
highp float random (highp vec2 st) {
    return fract(sin(dot(st.xy,vec2(12.9898,78.233))) * 43758.5453123);
}

highp vec2 hash(highp vec2 p) {
    uvec2 res = uvec2(p*4294967296.0);
    res.x ^= res.x >> 16;
    res.y ^= res.y >> 16;
    res.x *= 0x7feb352dU;
    res.y *= 0x7feb352dU;
    res.x ^= res.x >> 15;
    res.y ^= res.y >> 15;
    res.x *= 0x846ca68bU;
    res.y *= 0x846ca68bU;
    res.x ^= res.x >> 16;
    res.y ^= res.y >> 16;

    return vec2(res) / 4294967296.0;
}

highp vec2 random_gaussian_point(ivec2 offset, highp vec2 uv)
{
    ////////////////////////////////
    // Random disk with polar coordinates
    ////////////////////////////////
    // performance worsens by ~1-2 ms (in vienna) but prevents Moiré patterns
    // {
    // float r = random(uv + vec2(0.432823*offset.x, 0.282*offset.y))*aa_sample_dist;
    // float angle = random(uv + vec2(0.5484523*offset.x, 0.81054*offset.y)) * PI*2.0;

    // return vec2(r*cos(angle), r*sin(angle));
    // }

    ////////////////////////////////
    // Gaussian distribution with random start_index+mirrors
    ////////////////////////////////
    // performance worsens by ~1-2 ms (in vienna) but prevents Moiré patterns
    // {
    //     mediump vec2 rand = hash(uv);
    //     // how many different samples can we create from the random samples provided by cpp
    //     const lowp float sample_splits = float(num_random_samples / n_aa_samples);
    //     // get a random start_index
    //     lowp int start_index = int(n_aa_samples)*int(sample_splits*rand.x);

    //     // get the index using start_index and the current x,y position
    //     lowp int index = (start_index + int(offset.x * n_aa_samples_row_cols + offset.y));

    //     // randomly mirror in both x and y direction
    //     return vec2((step(0.5,rand.y)-0.5)*2.0)*random_samples[index];
    // }

    ////////////////////////////////
    // Gaussian distribution that stays the same over all pixels
    ////////////////////////////////
    // -> similar performance to uniform grids, but advantages of prefering the pixel center
    lowp int index = int(offset.x * n_aa_samples_row_cols + offset.y);
    return random_samples[index];

}

#if SDF_MODE == 0
void calculate_samples(inout DrawMeta meta, highp vec2 uv)
{
    if(n_aa_samples <= 1)
    {
        meta.aa_sample_positions[0] = uv;

        return;
    }
    highp vec2 min_cell = meta.grid_cell_float * cell_size_uv;
    highp vec2 max_cell = min_cell + cell_size_uv;
    min_cell -= vec2(aa_cell_overlap_uv);
    max_cell += vec2(aa_cell_overlap_uv);

    highp vec2 grad_u = vec2(meta.duvdx.x, meta.duvdy.x);
    highp vec2 grad_v = vec2(meta.duvdx.y, meta.duvdy.y);

#if SAMPLE_DISTRIBUTION == 0
    highp float aa_sample_dist_increments = aa_sample_dist / float(n_aa_samples_row_cols);
    highp vec2 start = vec2(-aa_sample_dist / 2.0 + aa_sample_dist_increments / 2.0);
#endif

    for (int x = 0; x < n_aa_samples_row_cols; ++x) {
        for (int y = 0; y < n_aa_samples_row_cols; ++y) {
            lowp int index = x * n_aa_samples_row_cols + y;

#if SAMPLE_DISTRIBUTION == 0
            meta.aa_sample_multipliers[index] = start + vec2(x,y) * vec2(aa_sample_dist_increments);
#else
            meta.aa_sample_multipliers[index] = random_gaussian_point(ivec2(x,y), uv);
#endif

            meta.aa_sample_positions[index].x = uv.x + dot(grad_u, meta.aa_sample_multipliers[index]);
            meta.aa_sample_positions[index].y = uv.y + dot(grad_v, meta.aa_sample_multipliers[index]);

            // clip to current cell
            meta.aa_sample_positions[index] = min(max_cell, max(min_cell, meta.aa_sample_positions[index]));
        }
    }
}

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

        highp float amount_dash_gap_pairs = ceil(line_length/dash_info.y);
        // + 0.01 -> small delta to remove artifacts if there shouldn't be any dashes
        highp float dash_period = cos(PI*h*amount_dash_gap_pairs*2.0)+cos((1.0-dash_info.x)*PI)+0.01;
        // tanh is used as a differentiable step function -> all values above 0 are mapped to +1, all below to -1
        // multiplication by big value ensures a quick transition at 0 +/- small delta
        highp float dashes = tanh((dash_period)*500000.0);

        highp float line_endings = 1.0;
        if(!round_line_caps)
        {
            if(data.vertices.additional_info.y)
                line_endings *= dot(normalize(data.e0), v0);
            if(data.vertices.additional_info.z)
                line_endings *= dot(normalize(-data.e0), v1);
        }
        line_endings = (tanh(line_endings*500000.0)+1.0) / 2.0;

        mask = line_endings*dashes;
        result = dot(pq0,pq0);
    }

    return sqrt(result)*poly_sign - (line_width * mask);

}

#else
// screenspace sdf



highp vec3 create_halfspace_with_normal(highp vec2 a, highp vec2 n)
{
    highp float distance = dot(a, n);

    return vec3(n, -distance);
}

highp vec3 create_halfspace(highp vec2 a, highp vec2 b)
{
    highp vec2 e = b - a;
    highp vec2 normal = normalize(vec2(-e.y, e.x));

    return create_halfspace_with_normal(a, normal);
}

highp vec3 create_round_cap(highp vec2 uv, highp vec2 point, highp float line_width)
{
    highp vec2 n = normalize(uv - point);
    highp float dist = dot(n, point + n*line_width);
    return vec3(n, -dist);
}

highp vec3 create_line_segment_end_halfspace(highp vec2 uv, VectorLayerData geom_data, highp vec2 n_line, highp float line_width, bool round_line_caps)
{
    // calculate the sign between current uv and the two normals for line endings
    highp vec2 n_line_end = vec2(n_line.y, -n_line.x);

    highp float dist_a = dot(uv-geom_data.a, -n_line_end);
    highp float dist_b = dot(uv-geom_data.b, n_line_end);

    // step gives us 0 or 1 -> depending which value is higher
    // mix gives us either first or second vec4 depending if step is 0 or 1
    // nearest_p.xy = coordinates
    // nearest_p.z = store choice -> we want to know if we chose a or b
    // nearest_p.w = dist to point
    highp vec4 nearest_p = mix(vec4(geom_data.a, -1, dist_a), vec4(geom_data.b, 1, dist_b), step(dist_a,dist_b));
    n_line_end *= nearest_p.z; // invert normal direction depending on which vertice we chose

    // we are evaluating inside the line if dist to point is negative
    bool inside_line = nearest_p.w < 0.0;

    // determine if the nearest point is a middle segment of a longer line or if we are at the end
    bool end_cap = bool(mix(float(geom_data.additional_info.y), float(geom_data.additional_info.z), step(0.0, nearest_p.z)));

    // logic table
    // butt  +  end_cap + !inside = butt
    // butt  +  end_cap +  inside = butt
    // butt  + !end_cap + !inside = round
    // butt  + !end_cap +  inside = square
    // round +  end_cap + !inside = round
    // round +  end_cap +  inside = square
    // round + !end_cap + !inside = round
    // round + !end_cap +  inside = square

    // 0u butt
    // 1u round
    // 2u square (round but inside line -> approximate round cap by square)
    //
    // (!round_line_caps && end_cap) -> we want a butt ending -> we negate it and multiply with previous type to force it to be 0u
    // for round vs square -> if we are inside we approximate with square (2u) otherwise we want it round (1u)
    lowp uint ending_type = (1u + uint(inside_line)) * uint(!(!round_line_caps && end_cap));

    // ending_type = 2u; // force specific cap (for debug)

    if(ending_type == 0u)
        return create_halfspace_with_normal(nearest_p.xy, n_line_end);
    else if(ending_type == 1u)
        return create_round_cap(uv, nearest_p.xy, line_width);
    else
        return create_halfspace_with_normal(nearest_p.xy + n_line_end * line_width, n_line_end);
}

// if we have dashes we are changing the a and b vertices to the nearest dash
// we also set line_caps in additional_info to true if we are within the line segment and not at the end
// -> this allows us to draw butt endings within dashed lines
void apply_dashes(highp vec2 uv, inout VectorLayerData geom_data, lowp vec2 dash_info)
{
    if(dash_info.x >= 0.99) // no dashes required
        return;

    highp vec2 e = geom_data.b - geom_data.a;

    // we want the squared distance since t calculation needs it
    highp float squared_dist = dot(e,e);

    // value between 0 and 1, depending on where on the line we are (start point geom_data.a)
    highp float t = clamp(dot(uv - geom_data.a, e)/squared_dist, 0.0, 1.0);

    // how many dash_gap pairs can we fit
    highp float amount_dash_gap_pairs = ceil(sqrt(squared_dist)/dash_info.y);

    // normed to [0,1] range how long is dash_gap for this line segment
    // NOTE: every line segment has slightly different dash gap sizes
    highp float dash_gap_pair_size = 1.0 / amount_dash_gap_pairs;

    // which dash_gap_pair index are we on
    highp float dash_size = dash_gap_pair_size * dash_info.x/2.0;
    lowp float dash_gap_index = floor(1.0 + (t + dash_size) / dash_gap_pair_size) - 1.0;

    highp float t0 = max(0.0, dash_gap_index * dash_gap_pair_size - dash_size);
    highp float t1 = min(1.0, dash_gap_index * dash_gap_pair_size + dash_size);

    // calculate new vertices from the dashes
    // important calculate b before a -> since a overrides value
    geom_data.b = geom_data.a + t1*e;
    geom_data.a = geom_data.a + t0*e;

    // force line_cap to true if we are within a dash (not at the line end)
    // if line_cap was set and we are at the end of the line -> we do need to keep the line_cap set from preprocessor
    // apparently |= does not work for bools in glsl?
    geom_data.additional_info.y = geom_data.additional_info.y || (dash_gap_index > 0.0);
    geom_data.additional_info.z = geom_data.additional_info.z || (dash_gap_index < amount_dash_gap_pairs);
}



void halfspace_uv_to_fragspace(inout highp vec3 halfspace, highp mat3x3 matrix)
{
    halfspace = matrix * halfspace;
    halfspace /= length(halfspace.xy);
}

// orders the halfspaces acccording to distance
void order_halfspace_distance(highp vec3 halfspaces[3], out highp int halfspace_order[3])
{
    halfspace_order[0] = 0;
    halfspace_order[1] = 1;
    halfspace_order[2] = 2;

    // compare halfspaces[0] and halfspaces[1]
    int swap = int(halfspaces[halfspace_order[0]].z < halfspaces[halfspace_order[1]].z);
    int temp = halfspace_order[0];
    halfspace_order[0] = swap * halfspace_order[1] + (1 - swap) * temp;
    halfspace_order[1] = swap * temp + (1 - swap) * halfspace_order[1];

    // compare halfspaces[1] and halfspaces[2]
    swap = int(halfspaces[halfspace_order[1]].z < halfspaces[halfspace_order[2]].z);
    temp = halfspace_order[1];
    halfspace_order[1] = swap * halfspace_order[2] + (1 - swap) * temp;
    halfspace_order[2] = swap * temp + (1 - swap) * halfspace_order[2];

    // compare halfspaces[0] and halfspaces[1] again
    swap = int(halfspaces[halfspace_order[0]].z < halfspaces[halfspace_order[1]].z);
    temp = halfspace_order[0];
    halfspace_order[0] = swap * halfspace_order[1] + (1 - swap) * temp;
    halfspace_order[1] = swap * temp + (1 - swap) * halfspace_order[1];
}

// adapted from https://computergraphics.stackexchange.com/a/13665
lowp uint max_index(highp float x, highp float y, highp float z)
{
   return uint((y>z)&&(y>x)) + (uint((z>y)&&(z>x))*2u);
}

highp float halfspace_coverage(highp float kernel_size, highp float distance)
{
    // linear interpolation (currently kernel_size is alway 1 here)
    return clamp((1.0 - (distance * 0.5 + 0.5)), 0.0, 1.0);
    // smooth interpolation
    // return smoothstep(kernel_size,-kernel_size, distance);
}

highp float calculate_coverage(highp vec3 halfspaces[3], highp int halfspace_order[3], highp float kernel_size, bool inner_edge)
{
    highp float d0 = 0.0;
    if(inner_edge)
        d0 = step(halfspaces[halfspace_order[0]].z, 0.0);
    else
        d0 = halfspace_coverage(kernel_size, halfspaces[halfspace_order[0]].z);
    highp float d1 = halfspace_coverage(kernel_size, halfspaces[halfspace_order[1]].z);
    highp float d2 = halfspace_coverage(kernel_size, halfspaces[halfspace_order[2]].z);

    // determine if we need to subtract or multiply remaining two half spaces
    // -> this depends if the normal is orthogonal or not to normal of nearest halfspace
    highp float perpendicular_multiplications = 1.0;
    highp float paralell_subractions = 0.0;

    highp float dot_01 = dot(halfspaces[halfspace_order[0]].xy, halfspaces[halfspace_order[1]].xy);
    highp float dot_02 = dot(halfspaces[halfspace_order[0]].xy, halfspaces[halfspace_order[2]].xy);

    if(abs(dot_01) < 0.5) // perpendicular
    {
        perpendicular_multiplications = d1;
    }
    else
    {
        // TODO is case:
        // - triangle normal of a triangle -> both normals look in same direction
        // correctly handled?
        paralell_subractions = (1.0 - d1);
    }
    // TODO for the second case we only want to use it if we choose the other method as the first case
    // -> min and max are here to prevent this for now but there should be better method where if 1 sets mult 2 checks only for sub and only sets this
    if(abs(dot_02) < 0.5) // perpendicular
    {
        perpendicular_multiplications = min(perpendicular_multiplications, d2);
    }
    else
    {
        // TODO same as above
        paralell_subractions = max(paralell_subractions, (1.0 - d2));
    }


    return (d0 - paralell_subractions) * perpendicular_multiplications; // no flickering
    // return (d0 - paralell_subractions);  // less flickering
    // return d0; // loads of flickering
}

// we only want to compute the coverage of the halfspace with the largest distance
// if it is outside it depends how much outside it is (if it is not a triangle with an inner_edge)
// if it is inside we only care for the closest halfspace distance
// biggest problem with this is, that thin lines will have aliasing artifacts (because it depends where we sample the thin line)
highp float calculate_coverage_simple(highp float d, highp float kernel_size, bool inner_edge)
{
    if(inner_edge)
        return step(d, 0.0);
    else
        return halfspace_coverage(kernel_size, d);
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
    // TODO I think it should not be necessary to clamp the zoom offset anymore
    lowp int zoom_offset_lower = max(int(floor(zoom_offset - 1.0)), -max_offset_levels + 1);
    lowp int zoom_offset_higher = max(int(floor(zoom_offset - 0.0)), -max_offset_levels + 1);

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
    mediump float zoomed_tile_extent_lower = inv_tile_extent * pow(0.5, float(zoom_offset_lower + 1));
    mediump float zoomed_tile_extent_higher = inv_tile_extent * pow(0.5, float(zoom_offset_higher + 1));

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

#if SDF_MODE == 0
void alpha_blend(inout lowp vec4 pixel_color, LayerStyle style, highp uint intersections)
{
    // we store which sample has hit the geometry -> if two geometries hit the same sample we only store one hit
    // this should make multisample anti-aliasing a bit better
    lowp float intersection_percentage = hit_percentage(intersections);

    pixel_color = pixel_color + ((1.0-pixel_color.a) * style.color * intersection_percentage);
}
#else
void alpha_blend(inout lowp vec4 pixel_color, LayerStyle style, highp float intersection_percentage)
{


    pixel_color = pixel_color + ((1.0-pixel_color.a) * style.color * intersection_percentage);
}

// frag space -> origin is fragment (~pixel) center. going 0.5 units to left, right, up or down -> you reached the border of the fragment
// additionally since we are transforming normals we need to inverse and transpose the matrix
highp mat3x3 create_uv2fragspace_normal_matrix(in highp vec3 normal, in highp uint zoom_level, in highp vec3 ws_position, in highp vec2 uv_position)
{
    highp float scale = tile_size(zoom_level);

    highp vec3 x_axis = vec3(1, 0, 0);
    highp vec3 y_axis = vec3(0, 1, 0);
    highp vec3 z_axis = vec3(0, 0, 1);
    highp vec3 scaled_Rx = cross(vec3(0, 1, 0), normal);
    highp vec3 scaled_Ry = cross(normal, vec3(1, 0, 0));
    // highp vec3 scaled_Rz = scale * normal;

    highp mat3x4 uv2world = mat3x4(
                vec4(scaled_Rx, 0.0),
                vec4(scaled_Ry, 0.0),
                vec4(ws_position, 1.0)) * mat3(vec3(scale / dot(scaled_Rx, x_axis), 0, 0), vec3(0, -scale / dot(scaled_Ry, y_axis), 0), vec3(0, 0, 1)) * mat3x3(vec3(1, 0, 0), vec3(0, 1, 0), vec3(-uv_position, 1));

    // glsl is column major, but we want to remove the 3rd row,
    // so we transpose, remove the 3rd col and transpose back.
    // can be done faster by direct element access
    // even the matmul can be reduced (but probably the compiler does most of that (?))

    highp mat4x3 uv2clipspace_t = transpose(camera.view_proj_matrix * uv2world);

    // scale to [-screen_size,screen_size] -> translate to [0,screen_size] -> translate to [0+current_frag, screensize+current_frag]
    // result -> origin is at the current frag_coord scale is screen_size
    // frag_coord
    highp mat3x3 clip2fragspace =  mat3x3(  vec3(camera.viewport_size.x/2.0,               0,         0),
                                            vec3(0.0,             camera.viewport_size.y/2.0,         0),
                                            vec3(   camera.viewport_size/2.0-gl_FragCoord.xy,         1)
                                         );

    return inverse(transpose(clip2fragspace * transpose(mat3(uv2clipspace_t[0], uv2clipspace_t[1], uv2clipspace_t[3]))));
}

#endif

#if SDF_MODE == 0
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

     return false;

}
#else
bool draw_layer(inout lowp vec4 pixel_color, inout highp float intersection_percentage, inout LayerStyle style, highp vec2 screenspace_coord, highp vec2 uv, highp uint i, DrawMeta meta)
{
    highp uvec2 raw_geom_data = fetch_raw_geometry_data(meta.sampler_buffer_index, i, meta.texture_layer);
    highp uint style_index = unpack_style_index(raw_geom_data) + meta.tile_zoom;

    if (style_index != style.index) {
        // we changed style -> blend previous style, reset layer infos and parse the new style

        alpha_blend(pixel_color, style, intersection_percentage);
        intersection_percentage = 0.0;
        if (pixel_color.a > full_threshold) {
            return true;
        }

        parse_style(style, style_index, meta.zoom_offset, meta.zoom_blend, meta.ortho_color, meta.cos_smoothing_factor, is_polygon(raw_geom_data)); // mix floating zoom levels; for polygons, mul poly color with surface shading texture color
    }

    { // screenspace transformation nehab
        VectorLayerData geom_data = unpack_data(raw_geom_data, meta.grid_cell_float);

        // three half spaces
        // half space definition: xy=normalized normal; z=distance between origin and nearest point on line
        // distance negative implies that we are within the shape
        highp vec3 halfspaces[3];
        halfspaces[0] = create_halfspace(geom_data.a, geom_data.b);

        bool inner_edge[3];


        if(!geom_data.is_polygon)
        {
            inner_edge[0] = false;
            inner_edge[1] = false;
            inner_edge[2] = false;

            // apply_dashes(uv, geom_data, style.dash_info);

            halfspaces[2] = create_line_segment_end_halfspace(uv, geom_data, halfspaces[0].xy, style.line_width, style.round_line_caps);

            halfspaces[1] = -halfspaces[0];
            halfspaces[0].z -= style.line_width;
            halfspaces[1].z -= style.line_width;


            // // TODO calculate dashes and change a and b location


            // {// DEBUG value testing
            //     halfspace_uv_to_fragspace(halfspaces[0], meta.uv2fragspace_normal_matrix);
            //     halfspace_uv_to_fragspace(halfspaces[1], meta.uv2fragspace_normal_matrix);
            //     halfspace_uv_to_fragspace(halfspaces[2], meta.uv2fragspace_normal_matrix);

            //     // determine if if the current fragment is within the negative or the positive side of all halfspaces
            //     // if all negative -> we are within the line
            //     highp float test0 = step(halfspaces[0].z,0.0);
            //     highp float test1 = step(halfspaces[1].z,0.0);
            //     highp float test2 = step(halfspaces[2].z,0.0);
            //     highp float test = test0*test1*test2;

            //     pixel_color = vec4(vec3(test),1.0);
            //     return true;
            // }
        }
        else
        {
            halfspaces[0] = create_halfspace(geom_data.b, geom_data.a);
            halfspaces[1] = create_halfspace(geom_data.c, geom_data.b);
            halfspaces[2] = create_halfspace(geom_data.a, geom_data.c);

            inner_edge[0] = geom_data.additional_info.x;
            inner_edge[1] = geom_data.additional_info.y;
            inner_edge[2] = geom_data.additional_info.z;

        }


        // convert half spaces into fragspace
        halfspace_uv_to_fragspace(halfspaces[0], meta.uv2fragspace_normal_matrix);
        halfspace_uv_to_fragspace(halfspaces[1], meta.uv2fragspace_normal_matrix);
        halfspace_uv_to_fragspace(halfspaces[2], meta.uv2fragspace_normal_matrix);

#if COVERAGE_SIMPLE == 0
        // order halfspaces -> index 0 is always the nearest
        highp int halfspace_order[3];
        order_halfspace_distance(halfspaces, halfspace_order);

        highp float d = calculate_coverage(halfspaces, halfspace_order, 1.0, inner_edge[halfspace_order[0]]);
#else
        lowp uint biggest_halfspace_distance_index = max_index(halfspaces[0].z,halfspaces[1].z,halfspaces[2].z);
        highp float d = calculate_coverage_simple(halfspaces[biggest_halfspace_distance_index].z, 1.0, inner_edge[biggest_halfspace_distance_index]);
#endif

        intersection_percentage = max(d, intersection_percentage);

    }

    return false;
}
#endif

