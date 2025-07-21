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

struct Style_Data
{
    lowp vec4 fill_color;
    lowp vec4 outline_color;
    lowp float outline_width;
    lowp vec2 outline_dash;
};

struct Layer_Style
{
    highp uint style_index;
    lowp float layer_alpha;
    Style_Data lower_zoom_style;
    Style_Data higher_zoom_style;
};


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

const highp vec2 cell_size = vec2(float(tile_extent)) / vec2(grid_size);
const highp uint layer_mask = ((1u << sampler_offset) - 1u);
const highp uint bit_mask_ones = -1u;

const highp float full_threshold = 0.999;


highp float calculate_falloff(highp float dist, highp float from, highp float to) {
    return clamp(1.0 - (dist - from) / (to - from), 0.0, 1.0);
}

highp vec3 normal_by_fragment_position_interpolation() {
    highp vec3 dFdxPos = dFdx(var_pos_cws);
    highp vec3 dFdyPos = dFdy(var_pos_cws);
    return normalize(cross(dFdxPos, dFdyPos));
}

// https://iquilezles.org/articles/distfunctions2d/
highp float circle(highp vec2 uv, highp vec2 pos, highp float r) // DEBUG
{
    return distance(uv, pos) - r;
}

highp float sd_Line_Triangle( in highp vec2 p, in highp vec2 p0, in highp vec2 p1, in highp vec2 p2, bool triangle )
{

    highp vec2 e0 = p1-p0;
    highp vec2 v0 = p-p0;
    highp vec2 pq0 = v0 - e0*clamp( dot(v0,e0)/dot(e0,e0), 0.0, 1.0 );

    if(triangle)
    {
        highp vec2 e1 = p2 - p1;
        highp vec2 e2 = p0 - p2;
        highp vec2 v1 = p  - p1;
        highp vec2 v2 = p  - p2;
        highp vec2 pq1 = v1 - e1*clamp( dot(v1,e1)/dot(e1,e1), 0.0, 1.0 );
        highp vec2 pq2 = v2 - e2*clamp( dot(v2,e2)/dot(e2,e2), 0.0, 1.0 );
        highp float s = sign( e0.x*e2.y - e0.y*e2.x );
        highp vec2 d0 = vec2(dot(pq0,pq0), s*(v0.x*e0.y-v0.y*e0.x));
        highp vec2 d1 = vec2(dot(pq1,pq1), s*(v1.x*e1.y-v1.y*e1.x));
        highp vec2 d2 = vec2(dot(pq2,pq2), s*(v2.x*e2.y-v2.y*e2.x));
        highp vec2 d = min(min(d0,d1),d2);

        return -sqrt(d.x)*sign(d.y); // TODO sqrt and length is basically the same instruction -> get them out of the if
    }
    else{
        return length( pq0 );
    }

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


VectorLayerData vertex_sample(lowp uint sampler_index, highp uint index, highp uint texture_layer)
{

    mediump ivec3 dict_px = ivec3(int(index & ((64u<<sampler_index)-1u)), int(index >> (6u+sampler_index)), texture_layer);

    // return unpack_data(texelFetch(geometry_buffer_sampler[sampler_index], dict_px, 0).rg);

    highp uvec2 data;

    switch (sampler_index) {
        case 0u:
            data = texelFetch(geometry_buffer_sampler_0, dict_px, 0).rg;
            break;
        case 1u:
            data = texelFetch(geometry_buffer_sampler_1, dict_px, 0).rg;
            break;
        case 2u:
            data = texelFetch(geometry_buffer_sampler_2, dict_px, 0).rg;
            break;
        default:
            data = texelFetch(geometry_buffer_sampler_3, dict_px, 0).rg;
    }

    return unpack_data(data);
}

mediump float float_zoom_interpolation()
{
    // 3d
    highp float dist_camera = length(var_pos_cws.xyz);
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

Style_Data parse_style(highp uint style_index) {
    Style_Data style;

    highp uvec4 style_data = texelFetch(styles_sampler, ivec2(to_dict_pixel_64(style_index)), 0);
    style.fill_color = vec4((style_data.r & 4278190080u) >> 24, (style_data.r & 16711680u) >> 16, (style_data.r & 65280u) >> 8, style_data.r & 255u) / vec4(255.0f);
    style.outline_color = vec4((style_data.g & 4278190080u) >> 24, (style_data.g & 16711680u) >> 16, (style_data.g & 65280u) >> 8, style_data.g & 255u) / vec4(255.0f);

    style.outline_width = float(style_data.b) / float(style_precision);
    // lowp vec2 outline_dash; // TODO

    return style;
}

/**
  * zoom_blend: value between 0-1, determines how much of the lower and how much of the higher style should be used
  */
void draw_layer(inout Layer_Style layer_style, inout lowp vec4 pixel_color, highp float zoom_blend)
{
    if(layer_style.style_index == -1u)
        return; // we currently have an invalid style -> do not draw anything

    // mix the previous layer color information with output
    // calculate color of current layer by blending the styles of lower and higher zoom step (depending on zoom_blend factor)
    lowp vec4 col = mix(layer_style.lower_zoom_style.fill_color, layer_style.higher_zoom_style.fill_color, zoom_blend);

    // merge current layer color with previous pixel color
    pixel_color = pixel_color + ((1.0-pixel_color.a)*col) * layer_style.layer_alpha;
}

bool check_layer(VectorLayerData geometry_data, inout Layer_Style layer_style, inout lowp vec4 pixel_color, highp float float_zoom_offset)
{
    // we need to make sure that a layer with < 1 opacity does only fill the correct amount of opacity
    // we therefore fill colors per layerstyle

    if(layer_style.style_index == geometry_data.style_index)
    { // new geometry is on same layer as previous one

        if(layer_style.layer_alpha >= full_threshold)
        {
            // we already fully filled the current layer -> go to the next layer
            return true;
        }
    }
    else
    {// we encountered a new layer

        if(layer_style.layer_alpha < full_threshold)
        {
            // we havent drawn the previous layer yet, but encountered a new layer
            // draw the previous layer
            draw_layer(layer_style, pixel_color, fract(float_zoom_offset));
        }

        lowp int zoom_offset_lower = max(int(floor(float_zoom_offset-1.0)), -mipmap_levels+1); // calculate an integer zoom offset for lower style and clamp

        // get and store new style info
        layer_style.style_index = geometry_data.style_index;
        layer_style.lower_zoom_style = parse_style(uint(int(geometry_data.style_index) + min(zoom_offset_lower,0)));

        // the outline_width is saved as tile_extent dependent
        // by dividing by tile_extent we get the width we want to draw
        // by further dividing the tile_extent by 2^zoom_offset, we reduce the tile_extent and increase the line width.
        // we have to increase the zoom_offset by one in order to use the same size as in the preprocessor
        mediump float zoomed_tile_extent = float(tile_extent) * pow(2.0, float(zoom_offset_lower+1));

        layer_style.lower_zoom_style.outline_width /= zoomed_tile_extent;


        // lowp int zoom_offset_higher = clamp(int(floor(float_zoom_offset-0.0)), -mipmap_levels+1, 0); // calculate an integer zoom offset for higher style and clamp
        lowp int zoom_offset_higher = max(int(floor(float_zoom_offset-0.0)), -mipmap_levels+1); // calculate an integer zoom offset for higher style and clamp

        layer_style.higher_zoom_style = parse_style(uint(int(geometry_data.style_index)+min(zoom_offset_higher,0)));

        mediump float zoomed_tile_extent_higher = float(tile_extent) * pow(2.0, float(zoom_offset_higher+1));
        layer_style.higher_zoom_style.outline_width /= zoomed_tile_extent_higher;



        // how much alpha per layer we accumulate
        layer_style.layer_alpha = 0.0;
    }

    return false;

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
    lowp vec3 debug_triangle_lines = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp vec3 debug_index_buffer_start = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp vec3 debug_index_buffer_size = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp vec3 debug_texture_layer = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp int debug_draw_calls = 0; // DEBUG
    lowp float debug_poly_over_line = 0.0; // DEBUG



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
    lowp vec3 ortho_color = texture(texture_sampler, vec3(ortho_uv, texture_layer_f)).rgb;
    ortho_color = mix(ortho_color, conf.material_color.rgb, conf.material_color.a);

    // VECTOR color
    decrease_zoom_level_until(tile_id, uv, texelFetch(instanced_texture_zoom_sampler_vector, ivec2(instance_id, 0), 0).x);


    highp float float_zoom = float_zoom_interpolation();
    // float_zoom = tile_id.z;     // DEBUG
    lowp uint mipmap_level = uint(clamp(int(ceil(float(tile_id.z)-float_zoom)), 0,mipmap_levels-1));

    highp uvec2 texture_layer = texelFetch(instanced_texture_array_index_sampler_vector, ivec2(instance_id, mipmap_level), 0).xy;

    decrease_zoom_level_until(tile_id, uv, tile_id.z - mipmap_level);
    highp float zoom_offset = float_zoom-float(tile_id.z);


    highp float aa_half_radius = calculate_aa_half_radius(uv);

    highp float line_influence = 0.0;

    if(texture_layer.x != bit_mask_ones && texture_layer.y != bit_mask_ones) // check for valid data
    {
        // acceleration_grid_sampler contains the offset and the number of triangles of the current grid cell
        highp vec2 grid_lookup = grid_size*uv;
        offset_size = to_offset_size(texelFetch(acceleration_grid_sampler, ivec3(int(grid_lookup.x), int(grid_lookup.y), texture_layer.x & layer_mask),0).r);

        // using the grid data we now want to traverse all triangles referenced in grid cell and draw them.
        if(offset_size.y != uint(0)) // only if we have data here
        {
            // lowp vec3 raw_grid = vec3(float(offset_size.y),0,0);// DEBUG
            lowp vec3 raw_grid = vec3(1,0,0);// DEBUG
            lowp ivec2 grid_cell = ivec2(int(grid_lookup.x), int(grid_lookup.y));
            lowp vec3 cells = color_from_id_hash(uint(grid_cell.x ^ grid_cell.y)); // DEBUG



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

            Layer_Style layer_style;
            layer_style.style_index = -1u;
            layer_style.lower_zoom_style = Style_Data(vec4(0.0), vec4(0.0), 0.0, vec2(0.0));
            layer_style.higher_zoom_style = layer_style.lower_zoom_style;
            layer_style.layer_alpha = 0.0;

            int polygon_was_drawn = 0;

            lowp float geometry_influence = 0.0; // how much is the current line/polygon visible
            for(highp uint i = offset_size.x; i < offset_size.x + offset_size.y; i++)
                // for(highp uint i = offset_size.x; i < offset_size.x + min(6u,offset_size.y); i++) // show only x layers
            // for(highp uint i = offset_size.x+ offset_size.y; i --> offset_size.x ; ) // reverse traversal
            {

                debug_draw_calls = debug_draw_calls + 1;
                VectorLayerData geometry_data = vertex_sample(sampler_buffer_index,i, texture_layer.y);


                bool check_next_geometry = check_layer(geometry_data, layer_style, pixel_color, zoom_offset);
                if(check_next_geometry) // current layer is full -> we want to move on to a geometry on a different layer
                    continue;


                highp float d = 0.0;

                lowp float tile_scale = float(scale_lines) * (1.0-float(geometry_data.is_polygon)) + float(scale_polygons) * float(geometry_data.is_polygon);
                highp vec2 cell_offset = vec2(grid_cell) * cell_size * tile_scale;

                highp vec2 v0 = (vec2(geometry_data.a) + cell_offset) / vec2(float(tile_extent) * tile_scale);
                highp vec2 v1 = (vec2(geometry_data.b) + cell_offset) / vec2(float(tile_extent) * tile_scale);
                highp vec2 v2 = (vec2(geometry_data.c) + cell_offset) / vec2(float(tile_extent) * tile_scale);


                lowp float thickness_lower = layer_style.lower_zoom_style.outline_width;
                lowp float thickness_higher = layer_style.higher_zoom_style.outline_width;
                lowp float thickness = mix(thickness_lower, thickness_higher, fract(float_zoom));


                d = sd_Line_Triangle(uv, v0, v1, v2, geometry_data.is_polygon) - thickness;


                highp float geometry_influence =  1.0 - smoothstep(-aa_half_radius, aa_half_radius, d);

                // polygon does not influence pixel at all -> we do not draw it
                if(geometry_influence <= 0.0)
                    continue;

                line_influence += (1.0-float(geometry_data.is_polygon)) * geometry_influence;

                if(geometry_data.is_polygon)
                  polygon_was_drawn += 1;
                else if(polygon_was_drawn > 0)
                  debug_poly_over_line = float(polygon_was_drawn) / 1.0;


                { // DEBUG -> triangle lines
                    highp float t_line = 0.0;
                    if(d <= 0.001 && d >= -0.001)
                        t_line = 0.2;
                    debug_triangle_lines += vec3(0.0f, t_line, 0.0f);
                }




                // set layer_alpha from the current geometry
                layer_style.layer_alpha = min(1.0, layer_style.layer_alpha + geometry_influence);

                // we do not need to check any other geometry on this layer -> layer is already fully filled
                // -> draw the current layer to see if we can exit the whole loop in next if
                if(layer_style.layer_alpha >= full_threshold)
                    draw_layer(layer_style, pixel_color, fract(zoom_offset));

                // current pixel is fully filled -> no need to check anything further
                if(pixel_color.a >= full_threshold)
                    break;

            }


            if(pixel_color.a < 1.0  && layer_style.layer_alpha < full_threshold)
            {
                // mix the last layer we parsed
                draw_layer(layer_style, pixel_color, fract(float_zoom));
            }
        }
    }

    // mix polygon color with background
    lowp float vector_surface_mix = 0.0;
    if(display_mode == 0)
    {
        // mixed surface/vectortile
        vector_surface_mix = 0.5 + min(1.0, line_influence) * 0.5;
    }
    else if(display_mode == 1)
    {
        // only vector tile
        vector_surface_mix = 1.0;
    }
    // else // -> only surface/ortho

    highp vec3 vector_color = vec3(pixel_color.rgb + ((1.0-pixel_color.a)*background_color));
    // texout_albedo = mix(ortho_color, vector_color, vector_surface_mix);
    texout_albedo = mix(ortho_color * vector_color, vector_color, min(1.0, line_influence));

    if (conf.overlay_mode > 199u && conf.overlay_mode < 300u) {
        lowp vec3 zoom_debug_color =  color_from_id_hash(uint(float_zoom));

        lowp vec3 overlay_color = vec3(0.0);
        switch(conf.overlay_mode) {
            case 200u: overlay_color = vec3(uv, 0.0);break;
            case 201u: overlay_color = debug_cacscade_layer;break;
            case 202u: overlay_color = debug_cell_size;break;
            case 203u: overlay_color = debug_triangle_lines;break;
            case 204u: overlay_color = debug_index_buffer_start;break;
            case 205u: overlay_color = debug_index_buffer_size;break;
            case 206u: overlay_color = debug_texture_layer;break;
            case 207u: overlay_color = vec3(pixel_color.a, 0.0, 0.0);break;
            case 211u: overlay_color = vec3(debug_poly_over_line, 0.0, 0.0);break;
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
