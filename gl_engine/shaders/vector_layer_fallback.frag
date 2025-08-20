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

#line 30031

uniform lowp sampler2DArray texture_sampler;
uniform lowp sampler2DArray fallback_texture_array;

uniform highp usampler2D instanced_texture_array_index_sampler;
uniform highp usampler2D instanced_texture_zoom_sampler;

uniform highp usampler2DArray acceleration_grid_sampler;
uniform highp usampler2D instanced_texture_array_index_sampler_vector;
uniform highp usampler2D instanced_texture_zoom_sampler_vector;

layout (location = 0) out lowp vec4 texout_albedo;

in highp vec2 var_uv;

uniform highp int instance_id;
uniform highp int tile_zoom;


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

void debug_calculate_cascade_layer(out lowp vec3 debug_cacscade_layer, lowp uint sampler_buffer_index)
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

void debug_calculate_cell_size(out lowp vec3 debug_cell_size, mediump uint offset_size)
{
    if(offset_size < 32u)
        debug_cell_size = vec3(0,1,0); // green
    else if(offset_size < 64u)
        debug_cell_size = vec3(1,1,0); // yellow
    else if(offset_size  < 128u)
        debug_cell_size = vec3(1,0.5,0); // orange
    else if(offset_size < 255u)
        debug_cell_size = vec3(1,0,0); // red
    else
        debug_cell_size = vec3(1,0,1); // purple -> should never happen -> unrecognized index
}


void main() {

    DrawMeta meta;


    lowp vec3 debug_cacscade_layer = vec3(0,0,0);// DEBUG -> which cascade is being used
    lowp vec3 debug_cell_size = vec3(0,0,0);// DEBUG -> how many triangles per cell
    lowp vec3 debug_index_buffer_start = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp vec3 debug_index_buffer_size = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp vec3 debug_texture_layer = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp int debug_draw_calls = 0; // DEBUG



    highp vec2 uv = var_uv;
    highp uvec2 offset_size = uvec2(0u);


    lowp vec4 pixel_color = vec4(0.0f, 0.0, 0.0f, 0.0f);


//     /////////////////////////
//     // ORTHO color -> note we wrap tile_id in uvec3 and use ortho_uv, to not change those variables for the vector color
//     highp vec2 ortho_uv = uv;
//     highp float texture_layer_f = float(texelFetch(instanced_texture_array_index_sampler, ivec2(instance_id, 0), 0).x);
//     meta.ortho_color = vec4(texture(texture_sampler, vec3(ortho_uv, texture_layer_f)).rgb, 1.0);
//     meta.ortho_color = mix(meta.ortho_color, conf.material_color, conf.material_color.a);

// #if VIEW_MODE == 2
//     meta.ortho_color = vec4(1.0);
// #endif
    meta.ortho_color = vec4(1.0);


    /////////////////////////
    // VECTOR color

    mediump float float_zoom = tile_zoom;     // DEBUG

    // calculate uv derivatives before we apply mipmapping -> otherwise we have discontinous areas at border
    meta.duvdx = dFdx(uv);
    meta.duvdy = dFdy(uv);

    highp uvec2 texture_layer = texelFetch(instanced_texture_array_index_sampler_vector, ivec2(instance_id, 0), 0).xy;

    meta.zoom_offset = float_zoom-float(tile_zoom);
    meta.zoom_blend = fract(meta.zoom_offset);
    meta.tile_zoom = uint(tile_zoom);




    // acceleration_grid_sampler contains the offset and the number of triangles of the current grid cell
    lowp ivec2 grid_cell = ivec2(grid_size*uv);
    meta.grid_cell_float = vec2(grid_cell);
    offset_size = to_offset_size(texelFetch(acceleration_grid_sampler, ivec3(grid_cell.x, grid_cell.y, texture_layer.x & layer_mask),0).r);



    /////////////////////////
    // anti-alialing
    meta.cos_smoothing_factor = 1;
    meta.aa_half_radius = calculate_aa_half_radius(uv);

#if SDF_MODE == 0
    calculate_sample_positions(meta.aa_sample_positions, uv, grid_cell);
#endif
#if SDF_MODE == 1
    calculate_sample_multipliers(meta.aa_sample_multipliers, uv, meta.duvdx, meta.duvdy, meta.grid_cell_float);
#endif

    lowp float line_percentage = 0.0;


    // using the grid data we now want to traverse all triangles referenced in grid cell and draw them.
    if(offset_size.y != uint(0)) // only if we have data here
    {
        // get the buffer index and extract the correct texture_layer.y
        lowp uint sampler_buffer_index = texture_layer.y >> sampler_offset;
        texture_layer.y = texture_layer.y & layer_mask;

        meta.texture_layer = texture_layer.y;
        meta.sampler_buffer_index = sampler_buffer_index;

        { // DEBUG
            debug_texture_layer = color_from_id_hash(texture_layer.y);
            debug_calculate_cascade_layer(debug_cacscade_layer, sampler_buffer_index);
            debug_calculate_cell_size(debug_cell_size, offset_size.y);
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

        for(highp uint i = offset_size.x; i < offset_size.x + min(255u,offset_size.y); i++) // show only x layers
        {
            debug_draw_calls++;
            if(draw_layer(line_percentage, pixel_color, intersection_percentage, style, uv, i, meta))
                break; // pixel is finished -> we can exit the loop early
        }

        // blend the last layer (this only does something if we did not break out of the loop early)
        alpha_blend(pixel_color, style, intersection_percentage);
    }



    // blend with background color and write pixel color to output
#if VIEW_MODE == 0
    texout_albedo = vec4(background_color * meta.ortho_color.rgb, 1.0);
#else
    // the alpha value is 1 if we encountered a line (with at least 0.2 percentage summed up) or 0 if only polygons have been encountered
    // -> important for ortho color mixing
    texout_albedo = vec4(pixel_color.rgb + ((1.0-pixel_color.a)*background_color), 1.0);// step(0.2,line_percentage));

     // texout_albedo = vec4(debug_index_buffer_size, 1.0);

    // texout_albedo = vec3(pixel_color.rgb) + ((1.0-pixel_color.a)*mix(fallback_color2, fallback_color, fract(float_zoom)).rgb);
#endif


}
