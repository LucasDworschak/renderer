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
uniform lowp sampler2DArray fallback_texture_array;

uniform highp usampler2D instanced_texture_array_index_sampler;
uniform highp usampler2D instanced_texture_zoom_sampler;

uniform highp usampler2DArray acceleration_grid_sampler;
uniform highp usampler2D instanced_texture_array_index_sampler_vector;
uniform highp usampler2D instanced_texture_zoom_sampler_vector;


uniform highp int max_vector_geometry;

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



highp vec3 normal_by_fragment_position_interpolation() {
    highp vec3 dFdxPos = dFdx(var_pos_cws);
    highp vec3 dFdyPos = dFdy(var_pos_cws);
    return normalize(cross(dFdxPos, dFdyPos));
}


mediump float float_zoom_interpolation()
{
    highp float dist_camera = length(var_pos_cws.xyz);

    const highp float sqrt2 = 1.414213562373095;
    const highp float cEarthCircumference = 40075016.685578486;
    const highp float tile_size = 256.0;

    highp float camera_factors = camera.viewport_size.y * 0.5 * camera.distance_scaling_factor;
    highp float static_factors = camera_factors * sqrt2 * cEarthCircumference / tile_size;
    highp float z = log2(static_factors / dist_camera / camera.error_threshold_px)+1.0;

    return clamp(z, 0.0, float(max_zoom));
}


mediump float calculate_cos_smoothing()
{
#if shallow_angle_signal_frequency == 1
    highp float cos_angle = dot(normalize(var_pos_cws.xyz),vec3(0.,0.,1.));
    highp float dist = smoothstep(2000.0, 500.0, length(var_pos_cws.xy)); // between 0-500m -> no cos smoothing; between 2000m and inf use cos smoothing; inbetween transition

    return mix(0.0, sqrt(sqrt(abs(cos_angle))), smoothstep(0.0,0.15,abs(cos_angle))) * (1.0-dist) + dist;

    // return sqrt(sqrt(abs(cos_angle))) * (1.0-dist) + dist;
#else
    return 1.0;
#endif
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
#if CURTAIN_DEBUG_MODE == 2
    if (is_curtain == 0.0) {
        discard;
    }
#endif

    DrawMeta meta;


    lowp vec3 debug_cacscade_layer = vec3(0,0,0);// DEBUG -> which cascade is being used
    lowp vec3 debug_cell_size = vec3(0,0,0);// DEBUG -> how many triangles per cell
    lowp vec3 debug_index_buffer_start = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp vec3 debug_index_buffer_size = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp vec3 debug_texture_layer = vec3(0.0f, 0.0, 0.0f); // DEBUG
    lowp int debug_draw_calls = 0; // DEBUG



    highp uvec3 tile_id = var_tile_id;
    highp vec2 uv = var_uv;
    highp uvec2 offset_size = uvec2(0u);


    lowp vec4 pixel_color = vec4(0.0f, 0.0, 0.0f, 0.0f);


    /////////////////////////
    // ORTHO color -> note we wrap tile_id in uvec3 and use ortho_uv, to not change those variables for the vector color
    highp vec2 ortho_uv = uv;
    highp uvec3 temp_tile_id = uvec3(tile_id);
    decrease_zoom_level_until(temp_tile_id, ortho_uv, texelFetch(instanced_texture_zoom_sampler, ivec2(instance_id, 0), 0).x);
    highp float texture_layer_f = float(texelFetch(instanced_texture_array_index_sampler, ivec2(instance_id, 0), 0).x);
    meta.ortho_color = vec4(texture(texture_sampler, vec3(ortho_uv, texture_layer_f)).rgb, 1.0);
    meta.ortho_color = mix(meta.ortho_color, conf.material_color, conf.material_color.a);

#if VIEW_MODE == 2
    meta.ortho_color = vec4(1.0);
#endif
    // meta.ortho_color = vec4(1.0);





    /////////////////////////
    // VECTOR color
    decrease_zoom_level_until(tile_id, uv, texelFetch(instanced_texture_zoom_sampler_vector, ivec2(instance_id, 0), 0).x);


    mediump float float_zoom = float_zoom_interpolation();
    // float_zoom = tile_id.z;     // DEBUG
    lowp uint mipmap_level = uint(clamp(int(ceil(float(tile_id.z)-float_zoom)), 0,mipmap_levels-1));

    // calculate uv derivatives before we apply mipmapping -> otherwise we have discontinous areas at border
    // and correct uv derivatives scaling with mipmaplevel
    meta.duvdx = dFdx(uv) * pow(0.5, float(mipmap_level));
    meta.duvdy = dFdy(uv) * pow(0.5, float(mipmap_level));

    highp uvec2 texture_layer = texelFetch(instanced_texture_array_index_sampler_vector, ivec2(instance_id, mipmap_level), 0).xy;

    decrease_zoom_level_until(tile_id, uv, tile_id.z - mipmap_level);

    meta.zoom_offset = float_zoom-float(tile_id.z);
    meta.zoom_blend = fract(meta.zoom_offset);
    meta.tile_zoom = tile_id.z;



    /////////////////////////
    // FALLBACK COLOR
    lowp vec4 fallback_color = vec4(texture(fallback_texture_array, vec3(uv, texture_layer.x & layer_mask),0));

    // TODO not quite sure if 1) we want to mix two fallback colors and 2) if this is the best way to do this
    highp uvec2 texture_layer2 = texelFetch(instanced_texture_array_index_sampler_vector, ivec2(instance_id, mipmap_level+1u), 0).xy;
    highp vec2 fallback2_uv = uv;
    highp uvec3 temp_tile_id2 = uvec3(tile_id);
    decrease_zoom_level_until(temp_tile_id2, fallback2_uv, tile_id.z - 1u);
    lowp vec4 fallback_color2 = vec4(texture(fallback_texture_array, vec3(fallback2_uv, texture_layer2.x & layer_mask),0));




    // acceleration_grid_sampler contains the offset and the number of triangles of the current grid cell
    lowp ivec2 grid_cell = ivec2(grid_size*uv);
    meta.grid_cell_float = vec2(grid_cell);
    offset_size = to_offset_size(texelFetch(acceleration_grid_sampler, ivec3(grid_cell.x, grid_cell.y, texture_layer.x & layer_mask),0).r);



    /////////////////////////
    // anti-alialing
    meta.cos_smoothing_factor = calculate_cos_smoothing();
    meta.cos_smoothing_factor = 1;
    meta.aa_half_radius = calculate_aa_half_radius(uv);

#if SDF_MODE == 0
    calculate_sample_positions(meta.aa_sample_positions, uv, grid_cell);
#endif
#if SDF_MODE == 1
    calculate_sample_multipliers(meta.aa_sample_multipliers, uv, meta.duvdx, meta.duvdy, meta.grid_cell_float);
#endif

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


        for(highp uint i = offset_size.x; i < offset_size.x + min(uint(max_vector_geometry),offset_size.y); i++) // show only x layers
        {
            debug_draw_calls++;
            if(draw_layer(pixel_color, intersection_percentage, style, uv, i, meta))
                break; // pixel is finished -> we can exit the loop early
        }

        // blend the last layer (this only does something if we did not break out of the loop early)
        alpha_blend(pixel_color, style, intersection_percentage);
    }



    // blend with background color and write pixel color to output
#if VIEW_MODE == 0
    texout_albedo = vec3(background_color * meta.ortho_color.rgb);
#else
    texout_albedo = vec3(pixel_color.rgb + ((1.0-pixel_color.a)*background_color * meta.ortho_color.rgb));

    lowp vec4 fallback_mixed_color = mix(fallback_color2, fallback_color, fract(float_zoom));
    // lowp vec3 fallback_ortho_color = mix(fallback_mixed_color.rgb*meta.ortho_color.rgb, fallback_mixed_color.rgb, fallback_mixed_color.a);
    lowp vec3 fallback_ortho_color = mix(fallback_mixed_color.rgb*meta.ortho_color.rgb, fallback_mixed_color.rgb, 0.0);

    texout_albedo = vec3(pixel_color.rgb) + ((1.0-pixel_color.a) * fallback_ortho_color);
    // texout_albedo = mix(fallback_color2, fallback_color, fract(float_zoom)).rgb;
#endif




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
