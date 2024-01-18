/*****************************************************************************
 * Alpine Terrain Renderer
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
uniform sampler2D font_sampler;
uniform sampler2D icon_sampler;

in highp vec2 texcoords;
flat in highp float opacity;

layout (location = 0) out lowp vec4 out_Color;

lowp vec3 fontColor = vec3(0.1f);
lowp vec3 outlineColor = vec3(0.9f);

void main() {

    if(texcoords.x < 2.0f)
    {
        mediump float outline_mask = texture(font_sampler, texcoords).g;
        mediump float font_mask = texture(font_sampler, texcoords).r;

        out_Color = vec4(mix(outlineColor, fontColor, font_mask), outline_mask);
    }
    else
    {
        out_Color = texture(icon_sampler, texcoords-vec2(10.0f,10.0f));
    }

    // apply opacity calculated in vertex shader (opacity from occlusion and distance)
    out_Color.w = out_Color.w * opacity;
}