/*****************************************************************************
 * Alpine Terrain Renderer
 * Copyright (C) 2022 Adam Celarek
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

#include "camera_config.glsl"

const float nearLabel = 100.0;
const float farLabel = 500000.0;


vec2 offset_mask[4] = vec2[4](vec2(0,0), vec2(0,1), vec2(1,1), vec2(1,0));

layout (location = 0) in vec4 pos;
layout (location = 1) in vec4 vtexcoords;
layout (location = 2) in vec3 label_position;

out highp vec2 texcoords;

uniform highp mat4 inv_view_rot;

uniform bool label_dist_scaling;

void main() {
    float dist = length(label_position - camera.position.xyz);
    // remove distance scaling of labels
    float scale = dist / (camera.viewport_size.y * 0.5 * camera.distance_scaling_factor);

    // apply "soft" distance scaling depending on near/far label values
    if(label_dist_scaling)
    {
        float dist_scale = 1.0 - ((dist - nearLabel) / (farLabel - nearLabel)) * 0.4f;
        scale *= (dist_scale * dist_scale);
    }

    // remove rotation from position -> since we want to always face the camera
    // and apply the scaling

    vec4 rotationless_pos = (inv_view_rot * vec4((pos.xy + pos.zw*offset_mask[gl_VertexID]) * scale,0.0,1.0));
    rotationless_pos /= rotationless_pos.w;

    // apply camera matrix and position the label depending on world/camera position
    gl_Position = camera.view_proj_matrix * (vec4((label_position - camera.position.xyz) + rotationless_pos.xyz, 1.0));

    // pass through
    texcoords = vtexcoords.xy + vtexcoords.zw * offset_mask[gl_VertexID];
}