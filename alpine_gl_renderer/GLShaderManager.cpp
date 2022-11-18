/*****************************************************************************
 * Alpine Terrain Builder
 * Copyright (C) 2022 alpinemaps.org
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

#include "GLShaderManager.h"

#include <QOpenGLContext>

#include "ShaderProgram.h"

static const char* const tileVertexShaderSource = R"(
  layout(location = 0) in highp float height;
  uniform highp mat4 matrix;
  uniform highp vec3 camera_position;
  uniform highp vec4 bounds[32];
  uniform int n_edge_vertices;
  out lowp vec2 uv;
  out highp vec3 camera_rel_pos;

  void main() {
    int row = gl_VertexID / n_edge_vertices;
    int col = gl_VertexID - (row * n_edge_vertices);
    int geometry_id = 0;
    float tile_width = (bounds[geometry_id].z - bounds[geometry_id].x) / float(n_edge_vertices - 1);
    float tile_height = (bounds[geometry_id].w - bounds[geometry_id].y) / float(n_edge_vertices - 1);

    camera_rel_pos = vec3(float(col) * tile_width + bounds[geometry_id].x,
                          float(n_edge_vertices - row - 1) * tile_width + bounds[geometry_id].y,
                          height * 65536.0 * 0.125 - camera_position.z);
    uv = vec2(float(col) / float(n_edge_vertices - 1), float(row) / float(n_edge_vertices - 1));
    gl_Position = matrix * vec4(camera_rel_pos, 1);
  })";

static const char* const tileFragmentShaderSource = R"(
  uniform highp vec3 camera_position;
  uniform sampler2D texture_sampler;
  in lowp vec2 uv;
  in highp vec3 camera_rel_pos;
  out lowp vec4 out_Color;
  void main() {
     vec3 origin = vec3(camera_position);
     vec4 ortho = texture(texture_sampler, uv);
     float dist = length(camera_rel_pos) / 1000;
     out_Color = vec4(dist, dist, dist, 1.0);
     //gl_FragDepth = gl_FragCoord.z;
  })";


//#include "UnityCG.cginc"
//#include "../Includes/Math.cginc"

//struct appdata {
//    float4 vertex : POSITION;
//    float4 uv : TEXCOORD0;
//};

//struct v2f {
//    float4 pos : SV_POSITION;
//    float2 uv : TEXCOORD0;
//    float3 viewVector : TEXCOORD1;
//};

//v2f vert (appdata v) {
//    v2f output;
//    output.pos = UnityObjectToClipPos(v.vertex);
//    output.uv = v.uv;
//    // Camera space matches OpenGL convention where cam forward is -z. In unity forward is positive z.
//    // (https://docs.unity3d.com/ScriptReference/Camera-cameraToWorldMatrix.html)
//    float3 viewVector = mul(unity_CameraInvProjection, float4(v.uv.xy * 2 - 1, 0, -1));
//    output.viewVector = mul(unity_CameraToWorld, float4(viewVector,0));
//    return output;
//}

//float2 squareUV(float2 uv) {
//    float width = _ScreenParams.x;
//    float height =_ScreenParams.y;
//    //float minDim = min(width, height);
//    float scale = 1000;
//    float x = uv.x * width;
//    float y = uv.y * height;
//    return float2 (x/scale, y/scale);
//}

//sampler2D _BlueNoise;
//sampler2D _MainTex;
//sampler2D _BakedOpticalDepth;
//sampler2D _CameraDepthTexture;
//float4 params;

//float3 dirToSun;

//float3 planetCentre;
//float atmosphereRadius;
//float oceanRadius;
//float planetRadius;

//// Paramaters
//int numInScatteringPoints;
//int numOpticalDepthPoints;
//float intensity;
//float4 scatteringCoefficients;
//float ditherStrength;
//float ditherScale;
//float densityFalloff;


//float densityAtPoint(float3 densitySamplePoint) {
//    float heightAboveSurface = length(densitySamplePoint - planetCentre) - planetRadius;
//    float height01 = heightAboveSurface / (atmosphereRadius - planetRadius);
//    float localDensity = exp(-height01 * densityFalloff) * (1 - height01);
//    return localDensity;
//}

//float opticalDepth(float3 rayOrigin, float3 rayDir, float rayLength) {
//    float3 densitySamplePoint = rayOrigin;
//    float stepSize = rayLength / (numOpticalDepthPoints - 1);
//    float opticalDepth = 0;

//    for (int i = 0; i < numOpticalDepthPoints; i ++) {
//        float localDensity = densityAtPoint(densitySamplePoint);
//        opticalDepth += localDensity * stepSize;
//        densitySamplePoint += rayDir * stepSize;
//    }
//    return opticalDepth;
//}

//float opticalDepthBaked(float3 rayOrigin, float3 rayDir) {
//    float height = length(rayOrigin - planetCentre) - planetRadius;
//    float height01 = saturate(height / (atmosphereRadius - planetRadius));

//    float uvX = 1 - (dot(normalize(rayOrigin - planetCentre), rayDir) * .5 + .5);
//    return tex2Dlod(_BakedOpticalDepth, float4(uvX, height01,0,0));
//}

//float opticalDepthBaked2(float3 rayOrigin, float3 rayDir, float rayLength) {
//    float3 endPoint = rayOrigin + rayDir * rayLength;
//    float d = dot(rayDir, normalize(rayOrigin-planetCentre));
//    float opticalDepth = 0;

//    const float blendStrength = 1.5;
//    float w = saturate(d * blendStrength + .5);

//    float d1 = opticalDepthBaked(rayOrigin, rayDir) - opticalDepthBaked(endPoint, rayDir);
//    float d2 = opticalDepthBaked(endPoint, -rayDir) - opticalDepthBaked(rayOrigin, -rayDir);

//    opticalDepth = lerp(d2, d1, w);
//    return opticalDepth;
//}

//float3 calculateLight(float3 rayOrigin, float3 rayDir, float rayLength, float3 originalCol, float2 uv) {
//    float blueNoise = tex2Dlod(_BlueNoise, float4(squareUV(uv) * ditherScale,0,0));
//    blueNoise = (blueNoise - 0.5) * ditherStrength;

//    float3 inScatterPoint = rayOrigin;
//    float stepSize = rayLength / (numInScatteringPoints - 1);
//    float3 inScatteredLight = 0;
//    float viewRayOpticalDepth = 0;

//    for (int i = 0; i < numInScatteringPoints; i ++) {
//        float sunRayLength = raySphere(planetCentre, atmosphereRadius, inScatterPoint, dirToSun).y;
//        float sunRayOpticalDepth = opticalDepthBaked(inScatterPoint + dirToSun * ditherStrength, dirToSun);
//        float localDensity = densityAtPoint(inScatterPoint);
//        viewRayOpticalDepth = opticalDepthBaked2(rayOrigin, rayDir, stepSize * i);
//        float3 transmittance = exp(-(sunRayOpticalDepth + viewRayOpticalDepth) * scatteringCoefficients);

//        inScatteredLight += localDensity * transmittance;
//        inScatterPoint += rayDir * stepSize;
//    }
//    inScatteredLight *= scatteringCoefficients * intensity * stepSize / planetRadius;
//    inScatteredLight += blueNoise * 0.01;

//         // Attenuate brightness of original col (i.e light reflected from planet surfaces)
//         // This is a hacky mess, TODO: figure out a proper way to do this
//    const float brightnessAdaptionStrength = 0.15;
//    const float reflectedLightOutScatterStrength = 3;
//    float brightnessAdaption = dot (inScatteredLight,1) * brightnessAdaptionStrength;
//    float brightnessSum = viewRayOpticalDepth * intensity * reflectedLightOutScatterStrength + brightnessAdaption;
//    float reflectedLightStrength = exp(-brightnessSum);
//    float hdrStrength = saturate(dot(originalCol,1)/3-1);
//    reflectedLightStrength = lerp(reflectedLightStrength, 1, hdrStrength);
//    float3 reflectedLight = originalCol * reflectedLightStrength;

//    float3 finalCol = reflectedLight + inScatteredLight;


//    return finalCol;
//}


//float4 frag (v2f i) : SV_Target
//{
//    float4 originalCol = tex2D(_MainTex, i.uv);
//    float sceneDepthNonLinear = SAMPLE_DEPTH_TEXTURE(_CameraDepthTexture, i.uv);
//    float sceneDepth = LinearEyeDepth(sceneDepthNonLinear) * length(i.viewVector);

//    float3 rayOrigin = _WorldSpaceCameraPos;
//    float3 rayDir = normalize(i.viewVector);

//    float dstToOcean = raySphere(planetCentre, oceanRadius, rayOrigin, rayDir);
//    float dstToSurface = min(sceneDepth, dstToOcean);

//    float2 hitInfo = raySphere(planetCentre, atmosphereRadius, rayOrigin, rayDir);
//    float dstToAtmosphere = hitInfo.x;
//    float dstThroughAtmosphere = min(hitInfo.y, dstToSurface - dstToAtmosphere);

//    if (dstThroughAtmosphere > 0) {
//        const float epsilon = 0.0001;
//        float3 pointInAtmosphere = rayOrigin + rayDir * (dstToAtmosphere + epsilon);
//        float3 light = calculateLight(pointInAtmosphere, rayDir, dstThroughAtmosphere - epsilon * 2, originalCol, i.uv);
//        return float4(light, 1);
//    }
//    return originalCol;
//}

static const char* const screenQuadVertexShaderSource = R"(
// https://stackoverflow.com/a/59739538
out highp vec2 texcoords; // texcoords are in the normalized [0,1] range for the viewport-filling quad part of the triangle
void main() {
    vec2 vertices[3]=vec2[3](vec2(-1.0, -1.0),
                             vec2(3.0, -1.0),
                             vec2(-1.0, 3.0));
    gl_Position = vec4(vertices[gl_VertexID], 0.0, 1.0);
    texcoords = 0.5 * gl_Position.xy + vec2(0.5);
})";

static const char* const screenQuadFragmentShaderSource = R"(
  in highp vec2 texcoords;
  uniform sampler2D texture_sampler;
  out lowp vec4 out_Color;
  void main() {
     out_Color = vec4(texcoords.xy, 0.0, 1.0) * 0.1 + texture(texture_sampler, texcoords);
  })";

static const char* const debugVertexShaderSource = R"(
  layout(location = 0) in vec4 a_position;
  uniform highp mat4 matrix;
  void main() {
    gl_Position = matrix * a_position;
  })";

static const char* const debugFragmentShaderSource = R"(
  out lowp vec4 out_Color;
  void main() {
     out_Color = vec4(1.0, 0.0, 0.0, 1.0);
  })";

QByteArray versionedShaderCode(const char* const src)
{
    QByteArray versionedSrc;

    if (QOpenGLContext::currentContext()->isOpenGLES())
        versionedSrc.append(QByteArrayLiteral("#version 300 es\n"));
    else
        versionedSrc.append(QByteArrayLiteral("#version 330\n"));

    versionedSrc.append(src);
    return versionedSrc;
}

GLShaderManager::GLShaderManager()
{
    m_tile_program = std::make_unique<ShaderProgram>(tileVertexShaderSource, tileFragmentShaderSource);
    m_debug_program = std::make_unique<ShaderProgram>(debugVertexShaderSource, debugFragmentShaderSource);
    m_screen_quad_program = std::make_unique<ShaderProgram>(screenQuadVertexShaderSource, screenQuadFragmentShaderSource);
}

GLShaderManager::~GLShaderManager() = default;

void GLShaderManager::bindTileShader()
{
    m_tile_program->bind();
}

void GLShaderManager::bindDebugShader()
{
    m_debug_program->bind();
}

ShaderProgram* GLShaderManager::tileShader() const
{
    return m_tile_program.get();
}

ShaderProgram* GLShaderManager::debugShader() const
{
    return m_debug_program.get();
}

ShaderProgram* GLShaderManager::screen_quad_program() const
{
    return m_screen_quad_program.get();
}

void GLShaderManager::release()
{
    m_tile_program->release();
}
