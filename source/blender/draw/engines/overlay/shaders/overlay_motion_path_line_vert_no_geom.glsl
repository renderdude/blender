/* SPDX-FileCopyrightText: 2022-2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#pragma USE_SSBO_VERTEX_FETCH(TriangleList, 6)

#include "common_view_clipping_lib.glsl"
#include "common_view_lib.glsl"

#define frameCurrent mpathLineSettings.x
#define frameStart mpathLineSettings.y
#define frameEnd mpathLineSettings.z
#define cacheStart mpathLineSettings.w

/* Project to screen space. */
vec2 proj(vec4 pos)
{
  return (0.5 * (pos.xy / pos.w) + 0.5) * sizeViewport.xy;
}

#define SET_INTENSITY(A, B, C, min, max) \
  (((1.0 - (float(C - B) / float(C - A))) * (max - min)) + min)

vec2 compute_dir(vec2 v0, vec2 v1)
{
  vec2 dir = normalize(v1 - v0 + 1e-8);
  dir = vec2(-dir.y, dir.x);
  return dir;
}

void do_vertex_shader(vec4 pos, int vertex_id, out vec2 out_sspos, out vec4 out_finalcolor)
{
  out_sspos = proj(pos);
  out_finalcolor = vec4(0.0);

  int frame = vertex_id + cacheStart;
  float intensity; /* how faint */
  vec3 blend_base = (abs(frame - frameCurrent) == 0) ?
                        colorCurrentFrame.rgb :
                        colorBackground.rgb; /* "bleed" CFRAME color to ease color blending. */
  bool use_custom_color = customColorPre.x >= 0.0;

  if (frame < frameCurrent) {
    if (use_custom_color) {
      out_finalcolor.rgb = customColorPre;
    }
    else {
      if (selected) {
        intensity = SET_INTENSITY(frameStart, frame, frameCurrent, 0.25, 0.75);
      }
      else {
        intensity = SET_INTENSITY(frameStart, frame, frameCurrent, 0.68, 0.92);
      }
      out_finalcolor.rgb = mix(colorWire.rgb, blend_base, intensity);
    }
  }
  else if (frame > frameCurrent) {
    if (use_custom_color) {
      out_finalcolor.rgb = customColorPost;
    }
    else {
      if (selected) {
        intensity = SET_INTENSITY(frameCurrent, frame, frameEnd, 0.25, 0.75);
      }
      else {
        intensity = SET_INTENSITY(frameCurrent, frame, frameEnd, 0.68, 0.92);
      }

      out_finalcolor.rgb = mix(colorBonePose.rgb, blend_base, intensity);
    }
  }
  else {
    if (use_custom_color) {
      out_finalcolor.rgb = colorCurrentFrame.rgb;
    }
    else {
      if (selected) {
        intensity = 0.92f;
      }
      else {
        intensity = 0.75f;
      }
      out_finalcolor.rgb = mix(colorCurrentFrame.rgb, blend_base, intensity);
    }
  }

  out_finalcolor.a = 1.0;
}

void main()
{
  /** Determine Output Primitive ID and relative vertex. */
  /* Index of the quad primitive. We generate one quad for each input line. */
  int quad_id = gl_VertexID / 6;

  /* Determine vertex within the quad (A, B, C)(A, C, D). */
  int quad_vertex_id = gl_VertexID % 6;
  /* Base index of the line primitive:
   *  - IF PrimType == LineList:  base_vertex_id = quad_id*2
   *  - IF PrimType == LineStrip: base_vertex_id = quad_id
   *
   * NOTE: Primitive is LineStrip for this shader. */
  int base_vertex_id = quad_id;

  /* Fetch attributes for self and neighboring vertex. */
  vec3 in_pos0 = vertex_fetch_attribute(base_vertex_id, pos, vec3);
  vec3 in_pos1 = vertex_fetch_attribute(base_vertex_id + 1, pos, vec3);

  vec4 out_pos0 = drw_view.winmat * (drw_view.viewmat * vec4(in_pos0, 1.0));
  vec4 out_pos1 = drw_view.winmat * (drw_view.viewmat * vec4(in_pos1, 1.0));

  /* Final calculations required for Geometry Shader alternative.
   * We need to calculate values for each vertex position to correctly determine the final output
   * position. */
  vec2 ssPos[2];
  vec4 finalColor_geom[2];

  do_vertex_shader(out_pos0, base_vertex_id, ssPos[0], finalColor_geom[0]);
  do_vertex_shader(out_pos1, base_vertex_id + 1, ssPos[1], finalColor_geom[1]);

  /* Geometry shader alternative -- Output is triangle-list consisting of 6 vertices.
   * Each vertex shader invocation is one vertex in the output primitive, so output
   * required ID. */
  vec2 t;
  vec2 edge_dir = compute_dir(ssPos[0], ssPos[1]) * sizeViewportInv;

  bool is_persp = (ProjectionMatrix[3][3] == 0.0);
  float line_size = float(lineThickness) * sizePixel;

  if (quad_vertex_id == 0) {
    view_clipping_distances(out_pos0.xyz);

    interp.color = finalColor_geom[0];
    t = edge_dir * (line_size * (is_persp ? out_pos0.w : 1.0));
    gl_Position = out_pos0 + vec4(t, 0.0, 0.0);
  }
  else if (quad_vertex_id == 1 || quad_vertex_id == 3) {
    view_clipping_distances(out_pos0.xyz);

    interp.color = finalColor_geom[0];
    t = edge_dir * (line_size * (is_persp ? out_pos0.w : 1.0));
    gl_Position = out_pos0 - vec4(t, 0.0, 0.0);
  }
  else if (quad_vertex_id == 2 || quad_vertex_id == 5) {
    view_clipping_distances(out_pos1.xyz);

    interp.color = finalColor_geom[1];
    t = edge_dir * (line_size * (is_persp ? out_pos1.w : 1.0));
    gl_Position = out_pos1 + vec4(t, 0.0, 0.0);
  }
  else if (quad_vertex_id == 4) {
    view_clipping_distances(out_pos1.xyz);

    interp.color = finalColor_geom[1];
    t = edge_dir * (line_size * (is_persp ? out_pos1.w : 1.0));
    gl_Position = out_pos1 - vec4(t, 0.0, 0.0);
  }
}
