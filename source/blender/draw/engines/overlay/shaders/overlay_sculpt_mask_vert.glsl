/* SPDX-FileCopyrightText: 2016-2022 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "common_view_clipping_lib.glsl"
#include "common_view_lib.glsl"

void main()
{
  vec3 world_pos = point_object_to_world(pos);
  gl_Position = point_world_to_ndc(world_pos);

  faceset_color = mix(vec3(1.0), fset, faceSetsOpacity);
  mask_color = 1.0 - (msk * maskOpacity);

  view_clipping_distances(world_pos);
}
