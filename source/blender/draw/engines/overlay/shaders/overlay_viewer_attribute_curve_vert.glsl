/* SPDX-FileCopyrightText: 2017-2022 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

#include "common_view_clipping_lib.glsl"
#include "common_view_lib.glsl"

void main()
{
  vec3 world_pos = point_object_to_world(pos);
  gl_Position = point_world_to_ndc(world_pos);
  finalColor = attribute_value;
}
