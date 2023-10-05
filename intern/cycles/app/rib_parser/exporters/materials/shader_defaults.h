#ifndef SHADER_DEFAULTS_H
#define SHADER_DEFAULTS_H

#include "graph/node_type.h"
#include "scene/shader_graph.h"
#include "scene/shader_nodes.h"
CCL_NAMESPACE_BEGIN

void fill(PrincipledBsdfNode *node)
{
  node->set_base_color(make_float3(0.18f, 0.18f, 0.18f));
  node->set_metallic(0.0f);
  node->set_subsurface_weight(0.0f);
  node->set_subsurface_scale(0.1f);
  node->set_subsurface_radius(make_float3(0.1f, 0.1f, 0.1f));
  node->set_subsurface_ior(1.4f);
  node->set_subsurface_anisotropy(0.0f);
  node->set_specular_ior_level(1.0f);
  node->set_roughness(0.25f);
  node->set_specular_tint(zero_float3());
  node->set_anisotropic(0.0f);
  node->set_sheen_weight(0.0f);
  node->set_sheen_roughness(0.5f);
  node->set_sheen_tint(one_float3());
  node->set_coat_weight(0.0f);
  node->set_coat_roughness(0.03f);
  node->set_coat_ior(1.5f);
  node->set_coat_tint(one_float3());
  node->set_ior(1.5f);
  node->set_transmission_weight(0.0f);
  node->set_anisotropic_rotation(0.0f);
  node->set_emission_color(one_float3());
  node->set_emission_strength(0.0f);
  node->set_alpha(1.0f);
}

void fill(DiffuseBsdfNode *node)
{
  node->set_color(make_float3(0.5, 0.5, 0.5));
}

ShaderNode *create_shader(const NodeType *node_type)
{
  ShaderNode *node = static_cast<ShaderNode *>(node_type->create(node_type));

  if (node->is_a(PrincipledBsdfNode::node_type)) {
    fill((PrincipledBsdfNode *)node);
  }
  else if (node->is_a(DiffuseBsdfNode::node_type)) {
    fill((DiffuseBsdfNode *)node);
  }

  return node;
}

CCL_NAMESPACE_END

#endif  // SHADER_DEFAULTS_H
