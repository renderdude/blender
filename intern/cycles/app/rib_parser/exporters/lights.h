#pragma once

#include "scene/light.h"
#include "scene/scene.h"

#include "app/rib_parser/scene_entities.h"

CCL_NAMESPACE_BEGIN

class RIBCyclesLight {
 public:
  RIBCyclesLight(Scene *scene) : _scene(scene) {}

  ~RIBCyclesLight() = default;

  void build(Instance_Definition_Scene_Entity const *inst_def, Instance_Scene_Entity &inst);

  BoundBox const &bounds() const
  {
    return _bounds;
  }

  array<Node *> get_used_shaders()
  {
    return _light->get_used_shaders();
  }

 protected:
  void initialize(Light_Scene_Entity &light_inst, Instance_Definition_Scene_Entity const *inst_def);
  void populate_shader_graph(Light_Scene_Entity &light_inst, bool initializing = false);
  void create_uv_map(Parsed_Parameter *param);

  Parsed_Parameter *compute_triangulated_uniform_primvar(const Parsed_Parameter *param);
  Parsed_Parameter *compute_triangulated_face_varying_primvar(const Parsed_Parameter *param);
  Shape_Scene_Entity reduce_geometry_by_faceset(Shape_Scene_Entity const &shape,
                                                vector<int> const &faceset);
  void normalize_emission(Object *instance);

 private:
  Scene *_scene = nullptr;
  Light *_light = nullptr;
  Object *_instance;
  Parsed_Parameter *_uv_param = nullptr;
  ProjectionTransform _geomTransform;
  vector<int3> triangles;
  BoundBox _bounds{BoundBox::empty};
  Shape_Scene_Entity _shape;
  bool needs_emission_normalization = false;
};

CCL_NAMESPACE_END
