#pragma once

#include "bvh/build.h"
#include "scene/hair.h"
#include "scene/scene.h"

#include "app/rib_parser/scene_entities.h"
#include "util/boundbox.h"
#include <string>
#include <unordered_map>

CCL_NAMESPACE_BEGIN

class RIBCyclesCurves {
 public:
  RIBCyclesCurves(Scene *scene) : _scene(scene) {}

  ~RIBCyclesCurves() = default;

  void build_instance_definition(Instance_Definition_Scene_Entity const *inst_def);
  void build_instance(Instance_Scene_Entity &inst);

  BoundBox const &bounds() const
  {
    return _bounds;
  }

  array<Node *> get_used_shaders() {
    return _geom->get_used_shaders();
  }

 protected:
  void initialize(std::string name);
  std::string initialize_instance(Instance_Scene_Entity &inst);
  void populate();
  void populate_widths();
  void populate_primvars();
  void populate_points();
  void populate_topology();

 private:
  Scene *_scene = nullptr;
  Hair *_geom = nullptr;
  std::unordered_map<std::string, Hair *> _instanced_geom;
  Object* _instance;
  ProjectionTransform _geomTransform;
  BoundBox _bounds{BoundBox::empty};
  Shape_Scene_Entity _shape;
};

CCL_NAMESPACE_END
