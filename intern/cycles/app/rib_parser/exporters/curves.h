#pragma once

#include "bvh/build.h"
#include "scene/scene.h"

#include "app/rib_parser/scene_entities.h"
#include "util/boundbox.h"
#include <string>
#include <unordered_map>

CCL_NAMESPACE_BEGIN

class RIBCyclesCurves {
 public:
  RIBCyclesCurves(Scene *scene,
                  vector<Instance_Scene_Entity> &inst,
                  Instance_Definition_Scene_Entity *inst_def)
      : _scene(scene), _inst_v(inst), _inst_def(inst_def)
  {
  }

  ~RIBCyclesCurves() = default;

  void export_curves();
  BoundBox const& bounds() const {return _bounds;}
  
 protected:
  void initialize(std::string name);
  void initialize_instance(int index);
  void populate(bool &rebuild);
  void populate_widths();
  void populate_primvars();
  void populate_points();
  void populate_topology();

 private:
  Scene *_scene = nullptr;
  vector<Instance_Scene_Entity> const &_inst_v;
  Instance_Definition_Scene_Entity const *_inst_def;
  Hair *_geom = nullptr;
  std::unordered_map<std::string, Hair *> _instanced_geom;
  vector<Object *> _instances;
  ProjectionTransform _geomTransform;
  BoundBox _bounds{BoundBox::empty};
  Shape_Scene_Entity _shape;
};

CCL_NAMESPACE_END
