#pragma once

#include "app/rib_parser/parsed_parameter.h"
#include "bvh/build.h"
#include "scene/mesh.h"
#include "scene/scene.h"

#include "app/rib_parser/scene_entities.h"
#include "util/boundbox.h"
#include <string>

CCL_NAMESPACE_BEGIN

class RIBCyclesMesh {
 public:
  RIBCyclesMesh(Scene *scene)
      : _scene(scene)
  {
  }

  ~RIBCyclesMesh() = default;

  void build_instance_definition(Instance_Definition_Scene_Entity const *inst_def);
  void build_instance(Instance_Scene_Entity &inst);

  void export_geometry();

  BoundBox const& bounds() const {return _bounds;}

  array<Node *> get_used_shaders() {
    return _geom->get_used_shaders();
  }

 protected:
  void initialize(std::string name);
  std::string initialize_instance(Instance_Scene_Entity &inst);
  void populate();
  void separate_face_varying_normals();
  void populate_normals();
  void populate_primvars();
  void populate_points();
  void populate_topology();
  void populate_shader_graph(bool initializing = false);
  void create_uv_map(Parsed_Parameter* param);

  Parsed_Parameter* compute_triangulated_uniform_primvar(const Parsed_Parameter* param);
  Parsed_Parameter* compute_triangulated_face_varying_primvar(const Parsed_Parameter* param);
  Shape_Scene_Entity reduce_geometry_by_faceset(Shape_Scene_Entity const& shape, vector<int> const& faceset);
  void normalize_emission(Object* instance);

 private:
  Scene *_scene = nullptr;
  Mesh *_geom = nullptr;
  Object * _instance;
  Parsed_Parameter* _uv_param = nullptr;
  ProjectionTransform _geomTransform;
  vector<int3> triangles;
  BoundBox _bounds{BoundBox::empty};
  Shape_Scene_Entity _shape;
  bool needs_emission_normalization = false;

  void compute_triangle_indices(const vector<int>& vertices,
                                const vector<int>& nvertices,
                                vector<int3> &indices);
};

CCL_NAMESPACE_END
