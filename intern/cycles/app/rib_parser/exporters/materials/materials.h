
#ifndef __EXPORTERS_MATERIALS_H__
#define __EXPORTERS_MATERIALS_H__

#include "app/rib_parser/exporters/materials/convert_lama_network.h"
#include "app/rib_parser/param_dict.h"
#include "app/rib_parser/parsed_parameter.h"
#include "scene/scene.h"
#include "scene/shader.h"
#include "scene/shader_graph.h"

#include "app/rib_parser/scene_entities.h"
#include "util/vector.h"
#include <vector>

CCL_NAMESPACE_BEGIN

class RIBCyclesMaterials {
 public:
  RIBCyclesMaterials(Scene *scene, Mapped_Vector_Dictionary osl_shader_group)
      : _scene(scene), _osl_shader_group(osl_shader_group)
  {
  }

  ~RIBCyclesMaterials()
  {
  }

  void export_materials();

 protected:
  void initialize();
  void update_connections(class RIBtoCyclesMapping *mapping,
                          ShaderGraph *shader_graph,
                          vector<Parsed_Parameter const *> &pv);
  void populate_shader_graph(Vector_Dictionary shader_graph);

  void add_default_renderman_inputs(Shader *shader);

 private:
  Scene *_scene = nullptr;
  Shader *_shader = nullptr;
  Mapped_Vector_Dictionary _osl_shader_group;
  std::unordered_map<std::string, class RIBtoCyclesMapping *> _nodes;
};

CCL_NAMESPACE_END
#endif  //__EXPORTERS_MATERIALS_H__
