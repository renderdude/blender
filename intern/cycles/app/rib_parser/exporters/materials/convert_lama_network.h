#ifndef CONVERT_LAMA_NETWORK_H
#define CONVERT_LAMA_NETWORK_H

#include "app/rib_parser/param_dict.h"
#include "app/rib_parser/parsed_parameter.h"
#include <utility>

CCL_NAMESPACE_BEGIN

using Vector_Dictionary = std::pair<std::string, std::vector<Parameter_Dictionary>>;

class LamaNetwork {
 public:
  LamaNetwork(Vector_Dictionary &shader_graph) : _shader_graph(shader_graph)
  {
  }
  ~LamaNetwork()
  {
  }

  Vector_Dictionary convert();

 private:
  bool generate_osl(std::string shader_name);
  void generate_mtlx_definition();
  void remove_external_nodes();
  void find_common_references();
  void find_parameters();
  void remap_parameters();
	void split_nodegraph();
	void match_renderman_definitions();
  std::string generate_parameters();
  std::string generate_nodegraph();
  void adjust_connections();

  std::string remapped_name(std::string node_name, Parsed_Parameter *param, std::string def_name)
  {
    bool found = false;
    std::string iface_name = def_name;
    if (_remapped_params.find(node_name) != _remapped_params.end()) {
      auto remap = _remapped_params[node_name];
      if (remap.find(param) != remap.end()) {
        iface_name = remap[param];
        found = true;
      }
    }

    // Check "common" if not found
    if (!found) {
      if (_remapped_params.find("common") != _remapped_params.end()) {
        auto remap = _remapped_params["common"];
        if (remap.find(param) != remap.end()) {
          iface_name = remap[param];
          found = true;
        }
      }
    }
    return iface_name;
  }

  bool _a_node_was_split = false;
  Vector_Dictionary &_shader_graph;
  Vector_Dictionary _lama_shader_graph;
  std::string _mtlx_def;
  Parameter_Dictionary _lama_surface;
  std::map<std::string, Parameter_Dictionary *> _non_lama_nodes;
  std::map<std::string, vector<Parsed_Parameter *>> _common_ext_refs;
  std::map<std::string, vector<Parsed_Parameter *>> _constants, _external_references,
      _internal_references;
  std::map<std::string, Parameter_Dictionary> _handle_to_params;
  std::map<std::string, std::map<Parsed_Parameter *, std::string>> _remapped_params;
  std::map<std::string, std::pair<std::string, bool>> _handle_to_lama;
};

CCL_NAMESPACE_END

#endif  // CONVERT_LAMA_NETWORK_H
