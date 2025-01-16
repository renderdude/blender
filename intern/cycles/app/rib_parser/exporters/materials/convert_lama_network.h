#pragma once

#include "app/rib_parser/param_dict.h"
#include "app/rib_parser/parsed_parameter.h"
#include <set>
#include <utility>

CCL_NAMESPACE_BEGIN

class LamaNetwork {
 public:
  LamaNetwork(Vector_Dictionary &shader_graph);
  ~LamaNetwork() = default;

  Vector_Dictionary convert();
  bool has_emission()
  {
    return _has_emission_node;
  }

 private:
  bool generate_osl(std::string shader_name);
  void generate_mtlx_definition();
  void remove_external_nodes();
  void find_common_references();
  void find_parameters();
  void remap_parameters();
  void split_nodegraph();
  void map_renderman_to_mtlx();
  void match_renderman_definitions();
  std::string generate_parameters();
  std::string generate_nodegraph();
  std::string rescale_parameters();
  void adjust_connections();

  // RenderMan to MaterialX helper routines
  void flip_int(std::string name, std::string param_name, Parameter_Dictionary &params, int value);
  void color_to_float(std::string name, std::string param_name);
  void to_vector(std::string name, std::string param_name);

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

  void update_remapped_name(std::string node_name, Parsed_Parameter *param, std::string new_name)
  {
    bool found = false;
    if (_remapped_params.find(node_name) != _remapped_params.end()) {
      auto &remap = _remapped_params[node_name];
      if (remap.find(param) != remap.end()) {
        remap[param] = new_name;
        found = true;
      }
    }

    // Check "common" if not found
    if (!found) {
      if (_remapped_params.find("common") != _remapped_params.end()) {
        auto &remap = _remapped_params["common"];
        if (remap.find(param) != remap.end()) {
          remap[param] = new_name;
          found = true;
        }
      }
    }
  }

  bool _a_node_was_split = false;
  bool _has_emission_node = false;
  Vector_Dictionary _shader_graph;
  Vector_Dictionary _lama_shader_graph;
  std::string _mtlx_def;
  Parameter_Dictionary _lama_surface;
  std::map<std::string, Parameter_Dictionary *> _non_lama_nodes;
  std::map<std::string, vector<Parsed_Parameter *>> _common_ext_refs, _constants,
      _external_references, _internal_references, _rescaled_parameters;
  std::map<std::string, Parameter_Dictionary> _handle_to_params;
  std::map<std::string, std::map<Parsed_Parameter *, std::string>> _remapped_params;
  std::map<std::string, std::pair<std::string, bool>> _handle_to_lama;
  std::map<std::string, std::set<std::string>> _prman_lama_params;
};

CCL_NAMESPACE_END
