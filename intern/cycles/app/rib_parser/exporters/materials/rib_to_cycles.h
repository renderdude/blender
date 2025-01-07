#pragma once

#include "app/rib_parser/param_dict.h"
#include "app/rib_parser/parsed_parameter.h"
#include "scene/osl.h"
#include "scene/scene.h"
#include "scene/shader_graph.h"
#include "scene/shader_nodes.h"
#include "util/vector.h"
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

CCL_NAMESPACE_BEGIN

class RIBtoCyclesMapping {
 public:
  using ParamMap = std::unordered_map<std::string, ustring>;
  RIBtoCyclesMapping(
      std::vector<std::string> nodeType,
      ParamMap paramMap,
      std::function<void(std::vector<ShaderNode *>)> remapFunc =
          [](std::vector<ShaderNode *>) {})
      : _nodeType(nodeType), _paramMap(std::move(paramMap)), _remapFunc(remapFunc)
  {
  }

  void set(ShaderGraph *graph, Scene *scene, std::unordered_map<std::string, RIBtoCyclesMapping *> *nodes)
  {
    _graph = graph;
    _scene = scene;
    _processed_nodes = nodes;
  }

  ustring nodeType(int index) const
  {
    return ustring(_nodeType[index]);
  }

  virtual std::string parameter_name(const std::string &name) const
  {
    // Simple mapping case
    const auto it = _paramMap.find(name);
    return it != _paramMap.end() ? it->second.string() : name;
  }

  virtual void update_parameters(Parameter_Dictionary const &parameters,
                                 vector<Parsed_Parameter const *> &connections);

  virtual bool create_shader_node(std::string const &shader, std::string const &path);

  virtual ShaderNode *node([[maybe_unused]] std::string name)
  {
    return _nodes.back();
  }

 protected:
  std::vector<ShaderNode *> _nodes;
  std::vector<std::string> _nodeType;
  ParamMap _paramMap;
  ShaderGraph *_graph;
  Scene *_scene;
  std::unordered_map<std::string, RIBtoCyclesMapping *> *_processed_nodes;
  std::function<void(std::vector<ShaderNode *>)> _remapFunc;
};

class RIBtoMultiNodeCycles : public RIBtoCyclesMapping {
 public:
  RIBtoMultiNodeCycles(
      std::vector<std::string> nodeType,
      ParamMap paramMap,
      ParamMap connectionMap,
      std::function<void(std::vector<ShaderNode *>)> remapFunc =
          [](std::vector<ShaderNode *>) {})
      : RIBtoCyclesMapping(nodeType, paramMap, remapFunc), _connectionMap(connectionMap)
  {
  }

  bool create_shader_node(std::string const &shader, std::string const &path) override;

  void update_parameters(Parameter_Dictionary const &parameters,
                         vector<Parsed_Parameter const *> &connections) override;

  ShaderNode *node(std::string name) override
  {
    return _node_map[name];
  }

 protected:
  ParamMap _connectionMap;
  std::map<std::string, ShaderNode *> _node_map;
};

// Specializations
class PxrNormalMaptoCycles : public RIBtoCyclesMapping {
 public:
  using RIBtoCyclesMapping::RIBtoCyclesMapping;

  void update_parameters(Parameter_Dictionary const &parameters,
                         vector<Parsed_Parameter const *> &connections) override;

 private:
  std::unordered_map<std::string, Parsed_Parameter *> _parameters;
};

class PxrRamptoCycles : public RIBtoCyclesMapping {
 public:
  using RIBtoCyclesMapping::RIBtoCyclesMapping;

  void update_parameters(Parameter_Dictionary const &parameters,
                         vector<Parsed_Parameter const *> &connections) override;

 private:
  std::unordered_map<std::string, Parsed_Parameter *> _parameters;
};

class PxrImageNormalMaptoCycles : public RIBtoMultiNodeCycles {
 public:
  using RIBtoMultiNodeCycles::RIBtoMultiNodeCycles;

  bool create_shader_node(std::string const &shader, std::string const &path) override;
};

class RIBtoCyclesTexture : public RIBtoCyclesMapping {
 public:
  using RIBtoCyclesMapping::RIBtoCyclesMapping;

  bool create_shader_node(std::string const &shader, std::string const &path) override;
};

class PxrSurfacetoPrincipled : public RIBtoCyclesMapping {
 public:
  using RIBtoCyclesMapping::RIBtoCyclesMapping;

  void update_parameters(Parameter_Dictionary const &parameters,
                         vector<Parsed_Parameter const *> &connections) override;

 private:
  std::unordered_map<std::string, Parsed_Parameter *> _parameters;
};

class PxrDisneytoPrincipled : public RIBtoCyclesMapping {
 public:
  using RIBtoCyclesMapping::RIBtoCyclesMapping;

  void update_parameters(Parameter_Dictionary const &parameters,
                         vector<Parsed_Parameter const *> &connections) override;

 private:
  std::unordered_map<std::string, Parsed_Parameter *> _parameters;
};

class PxrDisneyBsdftoPrincipled : public RIBtoCyclesMapping {
 public:
  using RIBtoCyclesMapping::RIBtoCyclesMapping;

  void update_parameters(Parameter_Dictionary const &parameters,
                         vector<Parsed_Parameter const *> &connections) override;

 private:
  std::unordered_map<std::string, Parsed_Parameter *> _parameters;
};

class PxrMarschnerHairtoPrincipled : public RIBtoCyclesMapping {
 public:
  using RIBtoCyclesMapping::RIBtoCyclesMapping;

  void update_parameters(Parameter_Dictionary const &parameters,
                         vector<Parsed_Parameter const *> &connections) override;

 private:
  std::unordered_map<std::string, Parsed_Parameter *> _parameters;
};

CCL_NAMESPACE_END
