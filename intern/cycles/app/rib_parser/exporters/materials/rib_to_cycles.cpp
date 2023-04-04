
#include "app/rib_parser/exporters/materials/rib_to_cycles.h"

CCL_NAMESPACE_BEGIN

bool create_shader_node(std::string const &nodeType,
                        std::string const &shader,
                        std::string const &path,
                        ShaderGraph *graph,
                        Scene *scene,
                        ShaderNode **node)
{
  std::string shader_name = shader;
  bool check_if_osl = false;
  bool result = true;
  if (const NodeType *node_type = NodeType::find(ustring(nodeType))) {
    *node = static_cast<ShaderNode *>(node_type->create(node_type));
  }
  else {
    check_if_osl = true;
    auto sn = nodeType;
    if (!sn.empty())
      shader_name = nodeType;
  }

  if (check_if_osl) {
    if (string_startswith(shader_name, "Pxr"))
      shader_name = shader_name + ".oso";
    if (string_endswith(shader_name, ".oso")) {
      if (path_is_relative(shader_name))
        shader_name = path_join(path, shader_name);
      *node = OSLShaderManager::osl_node(graph, scene->shader_manager, shader_name, "");
    }
    else {
      fprintf(stderr, "Could not create node '%s'", shader_name.c_str());
      result = false;
    }
  }
  return result;
}

const SocketType *find_socket(std::string input_name, ShaderNode *node)
{
  const SocketType *input = nullptr;
  for (const SocketType &socket : node->type->inputs) {
    if (string_iequals(socket.name.string(), input_name) || socket.ui_name == input_name) {
      input = &socket;
      break;
    }
  }

  return input;
}

bool RIBtoCyclesMapping::create_shader_node(std::string const &shader,
                                            std::string const &path,
                                            ShaderGraph *graph,
                                            Scene *scene)
{
  bool result = false;
  ShaderNode *node = nullptr;

  if (::ccl::create_shader_node(_nodeType[0], shader, path, graph, scene, &node)) {
    _nodes.push_back(node);
    result = true;
  }
  return result;
}

void RIBtoCyclesMapping::update_parameters(Parameter_Dictionary const &parameters,
                                           vector<Parsed_Parameter const *> &connections)
{
  for (const auto param : parameters.get_parameter_vector()) {
    // Check if the parameter is a connection, and defer processing
    // if it is
    if (param->storage == Container_Type::Reference)
      connections.push_back(param);
    else {
      // See if the parameter name is in Pixar terms, and needs to be converted
      const std::string input_name = parameter_name(param->name);

      // Find the input to write the parameter value to
      const SocketType *input = find_socket(input_name, _nodes.back());

      if (!input) {
        VLOG_WARNING << "Could not find parameter '" << param->name.c_str() << "' on node '"
                     << _nodes.back()->name.c_str() << "'\n";
        continue;
      }

      set_node_value(_nodes.back(), *input, param);
    }
  }
}

bool RIBtoMultiNodeCycles::create_shader_node(std::string const &shader,
                                              std::string const &path,
                                              ShaderGraph *graph,
                                              Scene *scene)
{
  std::string shader_name = shader;
  bool result = true;
  ShaderNode *node = nullptr;

  for (auto &node_type : _nodeType)
    if (::ccl::create_shader_node(node_type, shader, path, graph, scene, &node)) {
      _nodes.push_back(node);
      _node_map[node_type] = node;
    }
    else
      result = false;

  // Create the internal connection map
  for (auto pp : _connectionMap) {
    vector<string> key_tokens, map_tokens;
    string_split(key_tokens, pp.first, ":");
    string_split(map_tokens, pp.second.string(), ":");

    ShaderNode *map_node = _node_map[map_tokens[0]];

    // Find the input to connect to on the passed in node
    ShaderInput *input = nullptr;
    for (ShaderInput *in : map_node->inputs) {
      if (string_iequals(in->socket_type.name.string(), map_tokens[1])) {
        input = in;
        break;
      }
    }

    if (!input) {
      fprintf(stderr,
              "Ignoring connection on '%s.%s', input '%s' was not found\n",
              map_node->name.c_str(),
              map_tokens[0].c_str(),
              map_tokens[1].c_str());
      continue;
    }

    ShaderNode *key_node = _node_map[key_tokens[0]];
    ShaderOutput *output = nullptr;
    for (ShaderOutput *out : key_node->outputs) {
      if (string_iequals(out->socket_type.name.string(), key_tokens[1])) {
        output = out;
        break;
      }
    }

    if (!output) {
      fprintf(stderr,
              "Ignoring connection from '%s.%s' to '%s.%s', output '%s' was not found\n",
              key_tokens[0].c_str(),
              key_tokens[1].c_str(),
              map_node->name.c_str(),
              map_tokens[0].c_str(),
              map_tokens[1].c_str());
      continue;
    }

    graph->connect(output, input);
  }

  return result;
}

void RIBtoMultiNodeCycles::update_parameters(Parameter_Dictionary const &parameters,
                                             vector<Parsed_Parameter const *> &connections)
{
  for (const auto param : parameters.get_parameter_vector()) {
    // Check if the parameter is a connection, and defer processing
    // if it is
    if (param->storage == Container_Type::Reference)
      connections.push_back(param);
    else {
      // See if the parameter name is in Pixar terms, and needs to be converted
      const std::string input_name = parameter_name(param->name);

      if (input_name.find(":") != std::string::npos) {
        vector<string> tokens;
        string_split(tokens, input_name, ":");
        ShaderNode *node = _node_map[tokens[0]];

        // Find the input to write the parameter value to
        const SocketType *input = find_socket(tokens[1], node);

        if (!input) {
          VLOG_WARNING << "Could not find parameter '" << tokens[1] << "' on node '"
                       << node->name.c_str() << "'\n";
          continue;
        }

        set_node_value(node, *input, param);
      }
    }
  }
}

bool PxrNormalMaptoCycles::create_shader_node(std::string const &shader,
                                              std::string const &path,
                                              ShaderGraph *graph,
                                              Scene *scene)
{
  std::string shader_name = shader;
  bool result = true;
  ShaderNode *node = nullptr;

  for (auto &node_type : _nodeType)
    if (::ccl::create_shader_node(node_type, shader, path, graph, scene, &node)) {
      _nodes.push_back(node);
      _node_map[node_type] = node;
      if (node->is_a(ImageTextureNode::node_type)) {
        ImageTextureNode *itn = (ImageTextureNode *)node;
        itn->set_colorspace(ustring("Non-Color"));
      }
      else if (node->is_a(NormalMapNode::node_type)) {
        NormalMapNode *itn = (NormalMapNode *)node;
      }
    }
    else
      result = false;

  // Create the internal connection map
  for (auto pp : _connectionMap) {
    vector<string> key_tokens, map_tokens;
    string_split(key_tokens, pp.first, ":");
    string_split(map_tokens, pp.second.string(), ":");

    ShaderNode *map_node = _node_map[map_tokens[0]];

    // Find the input to connect to on the passed in node
    ShaderInput *input = nullptr;
    for (ShaderInput *in : map_node->inputs) {
      if (string_iequals(in->socket_type.name.string(), map_tokens[1])) {
        input = in;
        break;
      }
    }

    if (!input) {
      fprintf(stderr,
              "Ignoring connection on '%s.%s', input '%s' was not found\n",
              map_node->name.c_str(),
              map_tokens[0].c_str(),
              map_tokens[1].c_str());
      continue;
    }

    ShaderNode *key_node = _node_map[key_tokens[0]];
    ShaderOutput *output = nullptr;
    for (ShaderOutput *out : key_node->outputs) {
      if (string_iequals(out->socket_type.name.string(), key_tokens[1])) {
        output = out;
        break;
      }
    }

    if (!output) {
      fprintf(stderr,
              "Ignoring connection from '%s.%s' to '%s.%s', output '%s' was not found\n",
              key_tokens[0].c_str(),
              key_tokens[1].c_str(),
              map_node->name.c_str(),
              map_tokens[0].c_str(),
              map_tokens[1].c_str());
      continue;
    }

    graph->connect(output, input);
  }

  return result;
}

bool RIBtoCyclesTexture::create_shader_node(std::string const &shader,
                                            std::string const &path,
                                            ShaderGraph *graph,
                                            Scene *scene)
{
  bool result = false;
  ShaderNode *node = nullptr;

  if (::ccl::create_shader_node(_nodeType[0], shader, path, graph, scene, &node)) {
    _nodes.push_back(node);
    result = true;
    if (node->is_a(ImageTextureNode::node_type)) {
      ImageTextureNode *itn = (ImageTextureNode *)node;
      itn->set_colorspace(ustring("sRGB"));
    }
  }
  return result;
}

void PxrSurfacetoPrincipled::update_parameters(Parameter_Dictionary const &parameters,
                                               vector<Parsed_Parameter const *> &connections)
{
  // Exactly the same as the base except we create a map of the parameters for
  // later retrieval
  for (const auto param : parameters.get_parameter_vector()) {
    // Check if the parameter is a connection, and defer processing
    // if it is
    _parameters[param->name] = param;
    if (param->storage == Container_Type::Reference)
      connections.push_back(param);
    else {
      // See if the parameter name is in Pixar terms, and needs to be converted
      const std::string input_name = parameter_name(param->name);

      // Find the input to write the parameter value to
      const SocketType *input = find_socket(input_name, _nodes.back());

      if (!input) {
        VLOG_WARNING << "Could not find parameter '" << param->name.c_str() << "' on node '"
                     << _nodes.back()->name.c_str() << "'\n";
        continue;
      }

      set_node_value(_nodes.back(), *input, param);
    }
  }

  // Now handle the funny one-offs that require remapping
  Parsed_Parameter updated_param;
  Parsed_Parameter *param;
  const SocketType *input;

  // Transmission
  param = _parameters["refractionGain"];
  if (param->floats[0] > 0) {  // Some Trasmission is set
    updated_param.floats.clear();
    updated_param.type = "float";
    updated_param.add_float(param->floats[0]);
    input = find_socket("transmission", _nodes.back());
    set_node_value(_nodes.back(), *input, &updated_param);

    param = _parameters["glassRoughness"];
    updated_param.floats.clear();
    updated_param.type = "float";
    updated_param.add_float(param->floats[0]);
    input = find_socket("roughness", _nodes.back());
    set_node_value(_nodes.back(), *input, &updated_param);

    param = _parameters["glassIor"];
    updated_param.floats.clear();
    updated_param.type = "float";
    updated_param.add_float(param->floats[0]);
    input = find_socket("ior", _nodes.back());
    set_node_value(_nodes.back(), *input, &updated_param);

    updated_param.floats.clear();
    updated_param.type = "float";
    updated_param.add_float(0.);
    input = find_socket("specular", _nodes.back());
    set_node_value(_nodes.back(), *input, &updated_param);

    updated_param.floats.clear();
    updated_param.type = "color";
    updated_param.add_float(1);
    updated_param.add_float(1);
    updated_param.add_float(1);
    input = find_socket("base_color", _nodes.back());
    set_node_value(_nodes.back(), *input, &updated_param);
  }
  else {
    // diffuse gain
    float gain = _parameters["diffuseGain"]->floats[0];
    param = _parameters["diffuseColor"];
    if (param->storage != Container_Type::Reference) {
      updated_param.floats.clear();
      updated_param.type = "color";
      updated_param.add_float(gain * param->floats[0]);
      updated_param.add_float(gain * param->floats[1]);
      updated_param.add_float(gain * param->floats[2]);
      input = find_socket("base_color", _nodes.back());
      set_node_value(_nodes.back(), *input, &updated_param);
    }

    // Specular
    param = _parameters["specularFresnelMode"];
    if (param->ints[0] == 0) {  // Artistic Mode
      param = _parameters["specularFaceColor"];
      if (param->storage != Container_Type::Reference) {
        float lum = 0.2126 * param->floats[0] + 0.7152 * param->floats[1] +
                    0.0722 * param->floats[2];
        updated_param.floats.clear();
        updated_param.type = "float";
        updated_param.add_float(lum);
        input = find_socket("specular", _nodes.back());
        set_node_value(_nodes.back(), *input, &updated_param);
      }
    }
  }
}

CCL_NAMESPACE_END
