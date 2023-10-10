#include "app/rib_parser/exporters/materials/rib_to_cycles.h"
#include "app/rib_parser/error.h"
#include "app/rib_parser/exporters/materials/shader_defaults.h"
#include "app/rib_parser/parsed_parameter.h"
#include "util/color.h"

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
    *node = create_shader(node_type);
  }
  else {
    check_if_osl = true;
    auto sn = nodeType;
    if (!sn.empty()) {
      shader_name = nodeType;
    }
  }

  if (check_if_osl) {
    if (!string_endswith(shader_name, ".oso")) {
      shader_name = shader_name + ".oso";
    }
    if (path_is_relative(shader_name)) {
      shader_name = path_join(path, shader_name);
    }
    *node = OSLShaderManager::osl_node(graph, scene->shader_manager, shader_name, "");
    if (!node) {
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
  for (auto *const param : parameters.get_parameter_vector()) {
    // Check if the parameter is a connection, and defer processing
    // if it is
    if (param->storage == Container_Type::Reference) {
      connections.push_back(param);
    }
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

  for (auto &node_type : _nodeType) {
    if (::ccl::create_shader_node(node_type, shader, path, graph, scene, &node)) {
      _nodes.push_back(node);
      _node_map[node_type] = node;
    }
    else {
      result = false;
    }
  }

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
  for (auto *const param : parameters.get_parameter_vector()) {
    // Check if the parameter is a connection, and defer processing
    // if it is
    if (param->storage == Container_Type::Reference) {
      connections.push_back(param);
    }
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

void PxrNormalMaptoCycles::update_parameters(Parameter_Dictionary const &parameters,
                                             vector<Parsed_Parameter const *> &connections)
{
  // Exactly the same as the base except we create a map of the parameters for
  // later retrieval
  for (auto *const param : parameters.get_parameter_vector()) {
    // Check if the parameter is a connection, and defer processing
    // if it is
    _parameters[param->name] = param;
    if (param->storage == Container_Type::Reference) {
      connections.push_back(param);
    }
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

  updated_param.payload = vector<float>();
  param = _parameters["inputRGB"];
  if (param) {
    if (param->storage == Container_Type::Reference) {
    }
  }
}

void PxrRamptoCycles::update_parameters(Parameter_Dictionary const &parameters,
                                        vector<Parsed_Parameter const *> &connections)
{
  // Exactly the same as the base except we create a map of the parameters for
  // later retrieval
  for (auto *const param : parameters.get_parameter_vector()) {
    // Check if the parameter is a connection, and defer processing
    // if it is
    _parameters[param->name] = param;
    if (param->storage == Container_Type::Reference) {
      connections.push_back(param);
    }
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

  updated_param.payload = vector<float>();
  param = _parameters["rampType"];
  int rampType = 0;
  if (param) {
    rampType = param->ints()[0];
  }

  if (rampType != 0) {
    std::cerr << "Only PxrRamp with \"S-Ramp\" mapping is currently handled." << std::endl;
  }

  param = _parameters["colorRamp"];
  if (param) {
    int num_knots = param->ints()[0];
    auto knots = _parameters["colorRamp_Knots"]->floats();
    auto colors = _parameters["colorRamp_Colors"]->floats();
    Parsed_Parameter *interp = _parameters["colorRamp_Interpolation"];
    if (!interp || interp->strings()[0] != "linear") {
      std::cerr << "Only linear interpolation for now with \"colorRamp_Interpolation\""
                << std::endl;
    }

    Parsed_Parameter out_color(Parameter_Type::Color, "color", File_Loc());
    Parsed_Parameter out_alpha(Parameter_Type::Real, "alpha", File_Loc());
    float dt = 1.0 / 256.0;
    float eval_pt = 0.0f;
    int knot_index = 0;
    int index = 0;
    while (index < 256) {
      float upper_bound = knots[knot_index];
      float3 left, right;
      if (knot_index > 0 && (eval_pt >= knots[knot_index - 1] && eval_pt < knots[knot_index])) {
        left = make_float3(colors[3 * (knot_index - 1)],
                           colors[3 * (knot_index - 1) + 1],
                           colors[3 * (knot_index - 1) + 2]);
        right = make_float3(
            colors[3 * knot_index], colors[3 * knot_index + 1], colors[3 * knot_index + 2]);
      }
      else if (eval_pt < knots[knot_index]) {
        right = make_float3(
            colors[3 * knot_index], colors[3 * knot_index + 1], colors[3 * knot_index + 2]);
        left = right;
      }
      else if (eval_pt > knots[3]) {
        left = make_float3(
            colors[3 * knot_index], colors[3 * knot_index + 1], colors[3 * knot_index + 2]);
        right = left;
        upper_bound = 1.0f;
      }
      else {
        knot_index++;
        continue;  // Duplicate Knot point and we've strided over it
      }
      while (index < 256 && eval_pt < upper_bound) {
        float3 result = lerp(left, right, eval_pt);
        out_color.add_float(result[0]);
        out_color.add_float(result[1]);
        out_color.add_float(result[2]);
        out_alpha.add_float(1.0);
        eval_pt += dt;
        index++;
      }
      knot_index++;
    }
    assert(out_color.floats().size() == 768);
    const SocketType *input = find_socket("ramp", _nodes.back());
    set_node_value(_nodes.back(), *input, &out_color);
    input = find_socket("ramp_alpha", _nodes.back());
    set_node_value(_nodes.back(), *input, &out_alpha);
  }
}

bool PxrImageNormalMaptoCycles::create_shader_node(std::string const &shader,
                                                   std::string const &path,
                                                   ShaderGraph *graph,
                                                   Scene *scene)
{
  std::string shader_name = shader;
  bool result = true;
  ShaderNode *node = nullptr;

  for (auto &node_type : _nodeType) {
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
    else {
      result = false;
    }
  }

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
  for (auto *const param : parameters.get_parameter_vector()) {
    // Check if the parameter is a connection, and defer processing
    // if it is
    _parameters[param->name] = param;
    if (param->storage == Container_Type::Reference) {
      connections.push_back(param);
    }
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
  if (param && param->floats()[0] > 0) {  // Some Transmission is set
    // The payload is in an undefinded state, so force it to floats
    updated_param.payload = vector<float>();
    updated_param.type = Parameter_Type::Real;
    updated_param.add_float(param->floats()[0]);
    input = find_socket("transmission", _nodes.back());
    set_node_value(_nodes.back(), *input, &updated_param);

    param = _parameters["glassRoughness"];
    updated_param.floats().clear();
    updated_param.type = Parameter_Type::Real;
    updated_param.add_float(param->floats()[0]);
    input = find_socket("roughness", _nodes.back());
    set_node_value(_nodes.back(), *input, &updated_param);

    param = _parameters["glassIor"];
    updated_param.floats().clear();
    updated_param.type = Parameter_Type::Real;
    updated_param.add_float(param->floats()[0]);
    input = find_socket("ior", _nodes.back());
    set_node_value(_nodes.back(), *input, &updated_param);

    updated_param.floats().clear();
    updated_param.type = Parameter_Type::Real;
    updated_param.add_float(0.);
    input = find_socket("specular", _nodes.back());
    set_node_value(_nodes.back(), *input, &updated_param);

    updated_param.floats().clear();
    updated_param.type = Parameter_Type::Color;
    updated_param.add_float(1);
    updated_param.add_float(1);
    updated_param.add_float(1);
    input = find_socket("base_color", _nodes.back());
    set_node_value(_nodes.back(), *input, &updated_param);
  }
  else {
    // The payload is in an undefinded state, so force it to floats
    updated_param.payload = vector<float>();
    // diffuse gain
    float gain = 1.0;
    if ((param = _parameters["diffuseGain"])) {
      gain = param->floats()[0];
    }
    if ((param = _parameters["diffuseColor"])) {
      if (param->storage != Container_Type::Reference) {
        updated_param.floats().clear();
        updated_param.type = Parameter_Type::Color;
        updated_param.add_float(gain * param->floats()[0]);
        updated_param.add_float(gain * param->floats()[1]);
        updated_param.add_float(gain * param->floats()[2]);
        input = find_socket("base_color", _nodes.back());
        set_node_value(_nodes.back(), *input, &updated_param);
      }
    }

    // Specular
    if ((param = _parameters["specularFresnelMode"])) {
      if (param->ints()[0] == 0) {  // Artistic Mode
        if ((param = _parameters["specularFaceColor"])) {
          if (param->storage != Container_Type::Reference) {
            float lum = 0.2126 * param->floats()[0] + 0.7152 * param->floats()[1] +
                        0.0722 * param->floats()[2];
            updated_param.floats().clear();
            updated_param.type = Parameter_Type::Real;
            updated_param.add_float(lum);
            input = find_socket("specular_ior_level", _nodes.back());
            set_node_value(_nodes.back(), *input, &updated_param);
          }
        }
      }
    }
  }
}

void PxrDisneytoPrincipled::update_parameters(Parameter_Dictionary const &parameters,
                                              vector<Parsed_Parameter const *> &connections)
{
  // Exactly the same as the base except we create a map of the parameters for
  // later retrieval
  for (auto *const param : parameters.get_parameter_vector()) {
    // Check if the parameter is a connection, and defer processing
    // if it is
    _parameters[param->name] = param;
    if (param->storage == Container_Type::Reference) {
      connections.push_back(param);
    }
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

  updated_param.payload = vector<float>();
  if ((param = _parameters["emitColor"])) {
    if (param->storage != Container_Type::Reference) {
      updated_param.type = Parameter_Type::Real;
      float3 hsv = rgb_to_hsv(
          make_float3(param->floats()[0], param->floats()[1], param->floats()[2]));
      updated_param.add_float(hsv[2]);
      input = find_socket("emission_strength", _nodes.back());
      set_node_value(_nodes.back(), *input, &updated_param);
    }
  }
}

void PxrDisneyBsdftoPrincipled::update_parameters(Parameter_Dictionary const &parameters,
                                                  vector<Parsed_Parameter const *> &connections)
{
  // Exactly the same as the base except we create a map of the parameters for
  // later retrieval
  for (auto *const param : parameters.get_parameter_vector()) {
    // Check if the parameter is a connection, and defer processing
    // if it is
    _parameters[param->name] = param;
    if (param->storage == Container_Type::Reference) {
      connections.push_back(param);
    }
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

  updated_param.payload = vector<float>();
  if ((param = _parameters["emitColor"])) {
    if (param->storage != Container_Type::Reference) {
      updated_param.type = Parameter_Type::Real;
      float3 hsv = rgb_to_hsv(
          make_float3(param->floats()[0], param->floats()[1], param->floats()[2]));
      updated_param.add_float(hsv[2]);
      input = find_socket("emission_strength", _nodes.back());
      set_node_value(_nodes.back(), *input, &updated_param);
    }
  }
  if ((param = _parameters["diffTrans"])) {
    if (param->storage != Container_Type::Reference) {
      updated_param.floats().clear();
      updated_param.type = Parameter_Type::Real;
      // RenderMan diffTrans is in the range [0,2] so rescale back to [0,1]
      updated_param.add_float(param->floats()[0] / 2.0f);
      input = find_socket("transmission", _nodes.back());
      set_node_value(_nodes.back(), *input, &updated_param);
    }
  }
}

CCL_NAMESPACE_END
