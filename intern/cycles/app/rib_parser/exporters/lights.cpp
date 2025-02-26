#include "scene/light.h"
#include "scene/object.h"
#include "scene/shader_graph.h"
#include "scene/shader_nodes.h"

#include "app/rib_parser/exporters/lights.h"
#include "util/transform.h"
#include <OpenImageIO/fmath.h>

CCL_NAMESPACE_BEGIN

Transform convert_transform(const ProjectionTransform &matrix)
{
  return make_transform(-matrix.x[0],
                        -matrix.x[1],
                        -matrix.x[2],
                        matrix.x[3],
                        -matrix.y[0],
                        -matrix.y[1],
                        -matrix.y[2],
                        matrix.y[3],
                        -matrix.z[0],
                        -matrix.z[1],
                        -matrix.z[2],
                        matrix.z[3]);
}

void RIBCyclesLight::build(Instance_Definition_Scene_Entity const *inst_def,
                           Instance_Scene_Entity &inst)
{
  for (auto light_inst : inst_def->lights) {
    const float metersPerUnit = 1.;
    ProjectionTransform xform_obj = *inst.render_from_instance * *light_inst.render_from_light;
    Transform xform = transform_scale(make_float3(metersPerUnit)) * convert_transform(xform_obj);
    vector<Transform> motion = {xform};
    vector<DecomposedTransform> decomp(motion.size());
    transform_motion_decompose(decomp.data(), motion.data(), motion.size());

    initialize(light_inst, inst_def);

    _instance->set_tfm(xform);

    float3 strength = make_float3(1.0f, 1.0f, 1.0f);

    auto color = light_inst.parameters.get_one_color("lightColor", make_float3(1.0f, 1.0f, 1.0f));
    strength = make_float3(color[0], color[1], color[2]);

    float exposure = light_inst.parameters.get_one_float("exposure", 0.0);
    strength *= exp2(exposure);

    float intensity = light_inst.parameters.get_one_float("intensity", 1.0);
    strength *= intensity;

    // Cycles lights are normalized by default, so need to scale intensity if RMan light is not
    bool normalize = light_inst.parameters.get_one_int("areaNormalize", 0) == 1;
    _light->set_normalize(normalize);

    uint visibility = PATH_RAY_ALL_VISIBILITY;
    int rib_vis = inst.parameters.at("visibility").get_one_int("camera", 1);
    if (rib_vis == 0) {
      visibility &= ~PATH_RAY_CAMERA;
    }
    _instance->set_visibility(visibility);

    // Default to shadow casting until we have an example
    _light->set_cast_shadow(true);

    if (light_inst.light_type == "PxrDistantLight") {
      _light->set_angle(OIIO::radians(light_inst.parameters.get_one_float("angle", 0.526f)));
    }
    else if (light_inst.light_type == "PxrDiskLight") {
      const float size = light_inst.parameters.get_one_float("size", 1.f) * 2.0f;
      _light->set_sizeu(size);
      _light->set_sizev(size);
    }
    else if (light_inst.light_type == "PxrRectLight") {
// Size is set in the tfm
#if 1
      _light->set_sizeu(1.f);
      _light->set_sizev(1.f);
#else
      _light->set_sizeu(1.f * fabsf(decomp[0].z.w));
      _light->set_sizev(1.f * fabsf(decomp[0].z.w));
#endif
    }
    else if (light_inst.light_type == "PxrSphereLight") {
#if 1
      _light->set_size(0.5f);
#else
      _light->set_size(0.5f * fabsf(decomp[0].z.w));
#endif

      bool shaping = false;
      const float shapingConeAngle = light_inst.parameters.get_one_float("shapingConeAngle", -1.f);
      if (shapingConeAngle > 0) {
        _light->set_spot_angle(OIIO::radians(shapingConeAngle * 2.0f));
        shaping = true;
      }

      const float shapingConeSoftness = light_inst.parameters.get_one_float("shapingConeSoftness",
                                                                            -1.f);
      if (shapingConeSoftness > 0) {
        _light->set_spot_smooth(shapingConeSoftness);
        shaping = true;
      }

      _light->set_light_type(shaping ? LIGHT_SPOT : LIGHT_POINT);
    }

    _light->set_strength(strength);
    _light->set_is_enabled(true);

    populate_shader_graph(light_inst, false);

    if (_light->is_modified()) {
      _light->tag_update(_scene);
    }
  }
}

void RIBCyclesLight::initialize(Light_Scene_Entity &light_inst,
                                Instance_Definition_Scene_Entity const *inst_def)
{
  _scene->mutex.lock();
  _light = _scene->create_node<Light>();
  _scene->mutex.unlock();
  _light->name = inst_def->name;

  _scene->mutex.lock();
  _instance = _scene->create_node<Object>();
  _scene->mutex.unlock();
  _instance->set_geometry(_light);

  _instance->set_random_id(hash_uint2(hash_string(_light->name.c_str()), 0));

  if (light_inst.light_type == "PxrDomeLight") {
    _light->set_light_type(LIGHT_BACKGROUND);
  }
  else if (light_inst.light_type == "PxrDistantLight") {
    _light->set_light_type(LIGHT_DISTANT);
  }
  else if (light_inst.light_type == "PxrDiskLight") {
    _light->set_light_type(LIGHT_AREA);
    _light->set_ellipse(true);
    _light->set_size(1.0f);
  }
  else if (light_inst.light_type == "PxrRectLight") {
    _light->set_light_type(LIGHT_AREA);
    _light->set_ellipse(false);
    _light->set_size(1.0f);
  }
  else if (light_inst.light_type == "PxrSphereLight") {
    _light->set_light_type(LIGHT_POINT);
    _light->set_size(1.0f);
  }

  _light->set_use_mis(true);

  _scene->mutex.lock();
  Shader *const shader = _scene->create_node<Shader>();
  _scene->mutex.unlock();
  array<Node *> used_shaders;
  used_shaders.push_back_slow(shader);
  _light->set_used_shaders(used_shaders);

  // Create default shader graph
  populate_shader_graph(light_inst, true);
}

void RIBCyclesLight::populate_shader_graph(Light_Scene_Entity &light_inst,
                                           bool initializing)
{
  unique_ptr<ShaderGraph> graph = make_unique<ShaderGraph>();
  ShaderNode *outputNode = nullptr;

  if (light_inst.light_type == "PxrDomeLight") {
    BackgroundNode *bgNode = graph->create_node<BackgroundNode>();
    // Bake strength into shader graph, since only the shader is used for background lights
    bgNode->set_color(_light->get_strength());

    graph->connect(bgNode->output("Background"), graph->output()->input("Surface"));

    outputNode = bgNode;
  }
  else {
    EmissionNode *emissionNode = graph->create_node<EmissionNode>();
    emissionNode->set_color(one_float3());
    emissionNode->set_strength(1.0f);

    graph->connect(emissionNode->output("Emission"), graph->output()->input("Surface"));

    outputNode = emissionNode;
  }

  bool hasSpatialVarying = false;
  bool hasColorTemperature = false;

  if (!initializing) {
    const bool enableColorTemperature = bool(
        light_inst.parameters.get_one_int("enableColorTemperature", 0));

    if (enableColorTemperature) {
      float temperature = light_inst.parameters.get_one_float("colorTemperature", 0);
      if (temperature > 0) {
        BlackbodyNode *blackbodyNode = graph->create_node<BlackbodyNode>();
        blackbodyNode->set_temperature(temperature);

        if (light_inst.light_type == "PxrDomeLight") {
          VectorMathNode *mathNode = graph->create_node<VectorMathNode>();
          mathNode->set_math_type(NODE_VECTOR_MATH_MULTIPLY);
          mathNode->set_vector2(_light->get_strength());

          graph->connect(blackbodyNode->output("Color"), mathNode->input("Vector1"));
          graph->connect(mathNode->output("Vector"), outputNode->input("Color"));
        }
        else {
          graph->connect(blackbodyNode->output("Color"), outputNode->input("Color"));
        }

        hasColorTemperature = true;
      }
    }

    std::string shapingIesFile = light_inst.parameters.get_one_string("shapingIesFile", "");
    if (!shapingIesFile.empty()) {

      TextureCoordinateNode *coordNode = graph->create_node<TextureCoordinateNode>();
      coordNode->set_ob_tfm(_instance->get_tfm());
      coordNode->set_use_transform(true);

      IESLightNode *iesNode = graph->create_node<IESLightNode>();
      iesNode->set_filename(ustring(shapingIesFile));

      graph->connect(coordNode->output("Normal"), iesNode->input("Vector"));
      graph->connect(iesNode->output("Fac"), outputNode->input("Strength"));

      hasSpatialVarying = true;
    }

    std::string texture_file = light_inst.parameters.get_one_string("lightColorMap", "");
    if (!texture_file.empty()) {

      ImageSlotTextureNode *textureNode = nullptr;
      if (light_inst.light_type == "PxrDomeLight") {
        Transform tfm = _instance->get_tfm();
        transform_set_column(&tfm, 3, zero_float3());  // Remove translation

        TextureCoordinateNode *coordNode = graph->create_node<TextureCoordinateNode>();
        coordNode->set_ob_tfm(tfm);
        coordNode->set_use_transform(true);

        MappingNode *mapper = graph->create_node<MappingNode>();
        mapper->set_rotation(make_float3(0, 0, M_PI_2));
        graph->connect(coordNode->output("Object"), mapper->input("Vector"));

        textureNode = graph->create_node<EnvironmentTextureNode>();
        static_cast<EnvironmentTextureNode *>(textureNode)->set_filename(ustring(texture_file));

        graph->connect(mapper->output("Vector"), textureNode->input("Vector"));

        hasSpatialVarying = true;
      }
      else {
        GeometryNode *coordNode = graph->create_node<GeometryNode>();

        textureNode = graph->create_node<ImageTextureNode>();
        static_cast<ImageTextureNode *>(textureNode)->set_filename(ustring(texture_file));

        graph->connect(coordNode->output("Parametric"), textureNode->input("Vector"));
      }

      if (hasColorTemperature) {
        VectorMathNode *mathNode = graph->create_node<VectorMathNode>();
        mathNode->set_math_type(NODE_VECTOR_MATH_MULTIPLY);

        graph->connect(textureNode->output("Color"), mathNode->input("Vector1"));
        ShaderInput *const outputNodeInput = outputNode->input("Color");
        graph->connect(outputNodeInput->link, mathNode->input("Vector2"));
        graph->disconnect(outputNodeInput);
        graph->connect(mathNode->output("Vector"), outputNodeInput);
      }
      else if (light_inst.light_type == "PxrDomeLight") {
        VectorMathNode *mathNode = graph->create_node<VectorMathNode>();
        mathNode->set_math_type(NODE_VECTOR_MATH_MULTIPLY);
        mathNode->set_vector2(_light->get_strength());

        graph->connect(textureNode->output("Color"), mathNode->input("Vector1"));
        graph->connect(mathNode->output("Vector"), outputNode->input("Color"));
      }
      else {
        graph->connect(textureNode->output("Color"), outputNode->input("Color"));
      }
    }
  }

  Shader *const shader = _light->get_shader();
  shader->set_graph(std::move(graph));
  shader->tag_update((Scene *)_light->get_owner());

  shader->has_surface_spatial_varying = hasSpatialVarying;
}

CCL_NAMESPACE_END
