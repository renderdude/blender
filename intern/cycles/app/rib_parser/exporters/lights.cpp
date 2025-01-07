#include "scene/light.h"
#include "scene/shader_graph.h"
#include "scene/shader_nodes.h"

#include "app/rib_parser/exporters/lights.h"
#include "util/transform.h"

CCL_NAMESPACE_BEGIN


Light *initialize(Scene *scene,
                  Light_Scene_Entity &light_inst,
                  Instance_Definition_Scene_Entity *inst_def);
void populate_shader_graph(Light_Scene_Entity &light_inst,
                           Light *light,
                           bool initializing = false);

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

void export_lights(Scene *scene,
                   vector<Instance_Scene_Entity> &inst_v,
                   Instance_Definition_Scene_Entity *inst_def)
{
  for (auto &light_inst : inst_def->lights) {
    for (auto &inst : inst_v) {
      const float metersPerUnit = 1.;
      ProjectionTransform xform_obj = *inst.render_from_instance * *light_inst.render_from_light;
      Transform xform = transform_scale(make_float3(metersPerUnit)) * convert_transform(xform_obj);
      vector<Transform> motion = {xform};
      vector<DecomposedTransform> decomp(motion.size());
      transform_motion_decompose(decomp.data(), motion.data(), motion.size());

      Light *light = initialize(scene, light_inst, inst_def);

      light->set_tfm(xform);

      float3 strength = make_float3(1.0f, 1.0f, 1.0f);

      auto color = light_inst.parameters.get_one_color("lightColor",
                                                       make_float3(1.0f, 1.0f, 1.0f));
      strength = make_float3(color[0], color[1], color[2]);

      float exposure = light_inst.parameters.get_one_float("exposure", 0.0);
      strength *= exp2(exposure);

      float intensity = light_inst.parameters.get_one_float("intensity", 1.0);
      strength *= intensity;

      // Cycles lights are normalized by default, so need to scale intensity if RMan light is not
      bool normalize = light_inst.parameters.get_one_int("areaNormalize", 0) == 1;
      light->set_normalize(normalize);

      auto &visibility = inst.parameters["visibility"];
      light->set_use_camera(bool(visibility.get_one_int("camera", 0)));
      // Default to shadow casting until we have an example
      light->set_cast_shadow(true);

      if (light_inst.light_type == "PxrDistantLight") {
        light->set_angle(OIIO::radians(light_inst.parameters.get_one_float("angle", 0.526f)));
      }
      else if (light_inst.light_type == "PxrDiskLight") {
        const float size = light_inst.parameters.get_one_float("size", 1.f) * 2.0f;
        light->set_sizeu(size);
        light->set_sizev(size);
      }
      else if (light_inst.light_type == "PxrRectLight") {
        // Size is set in the tfm
        #if 1
        light->set_sizeu(1.f);
        light->set_sizev(1.f);
        #else
        light->set_sizeu(1.f * fabsf(decomp[0].z.w));
        light->set_sizev(1.f * fabsf(decomp[0].z.w));
        #endif
      }
      else if (light_inst.light_type == "PxrSphereLight") {
        #if 1
        light->set_size(0.5f);
        #else
        light->set_size(0.5f * fabsf(decomp[0].z.w));
        #endif

        bool shaping = false;
        const float shapingConeAngle = light_inst.parameters.get_one_float("shapingConeAngle",
                                                                           -1.f);
        if (shapingConeAngle > 0) {
          light->set_spot_angle(OIIO::radians(shapingConeAngle * 2.0f));
          shaping = true;
        }

        const float shapingConeSoftness = light_inst.parameters.get_one_float(
            "shapingConeSoftness", -1.f);
        if (shapingConeSoftness > 0) {
          light->set_spot_smooth(shapingConeSoftness);
          shaping = true;
        }

        light->set_light_type(shaping ? LIGHT_SPOT : LIGHT_POINT);
      }

      light->set_strength(strength);
      light->set_is_enabled(true);

      populate_shader_graph(light_inst, light);

      if (light->is_modified())
        light->tag_update(scene);
    }
  }
}

Light *initialize(Scene *scene,
                  Light_Scene_Entity &light_inst,
                  Instance_Definition_Scene_Entity *inst_def)
{
  Light *light = scene->create_node<ccl::Light>();
  light->name = inst_def->name;

  light->set_random_id(hash_uint2(hash_string(light->name.c_str()), 0));

  if (light_inst.light_type == "PxrDomeLight") {
    light->set_light_type(LIGHT_BACKGROUND);
  }
  else if (light_inst.light_type == "PxrDistantLight") {
    light->set_light_type(LIGHT_DISTANT);
  }
  else if (light_inst.light_type == "PxrDiskLight") {
    light->set_light_type(LIGHT_AREA);
    light->set_ellipse(true);
    light->set_size(1.0f);
  }
  else if (light_inst.light_type == "PxrRectLight") {
    light->set_light_type(LIGHT_AREA);
    light->set_ellipse(false);
    light->set_size(1.0f);
  }
  else if (light_inst.light_type == "PxrSphereLight") {
    light->set_light_type(LIGHT_POINT);
    light->set_size(1.0f);
  }

  light->set_use_mis(true);
  light->set_use_camera(false);

  Shader *const shader = scene->create_node<Shader>();
  light->set_shader(shader);

  // Create default shader graph
  populate_shader_graph(light_inst, light, true);

  return light;
}

void populate_shader_graph(Light_Scene_Entity &light_inst, Light *light, bool initializing)
{
  unique_ptr<ShaderGraph> graph = make_unique<ShaderGraph>();
  ShaderNode *outputNode = nullptr;

  if (light_inst.light_type == "PxrDomeLight") {
    BackgroundNode *bgNode = graph->create_node<BackgroundNode>();
    // Bake strength into shader graph, since only the shader is used for background lights
    bgNode->set_color(light->get_strength());

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
          mathNode->set_vector2(light->get_strength());

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
      coordNode->set_ob_tfm(light->get_tfm());
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
        Transform tfm = light->get_tfm();
        transform_set_column(&tfm, 3, zero_float3());  // Remove translation

        TextureCoordinateNode *coordNode = graph->create_node<TextureCoordinateNode>();
        coordNode->set_ob_tfm(tfm);
        coordNode->set_use_transform(true);

        MappingNode* mapper = graph->create_node<MappingNode>();
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
        mathNode->set_vector2(light->get_strength());

        graph->connect(textureNode->output("Color"), mathNode->input("Vector1"));
        graph->connect(mathNode->output("Vector"), outputNode->input("Color"));
      }
      else {
        graph->connect(textureNode->output("Color"), outputNode->input("Color"));
      }
    }
  }

  Shader *const shader = light->get_shader();
  shader->set_graph(std::move(graph));
  shader->tag_update((Scene *)light->get_owner());

  shader->has_surface_spatial_varying = hasSpatialVarying;
}

CCL_NAMESPACE_END
