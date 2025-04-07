#include <cassert>
#include <cmath>
#include <fcntl.h>
#include <sstream>
#include <sys/stat.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <sys/sysctl.h>
#include <utility>
#include <vector>

#include <double-conversion/double-conversion.h>

#include "app/rib_parser/exporters/curves.h"
#include "app/rib_parser/exporters/geometry.h"
#ifdef WITH_CYCLES_DISTRIBUTED
#include "distributed/distributed.h"
#endif
#include "kernel/types.h"
#include "scene/camera.h"
#include "scene/object.h"
#include "scene/shader.h"
// Have to include shader.h before background.h so that 'set_shader' uses the correct 'set'
// overload taking a 'Node *', rather than the one taking a 'bool'
#include "scene/background.h"
#include "scene/light.h"
#include "scene/scene.h"
#include "scene/shader_graph.h"
#include "scene/shader_nodes.h"
#include "session/session.h"
#include "util/log.h"
#include "util/projection.h"
#include "util/transform.h"
#include "util/vector.h"

#include "exporters/lights.h"

#include "app/cycles_xml.h"
#include "app/rib_parser/param_dict.h"
#include "app/rib_parser/parsed_parameter.h"
#include "error.h"
#include "ri_api.h"

CCL_NAMESPACE_BEGIN

// Similar to util/projection.h with phi/theta reversed to match renderman specification
float3 ri_spherical_to_direction(float theta, float phi)
{
  float sin_phi = sinf(phi);
  return make_float3(sin_phi * cosf(theta), cosf(phi), sin_phi * sinf(theta));
}

std::vector<std::string> split_string(std::string_view str, char ch)
{
  std::vector<std::string> strings;

  if (str.empty()) {
    return strings;
  }

  std::string_view::iterator begin = str.begin();
  while (true) {
    std::string_view::iterator end = begin;
    while (end != str.end() && *end != ch) {
      ++end;
    }

    strings.emplace_back(begin, end);
    if (end == str.end()) {
      break;
    }

    begin = end + 1;
  }

  return strings;
}

template<typename T = std::mt19937> auto random_generator() -> T
{
  auto constexpr seed_bytes = sizeof(typename T::result_type) * T::state_size;
  auto constexpr seed_len = seed_bytes / sizeof(std::seed_seq::result_type);
  auto seed = std::array<std::seed_seq::result_type, seed_len>();
  auto dev = std::random_device();
  std::generate_n(begin(seed), seed_len, std::ref(dev));
  auto seed_seq = std::seed_seq(begin(seed), end(seed));
  return T{seed_seq};
}

auto generate_random_alphanumeric_string(std::size_t len = 8) -> std::string
{
  static constexpr auto chars =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";
  thread_local auto rng = random_generator<>();
  auto dist = std::uniform_int_distribution{{}, std::strlen(chars) - 1};
  auto result = std::string(len, '\0');
  std::generate_n(begin(result), len, [&]() { return chars[dist(rng)]; });
  return result;
}

Ri::Ri(Options &options) : session(options.session.get())
{
#ifdef WITH_CYCLES_DISTRIBUTED
  if (options.is_distributed) {
    Parsed_Parameter_Vector params;
    Parsed_Parameter *param;

    // Distributed class pointer
    param = new Parsed_Parameter(Parameter_Type::String, "class", File_Loc());
    param->payload = vector<void *>();
    param->add_pointer(options.session->distributed);
    params.push_back(param);

    // Remote Directory
    param = new Parsed_Parameter(Parameter_Type::String, "remote_directory", File_Loc());
    param->payload = vector<std::string>();
    param->add_string(options.directory);
    params.push_back(param);

    Option("distributed", params, File_Loc());

    init_request_thread();
  }
#endif
}

void Ri::adjust_buffer_parameters(BufferParams *buffer)
{
  Camera_Scene_Entity &camera = _camera[film.camera_name];
  int left = (int)(buffer->full_width * camera.crop_window[0]);
  int right = (int)(buffer->full_width * camera.crop_window[1]);
  int bottom = (int)(buffer->full_height * camera.crop_window[2]);
  int top = (int)(buffer->full_height * camera.crop_window[3]);

  buffer->width = right - left;
  buffer->window_width = buffer->width;
  buffer->height = top - bottom;
  buffer->window_height = buffer->height;
  buffer->full_x = left;
  buffer->full_y = bottom;
}

void Ri::export_to_cycles()
{
  BoundBox scene_bounds{BoundBox::empty};

  export_options(filter, film, _camera[film.camera_name], sampler);

  // Force the instance jobs to finish
  for (auto *instance_job : _instance_use_jobs) {
    instance_job->get_result();
  }
  _instance_use_jobs.clear();

  // that, should force the instance definitions and shders to finish, so we should
  // be good to clear them here
  {
    std::lock_guard<std::mutex> lock(shader_mutex);
    for (auto shader : _shader_jobs) {
      RIBCyclesMaterials material = shader.second->get_result();
    }
  }
  _shader_jobs.clear();

  // Lifted from hydra/session.cpp
  Scene *const scene = session->scene.get();

  // Update background depending on presence of a background light
  if (scene->light_manager->need_update()) {
    ccl::Light *background_light = nullptr;
    bool have_lights = false;
    for (Object *object : scene->objects) {
      if (!object->get_geometry()->is_light()) {
        continue;
      }

      have_lights = true;

      ccl::Light *light = static_cast<ccl::Light *>(object->get_geometry());
      if (light->get_light_type() == LIGHT_BACKGROUND) {
        background_light = light;
        break;
      }
    }

    if (!background_light) {
      scene->background->set_shader(scene->default_background);
      scene->background->set_transparent(true);

      /* Set background color depending to non-zero value if there are no
       * lights in the scene, to match behavior of other renderers. */
      for (ShaderNode *node : scene->default_background->graph->nodes) {
        if (node->is_a(BackgroundNode::get_node_type())) {
          BackgroundNode *bgNode = static_cast<BackgroundNode *>(node);
          bgNode->set_color((have_lights) ? zero_float3() : make_float3(0.5f));
        }
      }
    }
    else {
      scene->background->set_shader(background_light->get_shader());
      scene->background->set_transparent(false);
    }

    scene->background->tag_update(scene);
  }
}

inline float radians(float deg)
{
  return (M_PI_F / 180) * deg;
}

inline float degrees(float rad)
{
  return (180 / M_PI_F) * rad;
}

void Ri::export_options([[maybe_unused]] Scene_Entity &filter,
                        Scene_Entity &film,
                        Camera_Scene_Entity &camera,
                        [[maybe_unused]] Scene_Entity &sampler)
{
  // Immediately create filter and film
  VLOG(1) << "Starting to create filter and film";
  // Do something with the filter
  // Filter filt = Filter::create(filter.name, filter.parameters, &filter.loc, alloc);

  // It's a little ugly to poke into the camera's parameters here, but we
  // have this circular dependency that Camera::create() expects a
  // Film, yet now the film needs to know the exposure time from
  // the camera....
  float exposureTime = camera.parameters.get_one_float("shutterclose", 1.f) -
                       camera.parameters.get_one_float("shutteropen", 0.f);
  if (exposureTime <= 0) {
    error_exit(&camera.loc,
               "The specified camera shutter times imply that the shutter "
               "does not open.  A black image will result.");
  }

  _display_name = film.parameters.get_one_string("filename", "");

  Camera *cam = session->scene->camera;

  int x_res = film.parameters.get_one_int("xresolution", 1280);
  int y_res = film.parameters.get_one_int("yresolution", 720);
  cam->set_full_width(x_res);
  cam->set_full_height(y_res);
  // cam->set_screen_size(x_res, y_res);

  const float metersPerUnit = 1.;

  Transform t = projection_to_transform(camera.camera_transform.render_from_camera());
  t.x.w *= metersPerUnit;
  t.y.w *= metersPerUnit;
  t.z.w *= metersPerUnit;

  cam->set_matrix(t);
  float near = camera.parameters.get_one_float("nearClip", -1.f);
  if (near > 0) {
    cam->set_nearclip(near * metersPerUnit);
  }
  float far = camera.parameters.get_one_float("farClip", -1.f);
  if (far >= 0) {
    cam->set_farclip(far * metersPerUnit);
  }

  // Set FOV
  float fov = camera.parameters.get_one_float("fov", 45.f);
  cam->set_fov(radians(fov));

  // Set Depth of Field
  float fstop = camera.parameters.get_one_float("fStop", 9.999e37);
  if (fstop < 1e30) {
    float focaldistance = camera.parameters.get_one_float("focalDistance", 0);
    float focallength = camera.parameters.get_one_float("focalLength", 0);
    int blades = camera.parameters.get_one_int("apertureNSides", 0);
    float apertureAngle = camera.parameters.get_one_float("apertureAngle", 0);
    cam->set_focaldistance(focaldistance);
    cam->set_aperturesize(focallength / (2.0 * fstop));
    cam->set_blades(blades);
    cam->set_bladesrotation(apertureAngle);
  }

  // Set screen window
  auto screen_window = camera.parameters.get_float_array("ScreenWindow");
  cam->set_viewplane_left(screen_window[0]);
  cam->set_viewplane_right(screen_window[1]);
  cam->set_viewplane_bottom(screen_window[2]);
  cam->set_viewplane_top(screen_window[3]);

  cam->need_flags_update = true;
  cam->update(session->scene.get());
  *session->scene->dicing_camera = *session->scene->camera;
}

int Ri::add_area_light(Scene_Entity light)
{
  std::lock_guard<std::mutex> lock(area_light_mutex);
  area_lights.push_back(std::move(light));
  return area_lights.size() - 1;
}

void Ri::add_animated_shape(Animated_Shape_Scene_Entity shape)
{
  std::lock_guard<std::mutex> lock(animated_shape_mutex);
  animated_shapes.push_back(std::move(shape));
}

void Ri::add_instance_definition(Instance_Definition_Scene_Entity instance)
{
  Instance_Definition_Scene_Entity *def = new Instance_Definition_Scene_Entity(
      std::move(instance));

  if (def->lights.empty()) {
    if (def->shapes[0].name == "curve") {
      auto create = [this](Instance_Definition_Scene_Entity *def) {
        RIBCyclesCurves *curve = new RIBCyclesCurves(session->scene.get());
        curve->build_instance_definition(def);
        return curve;
      };
      _curve_definition_jobs[def->name] = run_async(create, def);
    }
    else {
      auto create = [this](Instance_Definition_Scene_Entity *def) {
        RIBCyclesMesh *mesh = new RIBCyclesMesh(session->scene.get());
        mesh->build_instance_definition(def);
        delete def;
        return mesh;
      };
      _mesh_definition_jobs[def->name] = run_async(create, def);
    }
  }
  else {
    std::lock_guard<std::mutex> lock(instance_definition_mutex);
    instance_definitions[def->name] = def;
  }
}

void Ri::add_instance_use(Instance_Scene_Entity instance)
{
  Instance_Scene_Entity *inst = new Instance_Scene_Entity(std::move(instance));
  auto create = [this](Instance_Scene_Entity *inst) {
    if (!inst->material_name.empty()) {
      // Wait until the material has been added as a shader job
      bool exists = false;
      while (!exists) {
        if (_shader_jobs.find(inst->material_name) != _shader_jobs.end()) {
          RIBCyclesMaterials material = _shader_jobs[inst->material_name]->get_result();
          exists = true;
        }
      }
    }
    if (instance_definitions.find(inst->name) != instance_definitions.end()) {
      auto *inst_def = instance_definitions[inst->name];
      RIBCyclesLight light(session->scene.get());
      light.build(inst_def, *inst);
    }
    // it's a shape
    else {
      if (_mesh_definition_jobs.find(inst->name) != _mesh_definition_jobs.end()) {
        RIBCyclesMesh *mesh = _mesh_definition_jobs[inst->name]->get_result();
        _mesh_elements[inst->name]["unassigned"] = mesh;
      }
      if (_mesh_elements.find(inst->name) != _mesh_elements.end()) {
        auto mapped_meshes = _mesh_elements[inst->name];
        RIBCyclesMesh *mesh = nullptr;
        if (mapped_meshes.find("unassigned") != mapped_meshes.end()) {
          std::lock_guard<std::mutex> lock(instance_definition_mutex);
          // Processing first instance of the definition, so remove it from the list
          mesh = mapped_meshes["unassigned"];
          mapped_meshes.erase("unassigned");
        }
        else {
          mesh = mapped_meshes[inst->material_name];
          array<Node *> geom_shaders = mesh->get_used_shaders();
          if (geom_shaders.size() > 0 && !geom_shaders[0]->name.compare(inst->material_name)) {
            RIBCyclesMesh *new_mesh = new RIBCyclesMesh(*mesh);
            mesh = new_mesh;
          }
        }

        mesh->build_instance(*inst);
        if (mapped_meshes.find(mesh->get_used_shaders()[0]->name.string()) == mapped_meshes.end())
        {
          std::lock_guard<std::mutex> lock(instance_definition_mutex);
          mapped_meshes[mesh->get_used_shaders()[0]->name.string()] = mesh;
          assert(mesh->get_used_shaders()[0]->name.string() == inst->material_name);
        }
      }
      else {
        assert(_curve_elements.find(inst->name) != _curve_elements.end());
        if (_curve_definition_jobs.find(inst->name) != _curve_definition_jobs.end()) {
          RIBCyclesCurves *mesh = _curve_definition_jobs[inst->name]->get_result();
          _curve_elements[inst->name]["unassigned"] = mesh;
        }
        if (_curve_elements.find(inst->name) != _curve_elements.end()) {
          auto mapped_curves = _curve_elements[inst->name];
          RIBCyclesCurves *curve = nullptr;
          if (mapped_curves.find("unassigned") != mapped_curves.end()) {
            std::lock_guard<std::mutex> lock(instance_definition_mutex);
            // Processing first instance of the definition, so remove it from the list
            curve = mapped_curves["unassigned"];
            mapped_curves.erase("unassigned");
          }
          else {
            array<Node *> geom_shaders = curve->get_used_shaders();
            if (geom_shaders.size() > 0 && !geom_shaders[0]->name.compare(inst->material_name)) {
              RIBCyclesCurves *new_curve = new RIBCyclesCurves(*curve);
              curve = new_curve;
            }
          }

          curve->build_instance(*inst);
          if (mapped_curves.find(curve->get_used_shaders()[0]->name.string()) ==
              mapped_curves.end())
          {
            std::lock_guard<std::mutex> lock(instance_definition_mutex);
            mapped_curves[curve->get_used_shaders()[0]->name.string()] = curve;
            assert(curve->get_used_shaders()[0]->name.string() == inst->material_name);
          }
        }
      }
    }
    return true;
  };

  _instance_use_jobs.push_back(run_async(create, inst));
}

void Ri::add_shader(Vector_Dictionary shader)
{
  std::lock_guard<std::mutex> lock(shader_mutex);

  auto create = [this](Vector_Dictionary shader) {
    RIBCyclesMaterials material(session->scene.get(), shader);
    material.export_materials();
    return material;
  };
  _shader_jobs[shader.first] = run_async(create, shader);
}

// RI API Default Implementation
void Ri::ArchiveBegin([[maybe_unused]] const std::string &name,
                      [[maybe_unused]] Parsed_Parameter_Vector params,
                      [[maybe_unused]] File_Loc loc)
{
  std::cout << "ArchiveBegin is unimplemented" << std::endl;
}

void Ri::ArchiveEnd([[maybe_unused]] File_Loc loc)
{
  std::cout << "ArchiveEnd is unimplemented" << std::endl;
}

void Ri::AreaLightSource([[maybe_unused]] const std::string &name,
                         [[maybe_unused]] Parsed_Parameter_Vector params,
                         [[maybe_unused]] File_Loc loc)
{
  std::cout << "AreaLightSource is unimplemented." << std::endl;
}

void Ri::Atmosphere([[maybe_unused]] const std::string &name,
                    [[maybe_unused]] Parsed_Parameter_Vector params,
                    [[maybe_unused]] File_Loc loc)
{
  std::cout << "Atmosphere is unimplemented" << std::endl;
}

void Ri::Attribute(const std::string &target,
                   Parsed_Parameter_Vector params,
                   [[maybe_unused]] File_Loc loc)
{
  for (Parsed_Parameter *p : params) {
    p->may_be_unused = true;
    graphics_state.rib_attributes[target].push_back(p);
  }

  if (target == "Ri") {
    for (Parsed_Parameter *p : params) {
      if (p->name == "Orientation") {
        if (p->strings()[0] == "inside") {
          graphics_state.reverse_orientation = !graphics_state.reverse_orientation;
        }
      }
    }
  }
  else if (target == "identifier") {
    for (Parsed_Parameter *p : params) {
      if (p->name == "name") {
        graphics_state.id_string = p->strings()[0];
      }
      else if (p->name == "id") {
        graphics_state.identifier = p->ints()[0];
      }
    }
  }
}

void Ri::AttributeBegin(File_Loc loc)
{
  VERIFY_WORLD("AttributeBegin");
  osl_parameters.clear();
  pushed_graphics_states.push_back(graphics_state);
  push_stack.emplace_back('a', loc);
}

void Ri::AttributeEnd(File_Loc loc)
{
  VERIFY_WORLD("AttributeEnd");
  if (!osl_parameters.empty()) {
    osl_shader_group[_shader_id] = true;
    add_shader(std::make_pair(_shader_id, std::move(osl_parameters)));
  }

  // Issue error on unmatched _AttributeEnd_
  if (pushed_graphics_states.empty()) {
    error(&loc, "Unmatched AttributeEnd encountered. Ignoring it.");
    return;
  }

  // NOTE: must keep the following consistent with code in ObjectEnd
  graphics_state = std::move(pushed_graphics_states.back());
  pushed_graphics_states.pop_back();

  if (push_stack.back().first == 'o') {
    std::stringstream ss;
    ss << "Mismatched nesting: open ObjectBegin from ";
    ss << push_stack.back().second.to_string();
    ss << " at AttributeEnd";
    error_exit_deferred(&loc, ss.str());
  }
  else {
    CHECK_EQ(push_stack.back().first, 'a');
  }
  push_stack.pop_back();
}

void Ri::Basis([[maybe_unused]] Parsed_Parameter_Vector params, [[maybe_unused]] File_Loc loc)
{
  std::cout << "Basis is unimplemented" << std::endl;
}

void Ri::Begin([[maybe_unused]] const std::string &name, [[maybe_unused]] File_Loc loc)
{
  std::cout << "Begin is unimplemented" << std::endl;
}

void Ri::Bound([[maybe_unused]] float bound[6], [[maybe_unused]] File_Loc loc)
{
  std::cout << "Bound is unimplemented" << std::endl;
}

void Ri::Bxdf(const std::string &bxdf,
              const std::string &name,
              Parsed_Parameter_Vector params,
              File_Loc loc)
{
  VERIFY_WORLD("Bxdf");

  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::String, "shader_type", loc);
  param->add_string("bxdf");
  param->add_string(bxdf);
  param->add_string(name);
  params.push_back(param);

  Parameter_Dictionary dict(std::move(params));
  std::string material_id = dict.get_one_string("__materialid", "");
  if (graphics_state.remapped_material.find(material_id) != graphics_state.remapped_material.end())
  {
    material_id = graphics_state.remapped_material[material_id];
    param = dict.get_parameter("__materialid");
    param->strings()[0] = material_id;
  }
  if (!material_id.empty()) {
    _shader_id = material_id;
    if (osl_shader_group.find(material_id) == osl_shader_group.end()) {
      osl_parameters.push_back(dict);
    }
    else {
      graphics_state.current_material_name = material_id;
      graphics_state.current_material_index = -1;
    }
  }
  else if (string_startswith(bxdf, "Lama")) {
    osl_parameters.push_back(dict);
  }
}

void Ri::camera(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  VERIFY_OPTIONS("Camera");
  // Remove any class designator from the camera options
  // Specifically, Ri:<something>
  for (auto it = params.begin(); it != params.end(); it++) {
    std::vector<std::string> strings = split_string((*it)->name, ':');
    if (strings.size() > 1) {
      // Could extract the size and check we have that many values below
      (*it)->name = strings[1];
    }
  }

  auto options = _rib_state.options["Ri"];

  for (auto s : options) {
    auto *param = s.second;
    bool found = false;
    for (auto it = params.begin(); it != params.end(); it++) {
      if ((*it)->name == param->name) {
        found = true;
        break;
      }
    }

    if (!found) {
      params.push_back(param);
    }
  }

  Parameter_Dictionary dict(std::move(params));

  auto items = _rib_state.options.find("trace");
  if (items != _rib_state.options.end()) {
    auto world_origin = items->second.find("worldorigin");
    if (world_origin != items->second.end()) {
      if (world_origin->second->strings()[0] == "worldoffset") {
        auto *world_offset = items->second["worldoffset"];
        _rib_state.world_offset = make_float3(
            -world_offset->floats()[0], -world_offset->floats()[1], -world_offset->floats()[2]);
      }
    }
  }

  Transform_Set camera_from_world = graphics_state.ctm;
  Transform_Set world_from_camera = inverse(graphics_state.ctm);
  named_coordinate_systems["camera"] = inverse(camera_from_world);

  // Camera motion
  Camera_Transform camera_transform(world_from_camera[0], _rib_state.world_offset);
  render_from_world = camera_transform.render_from_world();

  _camera_name = name;
  _camera[name] = Camera_Scene_Entity("perspective",
                                      std::move(dict),
                                      loc,
                                      camera_transform,
                                      graphics_state.current_outside_medium);
  _camera[name].crop_window = _crop_window;
}

void Ri::Clipping([[maybe_unused]] float cnear,
                  [[maybe_unused]] float cfar,
                  [[maybe_unused]] File_Loc loc)
{
  std::cout << "Clipping is unimplemented" << std::endl;
}

void Ri::ClippingPlane([[maybe_unused]] float x,
                       [[maybe_unused]] float y,
                       [[maybe_unused]] float z,
                       [[maybe_unused]] float nx,
                       [[maybe_unused]] float ny,
                       [[maybe_unused]] float nz,
                       [[maybe_unused]] File_Loc loc)
{
  std::cout << "ClippingPlane is unimplemented" << std::endl;
}

void Ri::Color([[maybe_unused]] float r,
               [[maybe_unused]] float g,
               [[maybe_unused]] float b,
               [[maybe_unused]] File_Loc loc)
{
  std::cout << "Color is unimplemented" << std::endl;
}

void Ri::ConcatTransform([[maybe_unused]] const float transform[16], [[maybe_unused]] File_Loc loc)
{
  graphics_state.for_active_transforms([=](auto t) {
    ProjectionTransform projection = *(ProjectionTransform *)&transform[0];
    return t * projection_transpose(projection);
  });
}

void Ri::CoordinateSystem(std::string const &name, [[maybe_unused]] File_Loc loc)
{
  named_coordinate_systems[name] = graphics_state.ctm;
}

void Ri::CoordSysTransform(std::string const &name, File_Loc loc)
{
  if (named_coordinate_systems.find(name) != named_coordinate_systems.end()) {
    graphics_state.ctm = named_coordinate_systems[name];
  }
  else {
    std::stringstream ss;
    ss << "Couldn't find named coordinate system \"" << name << "\"";
    warning(&loc, ss.str());
  }
}

void Ri::Cone(
    float height, float radius, float thetamax, Parsed_Parameter_Vector params, File_Loc loc)
{
  thetamax = (thetamax > 360.0f ? 360.0f : (thetamax < 0.f ? 0.f : thetamax));
  int offset = 1;
  thetamax = thetamax * M_PI_F / 180.0f;

  std::vector<float> pts;
  std::vector<float> uvs;

  int n_slices = 25;

  for (int j = 0; j < n_slices + offset; ++j) {
    float theta = thetamax * float(j) / float(n_slices);
    float x = radius * std::cosf(theta);
    float y = radius * std::sinf(theta);
    pts.push_back(x);
    pts.push_back(y);
    pts.push_back(0.f);
    uvs.push_back(theta / thetamax);
    uvs.push_back(0.f);
    pts.push_back(0.f);
    pts.push_back(0.f);
    pts.push_back(height);
    uvs.push_back(theta / thetamax);
    uvs.push_back(1.f);
  }

  int poly_count = 0;
  std::vector<int> polys;

  for (int j = 0; j < 2 * n_slices; j += 2) {
    int jj = (j + 1) % (2 * n_slices + 1);
    polys.push_back(j);
    polys.push_back(jj);
    polys.push_back((jj + 2) % (2 * (n_slices + offset)));
    polys.push_back((j + 2) % (2 * (n_slices + offset)));
    poly_count++;
  }

  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::Integer, "vertices", loc);
  for (int i = 0; i < polys.size(); ++i) {
    param->add_int(polys[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nvertices", loc);
  for (int i = 0; i < poly_count; ++i) {
    param->add_int(4);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Point3, "P", loc);
  param->storage = Container_Type::Vertex;
  param->elem_per_item = 3;
  for (int i = 0; i < pts.size(); ++i) {
    param->add_float(pts[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Point2, "uv", loc);
  param->storage = Container_Type::Vertex;
  param->elem_per_item = 2;
  for (int i = 0; i < uvs.size(); ++i) {
    param->add_float(uvs[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nfaces", loc);
  param->add_int(poly_count);
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Boolean, "smooth", loc);
  param->add_bool(true);
  params.push_back(param);

  Shape("mesh", params, loc);
}

void Ri::CropWindow(float xmin, float xmax, float ymin, float ymax, [[maybe_unused]] File_Loc loc)
{
  if (xmin >= 0. && xmin <= 1. && xmin < xmax && xmin >= _crop_window[0] && xmax >= 0. &&
      xmax <= 1. && xmin < xmax && xmax <= _crop_window[1] && ymin >= 0. && ymin <= 1. &&
      ymin < ymax && ymin >= _crop_window[2] && ymax >= 0. && ymax <= 1. && ymin < ymax &&
      ymax <= _crop_window[3])
  {
    _crop_window[0] = xmin;
    _crop_window[1] = xmax;
    _crop_window[2] = ymin;
    _crop_window[3] = ymax;

    if (!_camera_name.empty()) {
      _camera[_camera_name].crop_window = _crop_window;
    }
  }
}

void Ri::Curves(const std::string &type,
                std::vector<int> nvertices,
                const std::string &wrap,
                Parsed_Parameter_Vector params,
                File_Loc loc)
{
  VERIFY_WORLD("Shape");

  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::Integer, "nvertices", loc);
  for (int nv : nvertices) {
    param->add_int(nv);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::String, "wrap", loc);
  param->add_string(wrap);
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::String, "type", loc);
  param->add_string(type);
  params.push_back(param);

  Shape("curve", params, loc);
}

void Ri::Cylinder(float radius,
                  float zmin,
                  float zmax,
                  float thetamax,
                  Parsed_Parameter_Vector params,
                  File_Loc loc)
{
  VERIFY_WORLD("Shape");

  thetamax = (thetamax > 360.0f ? 360.0f : (thetamax < 0.f ? 0.f : thetamax));
  int offset = 1;
  thetamax = thetamax * M_PI_F / 180.0f;

  std::vector<float> pts;
  std::vector<float> uvs;

  int n_slices = 25;

  for (int j = 0; j < n_slices + offset; ++j) {
    float theta = thetamax * float(j) / float(n_slices);
    float x = radius * std::cosf(theta);
    float y = radius * std::sinf(theta);
    pts.push_back(x);
    pts.push_back(y);
    pts.push_back(zmin);
    uvs.push_back(theta / thetamax);
    uvs.push_back(0.f);
    pts.push_back(x);
    pts.push_back(y);
    pts.push_back(zmax);
    uvs.push_back(theta / thetamax);
    uvs.push_back(1.f);
  }

  int poly_count = 0;
  std::vector<int> polys;

  for (int j = 0; j < 2 * n_slices; j += 2) {
    int jj = (j + 1) % (2 * n_slices + 1);
    polys.push_back(j);
    polys.push_back(jj);
    polys.push_back((jj + 2) % (2 * (n_slices + offset)));
    polys.push_back((j + 2) % (2 * (n_slices + offset)));
    poly_count++;
  }

  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::Integer, "vertices", loc);
  for (int i = 0; i < polys.size(); ++i) {
    param->add_int(polys[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nvertices", loc);
  for (int i = 0; i < poly_count; ++i) {
    param->add_int(4);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Point3, "P", loc);
  param->storage = Container_Type::Vertex;
  param->elem_per_item = 3;
  for (int i = 0; i < pts.size(); ++i) {
    param->add_float(pts[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Point2, "uv", loc);
  param->storage = Container_Type::Vertex;
  param->elem_per_item = 2;
  for (int i = 0; i < uvs.size(); ++i) {
    param->add_float(uvs[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nfaces", loc);
  param->add_int(poly_count);
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Boolean, "smooth", loc);
  param->add_bool(true);
  params.push_back(param);

  Shape("mesh", params, loc);
}

void Ri::Declare([[maybe_unused]] const std::string &name,
                 [[maybe_unused]] const std::string &declaration,
                 [[maybe_unused]] File_Loc loc)
{
  std::cout << "Declare is unimplemented" << std::endl;
}

void Ri::DepthOfField([[maybe_unused]] float fstop,
                      [[maybe_unused]] float focallength,
                      [[maybe_unused]] float focaldistance,
                      [[maybe_unused]] File_Loc loc)
{
  std::cout << "DepthOfField is unimplemented" << std::endl;
}

void Ri::Detail([[maybe_unused]] float bound[6], [[maybe_unused]] File_Loc loc)
{
  std::cout << "Detail is unimplemented" << std::endl;
}

void Ri::DetailRange([[maybe_unused]] float offlow,
                     [[maybe_unused]] float onlow,
                     [[maybe_unused]] float onhigh,
                     [[maybe_unused]] float offhigh,
                     [[maybe_unused]] File_Loc loc)
{
  std::cout << "DetailRange is unimplemented" << std::endl;
}

void Ri::Disk(
    float height, float radius, float thetamax, Parsed_Parameter_Vector params, File_Loc loc)
{
  VERIFY_WORLD("Shape");

  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::Real, "radius", loc);
  param->add_float(radius);
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Real, "height", loc);
  param->add_float(height);
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Real, "phimax", loc);
  param->add_float(thetamax);
  params.push_back(param);

  Shape("disk", params, loc);
}

void Ri::Displacement([[maybe_unused]] const std::string &displace,
                      [[maybe_unused]] const std::string &name,
                      [[maybe_unused]] Parsed_Parameter_Vector params,
                      [[maybe_unused]] File_Loc loc)
{
  std::cout << "Displacement is unimplemented" << std::endl;
}

void Ri::Display([[maybe_unused]] const std::string &name,
                 [[maybe_unused]] const std::string &type,
                 [[maybe_unused]] const std::string &mode,
                 [[maybe_unused]] Parsed_Parameter_Vector params,
                 [[maybe_unused]] File_Loc loc)
{
  VERIFY_OPTIONS("Film");
  // If the first char of `name' is a '+', then it's an additional
  // channel to render. Ignore it for now.
  if (name[0] == '+') {
    fprintf(stdout, "Ignoring additional display channel, %s\n", name.c_str());
    return;
  }

  Parsed_Parameter_Vector new_params;

  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::String, "filename", loc);
  param->add_string(name);
  new_params.push_back(param);

  auto resolution = _rib_state.options["Ri"]["FormatResolution"]->ints();
  param = new Parsed_Parameter(Parameter_Type::Integer, "xresolution", loc);
  param->add_int((int)(resolution[0]));
  new_params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "yresolution", loc);
  param->add_int((int)(resolution[1]));
  new_params.push_back(param);

  Parameter_Dictionary dict(std::move(new_params));
  std::string camera_name = dict.get_one_string("camera_name", "");
  if (camera_name.empty()) {
    camera_name = _camera_name;
  }
  else {
    _camera_name = camera_name;
  }

  film = Display_Scene_Entity("rgb", std::move(dict), loc, camera_name);
}

void Ri::DisplayChannel([[maybe_unused]] const std::string &name,
                        [[maybe_unused]] Parsed_Parameter_Vector params,
                        [[maybe_unused]] File_Loc loc)
{
  std::cout << "DisplayChannel " << name << " is unimplemented" << std::endl;
}

void Ri::DisplayFilter([[maybe_unused]] const std::string &name,
                       [[maybe_unused]] const std::string &type,
                       [[maybe_unused]] Parsed_Parameter_Vector params,
                       [[maybe_unused]] File_Loc loc)
{
  std::cout << "DisplayFilter " << name << " is unimplemented" << std::endl;
}

void Ri::Else([[maybe_unused]] File_Loc loc)
{
  std::cout << "Else is unimplemented" << std::endl;
}

void Ri::ElseIf([[maybe_unused]] const std::string &condition, [[maybe_unused]] File_Loc loc)
{
  std::cout << "ElseIf is unimplemented" << std::endl;
}

void Ri::End([[maybe_unused]] File_Loc loc)
{
  std::cout << "End is unimplemented" << std::endl;
}

void Ri::Exposure([[maybe_unused]] float gain,
                  [[maybe_unused]] float gamma,
                  [[maybe_unused]] File_Loc loc)
{
  std::cout << "Exposure is unimplemented" << std::endl;
}

void Ri::Exterior([[maybe_unused]] const std::string &name,
                  [[maybe_unused]] Parsed_Parameter_Vector params,
                  [[maybe_unused]] File_Loc loc)
{
  std::cout << "Exterior is unimplemented" << std::endl;
}

void Ri::Format([[maybe_unused]] int xresolution,
                [[maybe_unused]] int yresolution,
                [[maybe_unused]] float pixelaspectratio,
                [[maybe_unused]] File_Loc loc)
{
  std::cout << "Format is unimplemented" << std::endl;
}

void Ri::FrameAspectRatio([[maybe_unused]] float frameratio, [[maybe_unused]] File_Loc loc)
{
  std::cout << "FrameAspectRatio is unimplemented" << std::endl;
}

void Ri::FrameBegin([[maybe_unused]] int number, [[maybe_unused]] File_Loc loc)
{
  std::cout << "FrameBegin is unimplemented" << std::endl;
}

void Ri::FrameEnd([[maybe_unused]] File_Loc loc)
{
  std::cout << "FrameEnd is unimplemented" << std::endl;
}

void Ri::GeneralPolygon([[maybe_unused]] std::vector<int> nvertices,
                        [[maybe_unused]] Parsed_Parameter_Vector params,
                        [[maybe_unused]] File_Loc loc)
{
  std::cout << "GeneralPolygon is unimplemented" << std::endl;
}

void Ri::GeometricApproximation([[maybe_unused]] const std::string &type,
                                [[maybe_unused]] float value,
                                [[maybe_unused]] File_Loc loc)
{
  std::cout << "GeometricApproximation is unimplemented" << std::endl;
}

void Ri::Geometry([[maybe_unused]] const std::string &type,
                  [[maybe_unused]] Parsed_Parameter_Vector params,
                  [[maybe_unused]] File_Loc loc)
{
  std::cout << "Geometry is unimplemented" << std::endl;
}

void Ri::Hider([[maybe_unused]] const std::string &name,
               [[maybe_unused]] Parsed_Parameter_Vector params,
               [[maybe_unused]] File_Loc loc)
{
  std::cout << "Hider is unimplemented" << std::endl;
}

void Ri::HierarchicalSubdivisionMesh(const std::string &scheme,
                                     std::vector<int> n_vertices,
                                     std::vector<int> vertices,
                                     std::vector<std::string> tags,
                                     std::vector<int> nargs,
                                     std::vector<int> intargs,
                                     std::vector<float> floatargs,
                                     std::vector<std::string> stringargs,
                                     Parsed_Parameter_Vector params,
                                     File_Loc loc)
{
  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::String, "scheme", loc);
  param->add_string(scheme);
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "vertices", loc);
  for (int i = 0; i < vertices.size(); ++i) {
    param->add_int(vertices[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nvertices", loc);
  for (int i = 0; i < n_vertices.size(); ++i) {
    param->add_int(n_vertices[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nfaces", loc);
  param->add_int(n_vertices.size());
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::String, "tags", loc);
  for (int i = 0; i < tags.size(); ++i) {
    param->add_string(tags[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nargs", loc);
  for (int i = 0; i < nargs.size(); ++i) {
    param->add_int(nargs[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "intargs", loc);
  for (int i = 0; i < intargs.size(); ++i) {
    param->add_int(intargs[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Real, "floatargs", loc);
  for (int i = 0; i < floatargs.size(); ++i) {
    param->add_float(floatargs[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Real, "stringargs", loc);
  for (int i = 0; i < stringargs.size(); ++i) {
    param->add_string(stringargs[i]);
  }
  params.push_back(param);

  Shape("subdivision_mesh", params, loc);
}
void Ri::Hyperboloid([[maybe_unused]] Point3f point1,
                     [[maybe_unused]] Point3f point2,
                     [[maybe_unused]] float thetamax,
                     [[maybe_unused]] Parsed_Parameter_Vector params,
                     [[maybe_unused]] File_Loc loc)
{
  std::cout << "Hyperboloid is unimplemented" << std::endl;
}

void Ri::Identity([[maybe_unused]] File_Loc loc)
{
  graphics_state.for_active_transforms(
      []([[maybe_unused]] auto t) { return projection_identity(); });
}

void Ri::IfBegin([[maybe_unused]] const std::string &condition, [[maybe_unused]] File_Loc loc)
{
  std::cout << "IfBegin is unimplemented" << std::endl;
}

void Ri::IfEnd([[maybe_unused]] File_Loc loc)
{
  std::cout << "IfEnd is unimplemented" << std::endl;
}

void Ri::Illuminate([[maybe_unused]] const std::string &light,
                    [[maybe_unused]] bool onoff,
                    [[maybe_unused]] File_Loc loc)
{
  std::cout << "Illuminate is unimplemented" << std::endl;
}

void Ri::Imager([[maybe_unused]] const std::string &name,
                [[maybe_unused]] Parsed_Parameter_Vector params,
                [[maybe_unused]] File_Loc loc)
{
  std::cout << "Imager is unimplemented" << std::endl;
}

void Ri::Integrator([[maybe_unused]] const std::string &type,
                    [[maybe_unused]] const std::string &name,
                    [[maybe_unused]] Parsed_Parameter_Vector params,
                    [[maybe_unused]] File_Loc loc)
{
  std::cout << "Integrator is unimplemented" << std::endl;
}

void Ri::Interior([[maybe_unused]] const std::string &name,
                  [[maybe_unused]] Parsed_Parameter_Vector params,
                  [[maybe_unused]] File_Loc loc)
{
  std::cout << "Interior is unimplemented" << std::endl;
}

void Ri::Light(const std::string &name,
               const std::string &handle,
               Parsed_Parameter_Vector params,
               File_Loc loc)
{
  VERIFY_WORLD("Light");
  Parsed_Parameter_Vector light_params;

  bool double_sided = false;
  for (auto it = graphics_state.shape_attributes.begin();
       it != graphics_state.shape_attributes.end();
       it++)
  {
    if ((*it)->name == "sides") {
      double_sided = (*it)->floats()[0] == 2;
      break;
    }
  }

  Parsed_Parameter *param;
  // Lights are single-sided by default
  if (double_sided) {
    param = new Parsed_Parameter(Parameter_Type::Boolean, "twosided", loc);
    param->may_be_unused = true;
    param->add_bool(true);
    params.push_back(param);
  }

  ProjectionTransform const *render_from_object = transform_cache.lookup(Render_From_Object(0));

  Light_Scene_Entity entity(
      handle,
      name,
      Parameter_Dictionary(std::move(params), graphics_state.light_attributes),
      loc,
      render_from_object);

  if (active_instance_definition) {
    if (name == "PxrCylinderLight" || name == "PxrSphereLight") {
      if (name == "PxrCylinderLight") {
        // A cylinder light lies along the X-axis, and is half the default cylinder size
        Rotate(90, 0, 1, 0, loc);
        Rotate(90, 0, 0, 1, loc);
        Cylinder(0.5, -0.5, 0.5, 360, params, loc);
      }
      else if (name == "PxrSphereLight") {
        // A Sphere light lies along the X-axis, and is half the default sphere size
        Rotate(90, 0, 1, 0, loc);
        Rotate(90, 0, 0, 1, loc);
        Sphere(0.5, -0.5, 0.5, 360, params, loc);
      }
      else {
        // Can't reach here yet
      }

      _light_material = new Parsed_Parameter_Vector;
      Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::String, "__materialid", loc);
      param->may_be_unused = true;
      param->add_string(handle);
      _light_material->push_back(param);

      // Check if the light has a materialid parameter
      if (entity.parameters.get_one_string("__materialid", "").empty()) {
        auto light_params = entity.parameters.get_parameter_vector();
        light_params.push_back(param);
        entity = Light_Scene_Entity(
            handle,
            name,
            Parameter_Dictionary(std::move(light_params), graphics_state.light_attributes),
            loc,
            render_from_object);
      }
      _lights[handle] = std::move(entity);
    }
    else {
      active_instance_definition->entity.lights.push_back(std::move(entity));
    }
  }
  else {
    std::string material_id = entity.parameters.get_one_string("__materialid", "");
    if (!material_id.empty()) {
      // If we're processing the light instance (only occurs with a mesh light),
      // we need to create an emission shader that's picked up by the geometry
      // at a later point.
      if (_lights.find(handle) != _lights.end()) {
        auto it = osl_shader_group.find(material_id);
        if (it != osl_shader_group.end()) {
          osl_shader_group.erase(it);
        }
        auto light_params = _lights[handle].parameters.get_parameter_vector();
        Parameter_Dictionary light_dict(light_params, graphics_state.light_attributes);
        float strength = 1.0f;

        float exposure = light_dict.get_one_float("exposure", 0.0);
        strength *= exp2(exposure);

        float intensity = light_dict.get_one_float("intensity", 1.0);
        strength *= intensity;

        param = new Parsed_Parameter(Parameter_Type::Real, "strength", loc);
        param->may_be_unused = true;
        param->add_float(strength);
        light_params.push_back(param);

        std::stringstream ss;
        ss << material_id << "_" << graphics_state.identifier;
        graphics_state.remapped_material[material_id] = ss.str();
        material_id = ss.str();
        param = light_dict.get_parameter("__materialid");
        param->strings()[0] = material_id;

        Bxdf(name, handle, light_params, loc);

        param = new Parsed_Parameter(Parameter_Type::Boolean, "areaNormalize", loc);
        param->may_be_unused = true;
        param->add_bool((light_dict.get_one_int("areaNormalize", 0) == 1));
        graphics_state.rib_attributes["Ri"].push_back(param);

        graphics_state.current_material_name = material_id;
        graphics_state.current_material_index = -1;
      }
      else {
        _lights[handle] = std::move(entity);
      }
    }
  }
}

void Ri::LightSource([[maybe_unused]] const std::string &name,
                     [[maybe_unused]] Parsed_Parameter_Vector params,
                     [[maybe_unused]] File_Loc loc)
{
  std::cout << "LightSource is unimplemented." << std::endl;
}

void Ri::MakeBrickMap([[maybe_unused]] std::vector<std::string> ptcnames,
                      [[maybe_unused]] const std::string &bkmname,
                      [[maybe_unused]] Parsed_Parameter_Vector params,
                      [[maybe_unused]] File_Loc loc)
{
  std::cout << "MakeBrickMap is unimplemented" << std::endl;
}

void Ri::MakeCubeFaceEnvironment([[maybe_unused]] const std::string &px,
                                 [[maybe_unused]] const std::string &nx,
                                 [[maybe_unused]] const std::string &py,
                                 [[maybe_unused]] const std::string &ny,
                                 [[maybe_unused]] const std::string &pz,
                                 [[maybe_unused]] const std::string &nz,
                                 [[maybe_unused]] const std::string &text,
                                 [[maybe_unused]] float fov,
                                 [[maybe_unused]] const std::string &filt,
                                 [[maybe_unused]] float swidth,
                                 [[maybe_unused]] float twidth,
                                 [[maybe_unused]] Parsed_Parameter_Vector params,
                                 [[maybe_unused]] File_Loc loc)
{
  std::cout << "MakeCubeFaceEnvironment is unimplemented" << std::endl;
}

void Ri::MakeLatLongEnvironment([[maybe_unused]] const std::string &picturename,
                                [[maybe_unused]] const std::string &texturename,
                                [[maybe_unused]] const std::string &filt,
                                [[maybe_unused]] float swidth,
                                [[maybe_unused]] float twidth,
                                [[maybe_unused]] Parsed_Parameter_Vector params,
                                [[maybe_unused]] File_Loc loc)
{
  std::cout << "MakeLatLongEnvironment is unimplemented" << std::endl;
}

void Ri::MakeShadow([[maybe_unused]] const std::string &picturename,
                    [[maybe_unused]] const std::string &texturename,
                    [[maybe_unused]] Parsed_Parameter_Vector params,
                    [[maybe_unused]] File_Loc loc)
{
  std::cout << "MakeShadow is unimplemented" << std::endl;
}

void Ri::MakeTexture([[maybe_unused]] const std::string &picturename,
                     [[maybe_unused]] const std::string &texturename,
                     [[maybe_unused]] const std::string &swrap,
                     [[maybe_unused]] const std::string &twrap,
                     [[maybe_unused]] const std::string &filt,
                     [[maybe_unused]] float swidth,
                     [[maybe_unused]] float twidth,
                     [[maybe_unused]] Parsed_Parameter_Vector params,
                     [[maybe_unused]] File_Loc loc)
{
  std::cout << "MakeTexture is unimplemented" << std::endl;
}

void Ri::Matte([[maybe_unused]] bool onoff, [[maybe_unused]] File_Loc loc)
{
  std::cout << "Matte is unimplemented" << std::endl;
}

void Ri::MotionBegin([[maybe_unused]] std::vector<float> times, [[maybe_unused]] File_Loc loc)
{
  std::cout << "MotionBegin is unimplemented" << std::endl;
}

void Ri::MotionEnd([[maybe_unused]] File_Loc loc)
{
  std::cout << "MotionEnd is unimplemented" << std::endl;
}

void Ri::NuPatch([[maybe_unused]] int nu,
                 [[maybe_unused]] int uorder,
                 [[maybe_unused]] float uknot[],
                 [[maybe_unused]] float umin,
                 [[maybe_unused]] float umax,
                 [[maybe_unused]] int nv,
                 [[maybe_unused]] int vorder,
                 [[maybe_unused]] float vknot[],
                 [[maybe_unused]] float vmin,
                 [[maybe_unused]] float vmax,
                 [[maybe_unused]] Parsed_Parameter_Vector params,
                 [[maybe_unused]] File_Loc loc)
{
  std::cout << "NuPatch is unimplemented" << std::endl;
}

void Ri::ObjectBegin(const std::string &name, File_Loc loc)
{
  VERIFY_WORLD("ObjectBegin");
  pushed_graphics_states.push_back(graphics_state);

  push_stack.emplace_back('o', loc);

  if (active_instance_definition) {
    error_exit_deferred(&loc, "ObjectBegin called inside of instance definition");
    return;
  }

  if (instance_names.find(name) != instance_names.end()) {
    std::stringstream ss;
    ss << name << ": trying to redefine an object instance";
    error_exit_deferred(&loc, ss.str());
    return;
  }
  instance_names.insert(name);
  active_instance_definition = new Active_Instance_Definition(name, loc);
}

void Ri::ObjectEnd(File_Loc loc)
{
  VERIFY_WORLD("ObjectEnd");
  if (!active_instance_definition) {
    error_exit_deferred(&loc, "ObjectEnd called outside of instance definition");
    return;
  }
  if (active_instance_definition->parent) {
    error_exit_deferred(&loc, "ObjectEnd called inside Import for instance definition");
    return;
  }

  // NOTE: Must keep the following consistent with AttributeEnd
  graphics_state = std::move(pushed_graphics_states.back());
  pushed_graphics_states.pop_back();

  if (push_stack.back().first == 'a') {
    std::stringstream ss;
    ss << "Mismatched nesting: open AttributeBegin from ";
    ss << push_stack.back().second.to_string();
    ss << " at ObjectEnd";
    error_exit_deferred(&loc, ss.str());
  }
  else {
    CHECK_EQ(push_stack.back().first, 'o');
  }
  push_stack.pop_back();

  if (--active_instance_definition->active_imports == 0) {
    add_instance_definition(std::move(active_instance_definition->entity));
    delete active_instance_definition;
  }

  active_instance_definition = nullptr;
}

void Ri::ObjectInstance(const std::string &name, File_Loc loc)
{
  VERIFY_WORLD("ObjectInstance");

  if (_light_material) {
    Light("PxrMeshLight", (*_light_material)[0]->strings()[0], *_light_material, loc);
    osl_shader_group[_shader_id] = true;
    delete _light_material;
    _light_material = nullptr;
  }

  Mapped_Parameter_Dictionary dict;
  for (auto &x : graphics_state.rib_attributes) {
    Parameter_Dictionary d(x.second);

    dict.emplace(x.first, d);
  }

  if (dict.find("identifier") == dict.end()) {
    Parsed_Parameter_Vector params;
    Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::String, "name", loc);
    param->may_be_unused = true;
    param->add_string(generate_random_alphanumeric_string());
    params.push_back(param);

    Parameter_Dictionary d(params);
    dict.emplace("identifier", d);
  }

  if (active_instance_definition) {
    error_exit_deferred(&loc, "ObjectInstance can't be called inside instance definition");
    return;
  }

  ProjectionTransform world_from_render = projection_inverse(render_from_world);

#if 0
      if ( CTM_Is_Animated() )
      {
         Animated_Transform animated_render_from_instance(
             Render_From_Object( 0 ) * world_from_render,
             graphics_state.transform_start_time,
             Render_From_Object( 1 ) * world_from_render,
             graphics_state.transform_end_time );

         // For very small changes, animated_render_from_instance may have both
         // xforms equal even if CTMIsAnimated() has returned true. Fall
         // through to create a regular non-animated instance in that case.
         if ( animated_render_from_instance.is_animated() )
         {
            instance_uses.push_back( Instance_Scene_Entity(
                name, loc, graphics_state.current_material_name,
                animated_render_from_instance ) );
            return;
         }
      }

      const class Transform* render_from_instance = scene->transform_cache.lookup(
          Render_From_Object( 0 ) * world_from_render );
#endif

  const ProjectionTransform *render_from_instance = transform_cache.lookup(Render_From_Object(0) *
                                                                           world_from_render);

  add_instance_use(Instance_Scene_Entity(name,
                                         loc,
                                         graphics_state.current_material_name,
                                         std::move(dict) /* RenderMan*/,
                                         render_from_instance));
}

void Ri::Option(const std::string &name,
                Parsed_Parameter_Vector params,
                [[maybe_unused]] File_Loc loc)
{
  auto items = _rib_state.options.find(name);
  if (items == _rib_state.options.end()) {
    _rib_state.options[name] = {};
    items = _rib_state.options.find(name);
  }

  for (auto *p : params) {
    // FIX_ME: If the name is already in the map, this will
    // result in a memory leak
    items->second[p->name] = p;
  }
}

void Ri::Opacity([[maybe_unused]] Parsed_Parameter_Vector params, [[maybe_unused]] File_Loc loc)
{
  std::cout << "Opacity is unimplemented" << std::endl;
}

void Ri::Orientation([[maybe_unused]] const std::string &orientation,
                     [[maybe_unused]] File_Loc loc)
{
  std::cout << "Orientation is unimplemented" << std::endl;
}

void Ri::Paraboloid([[maybe_unused]] float rmax,
                    [[maybe_unused]] float zmin,
                    [[maybe_unused]] float zmax,
                    [[maybe_unused]] float thetamax,
                    [[maybe_unused]] Parsed_Parameter_Vector params,
                    [[maybe_unused]] File_Loc loc)
{
  std::cout << "Paraboloid is unimplemented" << std::endl;
}

void Ri::Patch([[maybe_unused]] const std::string &type,
               [[maybe_unused]] Parsed_Parameter_Vector params,
               [[maybe_unused]] File_Loc loc)
{
  std::cout << "Patch is unimplemented" << std::endl;
}

void Ri::PatchMesh([[maybe_unused]] const std::string &type,
                   [[maybe_unused]] int nu,
                   [[maybe_unused]] const std::string &uwrap,
                   [[maybe_unused]] int nv,
                   [[maybe_unused]] const std::string &vwrap,
                   [[maybe_unused]] Parsed_Parameter_Vector params,
                   [[maybe_unused]] File_Loc loc)
{
  std::cout << "PatchMesh is unimplemented" << std::endl;
}

void Ri::Pattern(const std::string &name,
                 const std::string &handle,
                 Parsed_Parameter_Vector params,
                 File_Loc loc)
{
  VERIFY_WORLD("Pattern");

  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::String, "shader_type", loc);
  param->add_string("pattern");
  param->add_string(name);
  param->add_string(handle);
  params.push_back(param);

  Parameter_Dictionary dict(std::move(params));
  osl_parameters.push_back(dict);
}

void Ri::Perspective([[maybe_unused]] float fov, [[maybe_unused]] File_Loc loc)
{
  std::cout << "Perspective is unimplemented" << std::endl;
}

void Ri::PixelFilter([[maybe_unused]] const std::string &name,
                     [[maybe_unused]] Parsed_Parameter_Vector params,
                     [[maybe_unused]] File_Loc loc)
{
  std::cout << "PixelFilter is unimplemented." << std::endl;
}

void Ri::PixelSampleImager([[maybe_unused]] const std::string &name,
                           [[maybe_unused]] Parsed_Parameter_Vector params,
                           [[maybe_unused]] File_Loc loc)
{
  std::cout << "PixelSampleImager is unimplemented" << std::endl;
}

void Ri::PixelSamples([[maybe_unused]] float xsamples,
                      [[maybe_unused]] float ysamples,
                      [[maybe_unused]] File_Loc loc)
{
  std::cout << "PixelSamples is unimplemented" << std::endl;
}

void Ri::PixelVariance([[maybe_unused]] float variance, [[maybe_unused]] File_Loc loc)
{
  std::cout << "PixelVariance is unimplemented" << std::endl;
}

void Ri::Points([[maybe_unused]] int npoints,
                [[maybe_unused]] Parsed_Parameter_Vector params,
                [[maybe_unused]] File_Loc loc)
{
  std::cout << "Points is unimplemented" << std::endl;
}

void Ri::PointsGeneralPolygons(std::vector<int> n_loops,
                               std::vector<int> n_vertices,
                               std::vector<int> vertices,
                               Parsed_Parameter_Vector params,
                               File_Loc loc)
{
  bool single_only = true;

  // Only single loops for now
  for (int i = 0; single_only && i < n_loops.size(); i++) {
    single_only = ((int)n_loops[i] == 1);
  }

  if (single_only) {
    PointsPolygons(n_vertices, vertices, params, loc);
  }
}

void Ri::PointsPolygons(std::vector<int> n_vertices,
                        std::vector<int> vertices,
                        Parsed_Parameter_Vector params,
                        File_Loc loc)
{
  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::Integer, "vertices", loc);
  for (int i = 0; i < vertices.size(); ++i) {
    param->add_int(vertices[i]);
  }
  params.push_back(param);

  int nfaces = n_vertices.size();
  param = new Parsed_Parameter(Parameter_Type::Integer, "nvertices", loc);
  for (int i = 0; i < n_vertices.size(); ++i) {
    param->add_int(n_vertices[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nfaces", loc);
  param->add_int(nfaces);
  params.push_back(param);

  // Fix attributes whose storage class couldn't be determined at parsing
  // time
  for (auto &p : params) {
    if ((p->type == Parameter_Type::Point3 || p->type == Parameter_Type::Normal ||
         p->type == Parameter_Type::Color) &&
        p->storage == Container_Type::Constant && p->floats().size() > 3)
    {
      int num_vals = p->floats().size() / 3;
      if (num_vals == n_vertices.size()) {
        p->storage = Container_Type::Uniform;
      }
      else if (num_vals == vertices.size()) {
        p->storage = Container_Type::Vertex;
      }
    }
  }

  bool has_varying_normals = false;
  bool has_uvs = false;
  for (auto &p : params) {
    if (p->type == Parameter_Type::Normal) {
      if (p->storage == Container_Type::FaceVarying || p->storage == Container_Type::Varying ||
          p->storage == Container_Type::Vertex)
      {
        has_varying_normals = true;
      }
    }
    else if (p->type == Parameter_Type::Point2 && (p->name == "st" || p->name == "uv")) {
      has_uvs = true;
    }
  }

  if (has_varying_normals) {
    param = new Parsed_Parameter(Parameter_Type::Boolean, "smooth", loc);
    param->add_bool(true);
    params.push_back(param);
  }

  if (!has_uvs) {
    param = new Parsed_Parameter(Parameter_Type::Point2, "uv", loc);
    param->storage = Container_Type::FaceVarying;
    param->elem_per_item = 2;
    std::vector<float> tc = {0, 0, 1, 0, 1, 1, 0, 1};
    for (int i = 0; i < 2 * vertices.size(); i++) {
      param->add_float(0);
    }
    vector<float> &tex_coords = param->floats();
    int index = 0;
    for (auto i = 0; i < n_vertices.size(); i++) {
      for (auto j = 0; j < n_vertices[i]; j++) {
        tex_coords[index++] = tc[2 * j];
        tex_coords[index++] = tc[2 * j + 1];
      }
    }
    params.push_back(param);
  }

  Shape("mesh", params, loc);
}

void Ri::Polygon([[maybe_unused]] int nvertices,
                 [[maybe_unused]] Parsed_Parameter_Vector params,
                 [[maybe_unused]] File_Loc loc)
{
  std::cout << "Polygon is unimplemented" << std::endl;
}

void Ri::ProcDelayedReadArchive([[maybe_unused]] Parsed_Parameter_Vector params,
                                [[maybe_unused]] File_Loc loc)
{
  std::cout << "ProcDelayedReadArchive is unimplemented" << std::endl;
}

void Ri::ProcDynamicLoad([[maybe_unused]] Parsed_Parameter_Vector params,
                         [[maybe_unused]] File_Loc loc)
{
  std::cout << "ProcDynamicLoad is unimplemented" << std::endl;
}

void Ri::Procedural([[maybe_unused]] const std::string &proc_name,
                    [[maybe_unused]] Parsed_Parameter_Vector params,
                    [[maybe_unused]] File_Loc loc)
{
  std::cout << "Procedural is unimplemented" << std::endl;
}

void Ri::Procedural2([[maybe_unused]] const std::string &proc_name,
                     [[maybe_unused]] Parsed_Parameter_Vector params,
                     [[maybe_unused]] File_Loc loc)
{
  std::cout << "Procedural2 is unimplemented" << std::endl;
}

void Ri::ProcFree([[maybe_unused]] Parsed_Parameter_Vector params, [[maybe_unused]] File_Loc loc)
{
  std::cout << "ProcFree is unimplemented" << std::endl;
}

void Ri::ProcRunProgram([[maybe_unused]] Parsed_Parameter_Vector params,
                        [[maybe_unused]] File_Loc loc)
{
  std::cout << "ProcRunProgram is unimplemented" << std::endl;
}

void Ri::Projection(const std::string &name,
                    Parsed_Parameter_Vector params,
                    [[maybe_unused]] File_Loc loc)
{
  auto items = _rib_state.options.find("Ri");
  if (items == _rib_state.options.end()) {
    _rib_state.options["Ri"] = {};
    items = _rib_state.options.find(name);
  }

  for (auto *p : params) {
    // FIX_ME: If the name already is in the map, this will
    // result in a memory leak
    items->second[p->name] = p;
  }
}

void Ri::Quantize([[maybe_unused]] const std::string &type,
                  [[maybe_unused]] int one,
                  [[maybe_unused]] int min,
                  [[maybe_unused]] int max,
                  [[maybe_unused]] float ditheramplitude,
                  [[maybe_unused]] File_Loc loc)
{
  std::cout << "Quantize is unimplemented" << std::endl;
}

void Ri::ReadArchive([[maybe_unused]] const std::string &name,
                     [[maybe_unused]] Parsed_Parameter_Vector params,
                     [[maybe_unused]] File_Loc loc)
{
  std::cout << "ReadArchive is unimplemented" << std::endl;
}

void Ri::RelativeDetail([[maybe_unused]] float relativedetail, [[maybe_unused]] File_Loc loc)
{
  std::cout << "RelativeDetail is unimplemented" << std::endl;
}

void Ri::Resource([[maybe_unused]] const std::string &handle,
                  [[maybe_unused]] const std::string &type,
                  [[maybe_unused]] Parsed_Parameter_Vector params,
                  [[maybe_unused]] File_Loc loc)
{
  std::cout << "Resource is unimplemented" << std::endl;
}

void Ri::ResourceBegin([[maybe_unused]] File_Loc loc)
{
  std::cout << "ResourceBegin is unimplemented" << std::endl;
}

void Ri::ResourceEnd([[maybe_unused]] File_Loc loc)
{
  std::cout << "ResourceEnd is unimplemented" << std::endl;
}

void Ri::ReverseOrientation([[maybe_unused]] File_Loc loc)
{
  VERIFY_WORLD("ReverseOrientation");
  graphics_state.reverse_orientation = !graphics_state.reverse_orientation;
}

void Ri::Rotate(float angle, float ax, float ay, float az, [[maybe_unused]] File_Loc loc)
{
  graphics_state.for_active_transforms(
      [=](auto t) { return t * transform_rotate(DEG2RADF(angle), make_float3(ax, ay, az)); });
}

void Ri::Scale(float sx, float sy, float sz, [[maybe_unused]] File_Loc loc)
{
  graphics_state.for_active_transforms(
      [=](auto t) { return t * transform_scale(make_float3(sx, sy, sz)); });
}

void Ri::ScopedCoordinateSystem([[maybe_unused]] const std::string &name,
                                [[maybe_unused]] File_Loc loc)
{
  std::cout << "ScopedCoordinateSystem is unimplemented" << std::endl;
}

void Ri::ScreenWindow([[maybe_unused]] float left,
                      [[maybe_unused]] float right,
                      [[maybe_unused]] float bottom,
                      [[maybe_unused]] float top,
                      [[maybe_unused]] File_Loc loc)
{
  std::cout << "ScreenWindow is unimplemented" << std::endl;
}

void Ri::ShadingInterpolation([[maybe_unused]] const std::string &type,
                              [[maybe_unused]] File_Loc loc)
{
  std::cout << "ShadingInterpolation is unimplemented" << std::endl;
}

void Ri::ShadingRate([[maybe_unused]] float size, [[maybe_unused]] File_Loc loc)
{
  std::cout << "ShadingRate is unimplemented" << std::endl;
}

void Ri::Shutter([[maybe_unused]] float opentime,
                 [[maybe_unused]] float closetime,
                 [[maybe_unused]] File_Loc loc)
{
  std::cout << "Shutter is unimplemented" << std::endl;
}

void Ri::Sides(int nsides, File_Loc loc)
{
  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::Integer, "sides", loc);
  param->add_float(nsides);
  graphics_state.shape_attributes.push_back(param);
}

void Ri::Skew([[maybe_unused]] float angle,
              [[maybe_unused]] float dx1,
              [[maybe_unused]] float dy1,
              [[maybe_unused]] float dz1,
              [[maybe_unused]] float dx2,
              [[maybe_unused]] float dy2,
              [[maybe_unused]] float dz2,
              [[maybe_unused]] File_Loc loc)
{
  std::cout << "Skew is unimplemented" << std::endl;
}

void Ri::SolidBegin([[maybe_unused]] const std::string &type, [[maybe_unused]] File_Loc loc)
{
  std::cout << "SolidBegin is unimplemented" << std::endl;
}

void Ri::SolidEnd([[maybe_unused]] File_Loc loc)
{
  std::cout << "SolidEnd is unimplemented" << std::endl;
}

void Ri::Sphere(float radius,
                float zmin,
                float zmax,
                float thetamax,
                Parsed_Parameter_Vector params,
                File_Loc loc)
{
  VERIFY_WORLD("Shape");

  // Ensure thetamax in [0, 360]
  thetamax = (thetamax > 360.0f ? 360.0f : (thetamax < 0.f ? 0.f : thetamax));
  thetamax = thetamax * M_PI_F / 180.0f;
  zmin = min(max(-radius, zmin), radius);
  zmax = max(min(radius, zmax), -radius);
  float end_phi = M_PI_F;

  if (zmin > -radius) {
    end_phi = (M_PI_F / 2.0f) - std::asinf(zmin / radius);
  }

  float start_phi = 0;
  if (zmax < radius) {
    start_phi = (M_PI_F / 2.0f) - std::asinf(zmax / radius);
  }

  float delta_phi = end_phi - start_phi;

  std::vector<float> pts;
  std::vector<float2> uvs, top_uvs, body_uvs, bot_uvs;
  std::vector<float> norms;
  std::vector<int> top_ring;
  std::vector<int> bot_ring;
  std::vector<int> body;
  int point_index = 0;
  int n_slices = 25;
  int n_stacks = 24;

  // Allocate memory
  top_uvs.reserve(n_slices + 1);
  bot_uvs.reserve(n_slices + 1);
  body_uvs.reserve(n_slices * n_stacks);
  top_ring.reserve(n_slices);
  bot_ring.reserve(n_slices);
  pts.reserve(n_slices * (n_stacks + 1) * 3);
  norms.reserve(n_slices * (n_stacks + 1) * 3);

  // Create top vertices
  float phi = start_phi;
  for (int j = 0; j < n_slices; ++j) {
    float theta = thetamax * float(j) / float(n_slices);
    float3 dir = ri_spherical_to_direction(-theta, phi);
    norms.push_back(dir.x);
    norms.push_back(dir.y);
    norms.push_back(dir.z);
    pts.push_back(radius * dir.x);
    pts.push_back(radius * dir.y);
    pts.push_back(radius * dir.z);
    top_uvs.push_back(make_float2(theta / thetamax, (end_phi - phi) / delta_phi));
    top_ring.push_back(point_index++);
  }
  top_uvs.push_back(make_float2(1.f, (end_phi - phi) / delta_phi));

  // create body vertices
  for (int i = 0; i < n_stacks - 1; ++i) {
    float phi = start_phi + delta_phi * float(i + 1) / float(n_stacks);
    for (int j = 0; j < n_slices; ++j) {
      float theta = thetamax * float(j) / float(n_slices);
      float3 dir = ri_spherical_to_direction(-theta, phi);
      norms.push_back(dir.x);
      norms.push_back(dir.y);
      norms.push_back(dir.z);
      pts.push_back(radius * dir.x);
      pts.push_back(radius * dir.y);
      pts.push_back(radius * dir.z);
      body_uvs.push_back(make_float2(theta / thetamax, (end_phi - phi) / delta_phi));
      body.push_back(point_index++);
    }
    body_uvs.push_back(make_float2(1.f, (end_phi - phi) / delta_phi));
  }

  // create bottom vertices
  phi = end_phi;
  for (int j = 0; j < n_slices; ++j) {
    float theta = thetamax * float(j) / float(n_slices);
    float3 dir = ri_spherical_to_direction(-theta, phi);
    norms.push_back(dir.x);
    norms.push_back(dir.y);
    norms.push_back(dir.z);
    pts.push_back(radius * dir.x);
    pts.push_back(radius * dir.y);
    pts.push_back(radius * dir.z);
    bot_uvs.push_back(make_float2(theta / thetamax, (end_phi - phi) / delta_phi));
    bot_ring.push_back(point_index++);
  }
  bot_uvs.push_back(make_float2(1.f, (end_phi - phi) / delta_phi));

  int poly_count = 0;
  std::vector<int> polys;
  std::vector<int> poly_counts;

  uvs.reserve(n_slices * n_stacks * 4);
  polys.reserve(n_slices * n_stacks * 4);
  poly_counts.reserve((n_slices - 1) * n_stacks);

  // Create end cap tris/quads
  for (int i = 0; i < n_slices; ++i) {
    int i0 = i;
    int i1 = (i + 1);
    int i2 = (i + 1);
    int i3 = i;
    polys.push_back(body[i3]);
    polys.push_back(body[i2 % n_slices]);
    polys.push_back(top_ring[i1 % top_ring.size()]);
    polys.push_back(top_ring[i0]);
    uvs.push_back(body_uvs[i3]);
    uvs.push_back(body_uvs[i2]);
    uvs.push_back(top_uvs[i1]);
    uvs.push_back(top_uvs[i0]);
    poly_counts.push_back(4);
    poly_count = poly_count + 1;
#if 0
    std::cout << std::setw(3) << poly_count << ", " << std::setw(3) << poly_counts.back() << ", ";
    std::cout << std::setw(3) << i0 << ", " << std::setw(3) << i1 << ", " << std::setw(3) << i2
              << ", " << std::setw(3) << i3 << ", ";
    std::cout << std::setw(3) << top_ring[i0] << ", " << std::setw(3) << top_ring[i1] << ", "
              << std::setw(3) << body[i2] << ", " << std::setw(3) << body[i3] << std::endl;
#endif
  }

  for (int i = 0; i < n_slices; ++i) {
    int i0 = i + n_slices * (n_stacks - 2);
    int i1 = (i + 1) % n_slices + n_slices * (n_stacks - 2);
    int i2 = (i + 1) % bot_ring.size();
    int i3 = i % bot_ring.size();
    polys.push_back(bot_ring[i3]);
    polys.push_back(bot_ring[i2]);
    polys.push_back(body[i1]);
    polys.push_back(body[i0]);
    uvs.push_back(bot_uvs[i3]);
    uvs.push_back(bot_uvs[i2]);
    uvs.push_back(body_uvs[i1]);
    uvs.push_back(body_uvs[i0]);
    poly_count = poly_count + 1;
#if 0
    std::cout << std::endl;
    std::cout << std::setw(3) << poly_count << ", " << std::setw(3) << poly_counts.back() << ", ";
    std::cout << std::setw(3) << i0 << ", " << std::setw(3) << i1 << ", " << std::setw(3) << i2
              << ", " << std::setw(3) << i3 << ", ";
    std::cout << std::setw(3) << body[i0] << ", " << std::setw(3) << bot_ring[i3] << ", "
              << std::setw(3) << bot_ring[i2] << ", " << std::setw(3) << body[i1] << std::endl;
#endif
  }

  // add quads per stack / slice
  for (int j = 0; j < n_stacks - 2; j++) {
    int j0 = j * n_slices;
    int j1 = (j + 1) * n_slices;
    for (int i = 0; i < n_slices; i++) {
      int i0 = j0 + i;
      int i1 = j0 + (i + 1) % n_slices;
      int i2 = j1 + (i + 1) % n_slices;
      int i3 = j1 + i;
      polys.push_back(body[i3]);
      polys.push_back(body[i2]);
      polys.push_back(body[i1]);
      polys.push_back(body[i0]);
      poly_counts.push_back(4);
      poly_count = poly_count + 1;
#if 0
      std::cout << std::setw(3) << poly_count << ", " << std::setw(3) << poly_counts.back()
                << ", ";
      std::cout << std::setw(3) << i0 << ", " << std::setw(3) << i1 << ", " << std::setw(3) << i2
                << ", " << std::setw(3) << i3 << ", ";
      std::cout << std::setw(3) << body[i3] << ", " << std::setw(3) << body[i2] << ", "
                << std::setw(3) << body[i1] << ", " << std::setw(3) << body[i0] << std::endl;
#endif
    }
  }

  for (int j = 0; j < n_stacks - 2; j++) {
    int j0 = j * (n_slices + 1);
    int j1 = (j + 1) * (n_slices + 1);
    for (int i = 0; i < n_slices; i++) {
      int i0 = j0 + i;
      int i1 = j0 + (i + 1) % (n_slices + 1);
      int i2 = j1 + (i + 1) % (n_slices + 1);
      int i3 = j1 + i;
      uvs.push_back(body_uvs[i3]);
      uvs.push_back(body_uvs[i2]);
      uvs.push_back(body_uvs[i1]);
      uvs.push_back(body_uvs[i0]);
    }
  }

  Parsed_Parameter *param = new Parsed_Parameter(
      Parameter_Type::Integer, "vertices", loc, polys.size());
  for (int i = 0; i < polys.size(); ++i) {
    param->add_int(polys[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nvertices", loc, poly_counts.size());
  for (int i = 0; i < poly_counts.size(); ++i) {
    param->add_int(poly_counts[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Point3, "P", loc, pts.size());
  param->storage = Container_Type::Vertex;
  param->elem_per_item = 3;
  for (int i = 0; i < pts.size(); ++i) {
    param->add_float(pts[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Point2, "uv", loc, 2 * uvs.size());
  param->storage = Container_Type::FaceVarying;
  param->elem_per_item = 2;
  for (int i = 0; i < uvs.size(); ++i) {
    param->add_float(uvs[i].x);
    param->add_float(uvs[i].y);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nfaces", loc);
  param->add_int(poly_count);
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Boolean, "smooth", loc);
  param->add_bool(true);
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Normal, "N", loc, norms.size());
  param->storage = Container_Type::Varying;
  param->elem_per_item = 3;
  for (int i = 0; i < norms.size(); ++i) {
    param->add_float(norms[i]);
  }
  params.push_back(param);

  // rotate to match RenderMan orientation
  Rotate(90, 1, 0, 0, loc);
  Shape("mesh", params, loc);
}

void Ri::SubdivisionMesh(const std::string &scheme,
                         int nfaces,
                         std::vector<int> n_vertices,
                         std::vector<int> vertices,
                         std::vector<std::string> tags,
                         std::vector<int> nargs,
                         std::vector<int> intargs,
                         std::vector<float> floatargs,
                         Parsed_Parameter_Vector params,
                         File_Loc loc)
{
  Parsed_Parameter *param = new Parsed_Parameter(Parameter_Type::String, "scheme", loc);
  param->add_string(scheme);
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nfaces", loc);
  param->add_int(nfaces);
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "vertices", loc);
  for (int i = 0; i < vertices.size(); ++i) {
    param->add_int(vertices[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nvertices", loc);
  for (int i = 0; i < n_vertices.size(); ++i) {
    param->add_int(n_vertices[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::String, "tags", loc);
  for (int i = 0; i < tags.size(); ++i) {
    param->add_string(tags[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "nargs", loc);
  for (int i = 0; i < nargs.size(); ++i) {
    param->add_int(nargs[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Integer, "intargs", loc);
  for (int i = 0; i < intargs.size(); ++i) {
    param->add_int(intargs[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(Parameter_Type::Real, "floatargs", loc);
  for (int i = 0; i < floatargs.size(); ++i) {
    param->add_float(floatargs[i]);
  }
  params.push_back(param);

  Shape("subdivision_mesh", params, loc);
}

void Ri::Surface([[maybe_unused]] const std::string &name,
                 [[maybe_unused]] Parsed_Parameter_Vector params,
                 [[maybe_unused]] File_Loc loc)
{
  std::cout << "Surface is unimplemented" << std::endl;
}

void Ri::System([[maybe_unused]] const std::string &cmd, [[maybe_unused]] File_Loc loc)
{
  std::cout << "System is unimplemented" << std::endl;
}

void Ri::Texture([[maybe_unused]] const std::string &name,
                 [[maybe_unused]] const std::string &type,
                 [[maybe_unused]] const std::string &texname,
                 [[maybe_unused]] Parsed_Parameter_Vector params,
                 [[maybe_unused]] File_Loc loc)
{
  std::cout << "Texture is unimplemented." << std::endl;
}

void Ri::TextureCoordinates([[maybe_unused]] float s1,
                            [[maybe_unused]] float t1,
                            [[maybe_unused]] float s2,
                            [[maybe_unused]] float t2,
                            [[maybe_unused]] float s3,
                            [[maybe_unused]] float t3,
                            [[maybe_unused]] float s4,
                            [[maybe_unused]] float t4,
                            [[maybe_unused]] File_Loc loc)
{
  std::cout << "TextureCoordinates is unimplemented" << std::endl;
}

void Ri::Torus([[maybe_unused]] float majorrad,
               [[maybe_unused]] float minorrad,
               [[maybe_unused]] float phimin,
               [[maybe_unused]] float phimax,
               [[maybe_unused]] float thetamax,
               [[maybe_unused]] Parsed_Parameter_Vector params,
               [[maybe_unused]] File_Loc loc)
{
  std::cout << "Torus is unimplemented" << std::endl;
}

void Ri::transform(float const *transform, [[maybe_unused]] File_Loc loc)
{
  graphics_state.for_active_transforms([=]([[maybe_unused]] auto t) {
    // Stomp the current transform
    ProjectionTransform projection = *(ProjectionTransform *)&transform[0];
    return projection_transpose(projection);
  });
}

void Ri::TransformBegin(File_Loc loc)
{
  pushed_graphics_states.push_back(graphics_state);
  push_stack.emplace_back('t', loc);
}

void Ri::TransformEnd(File_Loc loc)
{
  // Issue error on unmatched _AttributeEnd_
  if (pushed_graphics_states.empty()) {
    error(&loc, "Unmatched TransformEnd encountered. Ignoring it.");
    return;
  }

  // NOTE: must keep the following consistent with code in ObjectEnd
  // We're treating a TransformBegin/End just like it's Attribute equivilent
  auto old_graphics_state = std::move(pushed_graphics_states.back());
  pushed_graphics_states.pop_back();
  // Keep any attributes that changed and revert just the transform
  graphics_state.ctm = old_graphics_state.ctm;

  if (push_stack.back().first == 'a') {
    std::stringstream ss;
    ss << "Mismatched nesting: open AttributeBegin from ";
    ss << push_stack.back().second.to_string();
    ss << " at TransformEnd";
    error_exit_deferred(&loc, ss.str());
  }
  else if (push_stack.back().first == 'o') {
    std::stringstream ss;
    ss << "Mismatched nesting: open ObjectBegin from ";
    ss << push_stack.back().second.to_string();
    ss << " at TransformEnd";
    error_exit_deferred(&loc, ss.str());
  }
  else {
    CHECK_EQ(push_stack.back().first, 't');
  }
  push_stack.pop_back();
}

void Ri::Translate(float dx, float dy, float dz, [[maybe_unused]] File_Loc loc)
{
  graphics_state.for_active_transforms(
      [=](auto t) { return t * transform_translate(make_float3(dx, dy, dz)); });
}

void Ri::TrimCurve([[maybe_unused]] int nloops,
                   [[maybe_unused]] int ncurves[],
                   [[maybe_unused]] int order[],
                   [[maybe_unused]] float knot[],
                   [[maybe_unused]] float min[],
                   [[maybe_unused]] float max[],
                   [[maybe_unused]] int n[],
                   [[maybe_unused]] float u[],
                   [[maybe_unused]] float v[],
                   [[maybe_unused]] float w[],
                   [[maybe_unused]] File_Loc loc)
{
  std::cout << "TrimCurve is unimplemented" << std::endl;
}

void Ri::WorldBegin(File_Loc loc)
{
  VERIFY_OPTIONS("WorldBegin");
  // Reset graphics state for _WorldBegin_
  current_block = Block_State::World_Block;
  for (int i = 0; i < Max_Transforms; ++i) {
    graphics_state.ctm[i] = projection_identity();
  }
  named_coordinate_systems["world"] = graphics_state.ctm;

  Parsed_Parameter_Vector params;
  Parsed_Parameter *param;

  auto options = _rib_state.options["Ri"];
  auto *opt_param = options["PixelFilterName"];
  if (opt_param != nullptr) {
    std::string name = opt_param->strings()[0];
    float xradius = 2, yradius = 2;

    auto *opt_width = options["PixelFilterWidth"];
    if (opt_width != nullptr) {
      xradius = opt_width->floats()[0];
      yradius = opt_width->floats()[1];
    }

    param = new Parsed_Parameter(Parameter_Type::Real, "xradius", loc);
    param->add_float(xradius);
    params.push_back(param);

    param = new Parsed_Parameter(Parameter_Type::Real, "yradius", loc);
    param->add_float(yradius);
    params.push_back(param);

    Parameter_Dictionary dict(std::move(params));
    filter = Scene_Entity(name, std::move(dict), loc);
  }

  opt_param = options["PixelVariance"];
  if (opt_param != nullptr) {
    float pixel_var = opt_param->floats()[0];

    param = new Parsed_Parameter(Parameter_Type::Real, "pixel_variance", loc);
    param->add_float(pixel_var);
    sampler.parameters.push_back(param);
  }

  options = _rib_state.options["searchpath"];
  opt_param = options["shader"];
  std::string paths = opt_param->strings()[0] + ":";

  options = _rib_state.options["hider"];
  opt_param = options["minsamples"];
  if (opt_param != nullptr) {
    int min_samples = opt_param->ints()[0];
    opt_param = options["maxsamples"];
    int max_samples = opt_param->ints()[0];

    param = new Parsed_Parameter(Parameter_Type::Integer, "min_samples", loc);
    param->add_int(min_samples);
    sampler.parameters.push_back(param);

    param = new Parsed_Parameter(Parameter_Type::Integer, "max_samples", loc);
    param->add_int(max_samples);
    sampler.parameters.push_back(param);
  }
}

void Ri::WorldEnd([[maybe_unused]] File_Loc loc)
{
  std::cout << "WorldEnd is unimplemented" << std::endl;
}

void Ri::end_of_files()
{
  // Do end of parsing operations
}

void Ri::add_default_search_paths(std::string filepath)
{
  Parsed_Parameter_Vector params;
  Parsed_Parameter *param;

  param = new Parsed_Parameter(Parameter_Type::String, "shader_default", File_Loc());
  param->payload = vector<std::string>();
  param->add_string(filepath);
  params.push_back(param);

  Option("searchpath", params, File_Loc());
}

void Ri::Shape(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  VERIFY_WORLD("Shape");

  Parameter_Dictionary dict(std::move(params), graphics_state.shape_attributes);

  int area_light_index = -1;
  if (!graphics_state.area_light_name.empty()) {
    area_light_index = add_area_light(Scene_Entity(graphics_state.area_light_name,
                                                   graphics_state.area_light_params,
                                                   graphics_state.area_light_loc));
  }
#if 0
  if (CTM_Is_Animated()) {
    Animated_Transform renderFromShape = Render_From_Object();
    const class Transform *identity = scene->transform_cache.lookup(papillon::Transform());

    Animated_Shape_Scene_Entity entity({name,
                                        std::move(dict),
                                        loc,
                                        renderFromShape,
                                        identity,
                                        graphics_state.reverse_orientation,
                                        graphics_state.current_material_index,
                                        graphics_state.current_material_name,
                                        area_light_index,
                                        graphics_state.current_inside_medium,
                                        graphics_state.current_outside_medium});

    if (active_instance_definition)
      active_instance_definition->entity.animated_shapes.push_back(std::move(entity));
    else
      scene->add_animated_shape(std::move(entity));
  }
  else {
#endif
  ProjectionTransform const *render_from_object = transform_cache.lookup(Render_From_Object(0));
  ProjectionTransform const *object_from_render = transform_cache.lookup(
      projection_inverse(*render_from_object));

  std::string material_id = dict.get_one_string("__materialid", "");

  Shape_Scene_Entity entity({name,
                             std::move(dict),
                             loc,
                             render_from_object,
                             object_from_render,
                             area_light_index,
                             graphics_state});
  if (active_instance_definition) {
    active_instance_definition->entity.shapes.push_back(std::move(entity));
  }
  else {
    shapes.push_back(std::move(entity));
  }
}

#ifdef WITH_CYCLES_DISTRIBUTED
void Ri::init_request_thread()
{
  auto options = Option("distributed");
  auto *opt_param = options["class"];
  if (opt_param != nullptr) {
    Distributed *distributed = static_cast<Distributed *>(opt_param->pointers()[0]);
    auto create = [this](Distributed *distributed) {
      bool done = false;
      do {
        auto status = distributed->inter_comm_world.probe(1, mpl::tag_t::any());
        switch (static_cast<int>(status.tag())) {
          case Distributed::request_checksum: {
            break;
          }
          case Distributed::request_file: {
            break;
          }
          default: {
            std::cerr << "Unknown message of type: " << static_cast<int>(status.tag())
                      << std::endl;
          }
        }
      } while (!done);
      return true;
    };
    _remote_request_job = run_async(create, distributed);
  }
}
#endif

CCL_NAMESPACE_END
