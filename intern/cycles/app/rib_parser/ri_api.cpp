#include <cmath>
#include <fcntl.h>
#include <sstream>
#include <sys/stat.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <sys/sysctl.h>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>

#include <double-conversion/double-conversion.h>
namespace bfs = boost::filesystem;

#include "kernel/types.h"
#include "scene/camera.h"
#include "scene/light.h"
#include "scene/scene.h"
#include "scene/shader_graph.h"
#include "scene/shader_nodes.h"
#include "session/session.h"
#include "util/log.h"
#include "util/math.h"
#include "util/projection.h"
#include "util/transform.h"

#include "exporters/geometry.h"
#include "exporters/lights.h"
#include "exporters/materials/materials.h"

#include "app/cycles_xml.h"
#include "app/rib_parser/parsed_parameter.h"
#include "error.h"
#include "ri_api.h"

CCL_NAMESPACE_BEGIN

static bfs::path search_directory;

float3 spherical_to_direction(float theta, float phi)
{
  float sin_phi = sinf(phi);
  return make_float3(sin_phi * cosf(theta), cosf(phi), sin_phi * sinf(theta));
}

std::vector<std::string> split_string(std::string_view str, char ch)
{
  std::vector<std::string> strings;

  if (str.empty())
    return strings;

  std::string_view::iterator begin = str.begin();
  while (true) {
    std::string_view::iterator end = begin;
    while (end != str.end() && *end != ch)
      ++end;

    strings.push_back(std::string(begin, end));
    if (end == str.end())
      break;

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
  RIBCyclesMaterials materials(session->scene, osl_shader_group);
  materials.export_materials();
  for (auto &inst : instance_uses) {
    auto inst_def = instance_definitions[inst.first];
    if (inst_def->lights.size() > 0)
      export_lights(session->scene, inst.second, inst_def);
    else if (inst_def->shapes.size() > 0) {
      RIBCyclesMesh mesh(session->scene, inst.second, inst_def);
      mesh.export_geometry();
      scene_bounds.grow(mesh.bounds());
    }
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

void Ri::export_options(Scene_Entity &filter,
                        Scene_Entity &film,
                        Camera_Scene_Entity &camera,
                        Scene_Entity &sampler)
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
  if (exposureTime <= 0)
    error_exit(&camera.loc,
               "The specified camera shutter times imply that the shutter "
               "does not open.  A black image will result.");

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
  if (near > 0)
    cam->set_nearclip(near * metersPerUnit);
  float far = camera.parameters.get_one_float("farClip", -1.f);
  if (far >= 0)
    cam->set_farclip(far * metersPerUnit);

  // Set FOV
  float fov = camera.parameters.get_one_float("fov", 45.f);
  cam->set_fov(radians(fov));

  // Set screen window
  auto screen_window = camera.parameters.get_float_array("ScreenWindow");
  cam->set_viewplane_left(screen_window[0]);
  cam->set_viewplane_right(screen_window[1]);
  cam->set_viewplane_bottom(screen_window[2]);
  cam->set_viewplane_top(screen_window[3]);

  cam->need_flags_update = true;
  cam->update(session->scene);
}

void Ri::add_light(Light_Scene_Entity light)
{
  // Medium light_medium = get_medium(light.medium, &light.loc);
  std::lock_guard<std::mutex> lock(light_mutex);

#if 0
  auto create = [this, light, light_medium]() {
    return Light::create(light.name,
                         light.parameters,
                         light.render_from_object.start_transform,
                         get_camera().get_camera_transform(),
                         light_medium,
                         &light.loc,
                         thread_allocators.get());
  };
  light_jobs.push_back(run_async(create));
#endif
}

int Ri::add_area_light(Scene_Entity light)
{
  std::lock_guard<std::mutex> lock(area_light_mutex);
  area_lights.push_back(std::move(light));
  return area_lights.size() - 1;
}

void Ri::add_shapes(p_std::span<Shape_Scene_Entity> s)
{
  std::lock_guard<std::mutex> lock(shape_mutex);
  std::move(std::begin(s), std::end(s), std::back_inserter(shapes));
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

  std::lock_guard<std::mutex> lock(instance_definition_mutex);
  instance_definitions[def->name] = def;
}

void Ri::add_instance_uses(p_std::span<Instance_Scene_Entity> in)
{
  std::lock_guard<std::mutex> lock(instance_use_mutex);
  std::move(std::begin(in), std::end(in), std::back_inserter(instances));
}

Ri::~Ri()
{
}

// RI API Default Implementation
void Ri::ArchiveBegin(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "ArchiveBegin is unimplemented" << std::endl;
}

void Ri::ArchiveEnd(File_Loc loc)
{
  std::cout << "ArchiveEnd is unimplemented" << std::endl;
}

void Ri::AreaLightSource(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "AreaLightSource is unimplemented." << std::endl;
}

void Ri::Atmosphere(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Atmosphere is unimplemented" << std::endl;
}

void Ri::Attribute(const std::string &target, Parsed_Parameter_Vector params, File_Loc loc)
{
  for (Parsed_Parameter *p : params) {
    p->may_be_unused = true;
    graphics_state.rib_attributes[target].push_back(p);
  }

  if (target == "Ri") {
    for (Parsed_Parameter *p : params) {
      if (p->name == "Orientation")
        if (p->strings[0] == "inside")
          graphics_state.reverse_orientation = !graphics_state.reverse_orientation;
    }
  }
}

void Ri::AttributeBegin(File_Loc loc)
{
  VERIFY_WORLD("AttributeBegin");
  osl_parameters.clear();
  pushed_graphics_states.push_back(graphics_state);
  push_stack.push_back(std::make_pair('a', loc));
}

void Ri::AttributeEnd(File_Loc loc)
{
  VERIFY_WORLD("AttributeEnd");
  if (osl_parameters.size() > 0)
    osl_shader_group[_shader_id] = osl_parameters;

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
  else
    CHECK_EQ(push_stack.back().first, 'a');
  push_stack.pop_back();
}

void Ri::Basis(Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Basis is unimplemented" << std::endl;
}

void Ri::Begin(const std::string &name, File_Loc loc)
{
  std::cout << "Begin is unimplemented" << std::endl;
}

void Ri::Bound(float bound[6], File_Loc loc)
{
  std::cout << "Bound is unimplemented" << std::endl;
}

void Ri::Bxdf(const std::string &bxdf,
              const std::string &name,
              Parsed_Parameter_Vector params,
              File_Loc loc)
{
  VERIFY_WORLD("Bxdf");

  Parsed_Parameter *param = new Parsed_Parameter(loc);
  param->type = "string";
  param->name = "shader_type";
  param->add_string("bxdf");
  param->add_string(bxdf);
  param->add_string(name);
  params.push_back(param);

  Parameter_Dictionary dict(std::move(params));
  std::string material_id = dict.get_one_string("__materialid", "");
  if (material_id != "") {
    _shader_id = material_id;
    if (osl_shader_group.find(material_id) == osl_shader_group.end())
      osl_parameters.push_back(dict);
    else {
      graphics_state.current_material_name = material_id;
      graphics_state.current_material_index = -1;
    }
  }
  else if (string_startswith(bxdf, "Lama"))
    osl_parameters.push_back(dict);
}

void Ri::camera(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  VERIFY_OPTIONS("Camera");
  // Remove any class designator from the camera options
  // Specifically, Ri:<something>
  for (auto it = params.begin(); it != params.end(); it++) {
    std::vector<std::string> strings = split_string((*it)->name, ':');
    if (strings.size() > 1)
      // Could extract the size and check we have that many values below
      (*it)->name = strings[1];
  }

  auto options = _rib_state.options["Ri"];

  for (auto s : options) {
    auto *param = s.second;
    bool found = false;
    for (auto it = params.begin(); it != params.end(); it++)
      if ((*it)->name == param->name) {
        found = true;
        break;
      }

    if (!found)
      params.push_back(param);
  }

  Parameter_Dictionary dict(std::move(params));

  auto items = _rib_state.options.find("trace");
  if (items != _rib_state.options.end()) {
    auto world_origin = items->second.find("worldorigin");
    if (world_origin != items->second.end()) {
      if (world_origin->second->strings[0] == "worldoffset") {
        auto world_offset = items->second["worldoffset"];
        _rib_state.world_offset = make_float3(
            -world_offset->floats[0], -world_offset->floats[1], -world_offset->floats[2]);
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

void Ri::Clipping(float cnear, float cfar, File_Loc loc)
{
  std::cout << "Clipping is unimplemented" << std::endl;
}

void Ri::ClippingPlane(float x, float y, float z, float nx, float ny, float nz, File_Loc loc)
{
  std::cout << "ClippingPlane is unimplemented" << std::endl;
}

void Ri::Color(float r, float g, float b, File_Loc loc)
{
  std::cout << "Color is unimplemented" << std::endl;
}

void Ri::ConcatTransform(float transform[16], File_Loc loc)
{
  graphics_state.for_active_transforms([=](auto t) {
    ProjectionTransform projection = *(ProjectionTransform *)&transform[0];
    return t * projection_transpose(projection);
  });
}

void Ri::CoordinateSystem(std::string const &name, File_Loc loc)
{
  named_coordinate_systems[name] = graphics_state.ctm;
}

void Ri::CoordSysTransform(std::string const &name, File_Loc loc)
{
  if (named_coordinate_systems.find(name) != named_coordinate_systems.end())
    graphics_state.ctm = named_coordinate_systems[name];
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

  Parsed_Parameter *param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "vertices";
  for (int i = 0; i < polys.size(); ++i)
    param->add_int(polys[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nvertices";
  for (int i = 0; i < poly_count; ++i) {
    param->add_int(4);
  }
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->storage = Container_Type::Vertex;
  param->type = "point";
  param->name = "P";
  param->elem_per_item = 3;
  for (int i = 0; i < pts.size(); ++i) {
    param->add_float(pts[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->storage = Container_Type::Vertex;
  param->type = "float";
  param->name = "uv";
  param->elem_per_item = 2;
  for (int i = 0; i < uvs.size(); ++i)
    param->add_float(uvs[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nfaces";
  param->add_int(poly_count);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "bool";
  param->name = "smooth";
  param->add_bool(true);
  params.push_back(param);

  Shape("mesh", params, loc);
}

void Ri::CropWindow(float xmin, float xmax, float ymin, float ymax, File_Loc loc)
{
  if (xmin >= 0. && xmin <= 1. && xmin < xmax && xmin >= _crop_window[0] &&
      xmax >= 0. && xmax <= 1. && xmin < xmax && xmax <= _crop_window[1] &&
      ymin >= 0. && ymin <= 1. && ymin < ymax && ymin >= _crop_window[2] &&
      ymax >= 0. && ymax <= 1. && ymin < ymax && ymax <= _crop_window[3]) {
    _crop_window[0] = xmin;
    _crop_window[1] = xmax;
    _crop_window[2] = ymin;
    _crop_window[3] = ymax;

    if (_camera_name != "")
      _camera[_camera_name].crop_window = _crop_window;
  }
}

void Ri::Curves(const std::string &type,
                std::vector<int> nvertices,
                const std::string &wrap,
                Parsed_Parameter_Vector params,
                File_Loc loc)
{
  std::cout << "Curves is unimplemented" << std::endl;
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

  Parsed_Parameter *param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "vertices";
  for (int i = 0; i < polys.size(); ++i)
    param->add_int(polys[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nvertices";
  for (int i = 0; i < poly_count; ++i) {
    param->add_int(4);
  }
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->storage = Container_Type::Vertex;
  param->type = "point";
  param->name = "P";
  param->elem_per_item = 3;
  for (int i = 0; i < pts.size(); ++i) {
    param->add_float(pts[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->storage = Container_Type::Vertex;
  param->type = "float";
  param->name = "uv";
  param->elem_per_item = 2;
  for (int i = 0; i < uvs.size(); ++i)
    param->add_float(uvs[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nfaces";
  param->add_int(poly_count);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "bool";
  param->name = "smooth";
  param->add_bool(true);
  params.push_back(param);

  Shape("mesh", params, loc);
}

void Ri::Declare(const std::string &name, const std::string &declaration, File_Loc loc)
{
  std::cout << "Declare is unimplemented" << std::endl;
}

void Ri::DepthOfField(float fstop, float focallength, float focaldistance, File_Loc loc)
{
  std::cout << "DepthOfField is unimplemented" << std::endl;
}

void Ri::Detail(float bound[6], File_Loc loc)
{
  std::cout << "Detail is unimplemented" << std::endl;
}

void Ri::DetailRange(float offlow, float onlow, float onhigh, float offhigh, File_Loc loc)
{
  std::cout << "DetailRange is unimplemented" << std::endl;
}

void Ri::Disk(
    float height, float radius, float thetamax, Parsed_Parameter_Vector params, File_Loc loc)
{
  VERIFY_WORLD("Shape");

  Parsed_Parameter *param = new Parsed_Parameter(loc);
  param->type = "float";
  param->name = "radius";
  param->add_float(radius);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "float";
  param->name = "height";
  param->add_float(height);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "float";
  param->name = "phimax";
  param->add_float(thetamax);
  params.push_back(param);

  Shape("disk", params, loc);
}

void Ri::Displacement(const std::string &displace,
                      const std::string &name,
                      Parsed_Parameter_Vector params,
                      File_Loc loc)
{
  std::cout << "Displacement is unimplemented" << std::endl;
}

void Ri::Display(const std::string &name,
                 const std::string &type,
                 const std::string &mode,
                 Parsed_Parameter_Vector params,
                 File_Loc loc)
{
  VERIFY_OPTIONS("Film");
  // If the first char of `name' is a '+', then it's an additional
  // channel to render. Ignore it for now.
  if (name[0] == '+') {
    fprintf(stdout, "Ignoring additional display channel, %s\n", name.c_str());
    return;
  }

  Parsed_Parameter_Vector new_params;

  Parsed_Parameter *param = new Parsed_Parameter(loc);
  param->type = "string";
  param->name = "filename";
  param->add_string(name);
  new_params.push_back(param);

  auto resolution = _rib_state.options["Ri"]["FormatResolution"]->ints;
  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "xresolution";
  param->add_int((int)(resolution[0]));
  new_params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "yresolution";
  param->add_int((int)(resolution[1]));
  new_params.push_back(param);

  Parameter_Dictionary dict(std::move(new_params));
  std::string camera_name = dict.get_one_string("camera_name", "");
  if (camera_name.empty())
    camera_name = _camera_name;
  else
    _camera_name = camera_name;

  film = Display_Scene_Entity("rgb", std::move(dict), loc, camera_name);
}

void Ri::DisplayChannel(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "DisplayChannel " << name << " is unimplemented" << std::endl;
}

void Ri::DisplayFilter(const std::string &name,
                       const std::string &type,
                       Parsed_Parameter_Vector params,
                       File_Loc loc)
{
  std::cout << "DisplayFilter " << name << " is unimplemented" << std::endl;
}

void Ri::Else(File_Loc loc)
{
  std::cout << "Else is unimplemented" << std::endl;
}

void Ri::ElseIf(const std::string &condition, File_Loc loc)
{
  std::cout << "ElseIf is unimplemented" << std::endl;
}

void Ri::End(File_Loc loc)
{
  std::cout << "End is unimplemented" << std::endl;
}

void Ri::Exposure(float gain, float gamma, File_Loc loc)
{
  std::cout << "Exposure is unimplemented" << std::endl;
}

void Ri::Exterior(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Exterior is unimplemented" << std::endl;
}

void Ri::Format(int xresolution, int yresolution, float pixelaspectratio, File_Loc loc)
{
  std::cout << "Format is unimplemented" << std::endl;
}

void Ri::FrameAspectRatio(float frameratio, File_Loc loc)
{
  std::cout << "FrameAspectRatio is unimplemented" << std::endl;
}

void Ri::FrameBegin(int number, File_Loc loc)
{
  std::cout << "FrameBegin is unimplemented" << std::endl;
}

void Ri::FrameEnd(File_Loc loc)
{
  std::cout << "FrameEnd is unimplemented" << std::endl;
}

void Ri::GeneralPolygon(std::vector<int> nvertices, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "GeneralPolygon is unimplemented" << std::endl;
}

void Ri::GeometricApproximation(const std::string &type, float value, File_Loc loc)
{
  std::cout << "GeometricApproximation is unimplemented" << std::endl;
}

void Ri::Geometry(const std::string &type, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Geometry is unimplemented" << std::endl;
}

void Ri::Hider(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Hider is unimplemented" << std::endl;
}

void Ri::Hyperboloid(
    Point3f point1, Point3f point2, float thetamax, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Hyperboloid is unimplemented" << std::endl;
}

void Ri::Identity(File_Loc loc)
{
  graphics_state.for_active_transforms([](auto t) { return projection_identity(); });
}

void Ri::IfBegin(const std::string &condition, File_Loc loc)
{
  std::cout << "IfBegin is unimplemented" << std::endl;
}

void Ri::IfEnd(File_Loc loc)
{
  std::cout << "IfEnd is unimplemented" << std::endl;
}

void Ri::Illuminate(const std::string &light, bool onoff, File_Loc loc)
{
  std::cout << "Illuminate is unimplemented" << std::endl;
}

void Ri::Imager(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Imager is unimplemented" << std::endl;
}

void Ri::Integrator(const std::string &type,
                    const std::string &name,
                    Parsed_Parameter_Vector params,
                    File_Loc loc)
{
  std::cout << "Integrator is unimplemented" << std::endl;
}

void Ri::Interior(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
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
    if ((*it)->name == "sides") {
      double_sided = (*it)->floats[0] == 2;
      break;
    }

  Parsed_Parameter *param;
  // Lights are single-sided by default
  if (double_sided) {
    param = new Parsed_Parameter(loc);
    param->type = "bool";
    param->name = "twosided";
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
    if (name == "PxrCylinderLight") {
      // A cylinder light lies along the X-axis, and is half the default cylinder size
      Rotate(90, 0, 1, 0, loc);
      Rotate(90, 0, 0, 1, loc);
      Cylinder(0.5, -0.5, 0.5, 360, params, loc);

      _cylinder_light_material = new Parsed_Parameter_Vector;
      Parsed_Parameter *param = new Parsed_Parameter(loc);
      param->type = "string";
      param->name = "__materialid";
      param->may_be_unused = true;
      param->add_string(handle);
      _cylinder_light_material->push_back(param);

      // Check if the light has a materialid parameter
      if (entity.parameters.get_one_string("__materialid", "") == "") {
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
    else
      active_instance_definition->entity.lights.push_back(std::move(entity));
  }
  else {
    std::string material_id = entity.parameters.get_one_string("__materialid", "");
    if (material_id != "") {
      size_t len = entity.parameters.get_parameter_vector().size();
      // If we're processing the light instance (only occurs with a mesh light),
      // we need to create an emission shader that's picked up by the geometry
      // at a later point.
      if (len <= 2) {
        auto it = osl_shader_group.find(material_id);
        if (it != osl_shader_group.end())
          osl_shader_group.erase(it);
        auto light_params = _lights[handle].parameters.get_parameter_vector();
        Parameter_Dictionary light_dict(light_params, graphics_state.light_attributes);
        float strength = 1.0f;

        float exposure = light_dict.get_one_float("exposure", 1.0);
        strength *= exp2(exposure);

        float intensity = light_dict.get_one_float("intensity", 1.0);
        strength *= intensity;
        param = new Parsed_Parameter(loc);
        param->type = "float";
        param->name = "strength";
        param->may_be_unused = true;
        param->add_float(strength);
        light_params.push_back(param);
        Bxdf(name, handle, light_params, loc);
        graphics_state.current_material_name = material_id;
        graphics_state.current_material_index = -1;
      }
      else
        _lights[handle] = std::move(entity);
    }
  }
}

void Ri::LightSource(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "LightSource is unimplemented." << std::endl;
}

void Ri::MakeBrickMap(std::vector<std::string> ptcnames,
                      const std::string &bkmname,
                      Parsed_Parameter_Vector params,
                      File_Loc loc)
{
  std::cout << "MakeBrickMap is unimplemented" << std::endl;
}

void Ri::MakeCubeFaceEnvironment(const std::string &px,
                                 const std::string &nx,
                                 const std::string &py,
                                 const std::string &ny,
                                 const std::string &pz,
                                 const std::string &nz,
                                 const std::string &text,
                                 float fov,
                                 const std::string &filt,
                                 float swidth,
                                 float twidth,
                                 Parsed_Parameter_Vector params,
                                 File_Loc loc)
{
  std::cout << "MakeCubeFaceEnvironment is unimplemented" << std::endl;
}

void Ri::MakeLatLongEnvironment(const std::string &picturename,
                                const std::string &texturename,
                                const std::string &filt,
                                float swidth,
                                float twidth,
                                Parsed_Parameter_Vector params,
                                File_Loc loc)
{
  std::cout << "MakeLatLongEnvironment is unimplemented" << std::endl;
}

void Ri::MakeShadow(const std::string &picturename,
                    const std::string &texturename,
                    Parsed_Parameter_Vector params,
                    File_Loc loc)
{
  std::cout << "MakeShadow is unimplemented" << std::endl;
}

void Ri::MakeTexture(const std::string &picturename,
                     const std::string &texturename,
                     const std::string &swrap,
                     const std::string &twrap,
                     const std::string &filt,
                     float swidth,
                     float twidth,
                     Parsed_Parameter_Vector params,
                     File_Loc loc)
{
  std::cout << "MakeTexture is unimplemented" << std::endl;
}

void Ri::Matte(bool onoff, File_Loc loc)
{
  std::cout << "Matte is unimplemented" << std::endl;
}

void Ri::MotionBegin(std::vector<float> times, File_Loc loc)
{
  std::cout << "MotionBegin is unimplemented" << std::endl;
}

void Ri::MotionEnd(File_Loc loc)
{
  std::cout << "MotionEnd is unimplemented" << std::endl;
}

void Ri::NuPatch(int nu,
                 int uorder,
                 float uknot[],
                 float umin,
                 float umax,
                 int nv,
                 int vorder,
                 float vknot[],
                 float vmin,
                 float vmax,
                 Parsed_Parameter_Vector params,
                 File_Loc loc)
{
  std::cout << "NuPatch is unimplemented" << std::endl;
}

void Ri::ObjectBegin(const std::string &name, File_Loc loc)
{
  VERIFY_WORLD("ObjectBegin");
  pushed_graphics_states.push_back(graphics_state);

  push_stack.push_back(std::make_pair('o', loc));

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
  else
    CHECK_EQ(push_stack.back().first, 'o');
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

  if (_cylinder_light_material) {
    Light("PxrMeshLight",
          (*_cylinder_light_material)[0]->strings[0],
          *_cylinder_light_material,
          loc);
    osl_shader_group[_shader_id] = osl_parameters;
    delete _cylinder_light_material;
    _cylinder_light_material = nullptr;
  }

  Mapped_Parameter_Dictionary dict;
  for (auto &x : graphics_state.rib_attributes) {
    Parameter_Dictionary d(x.second);

    dict.emplace(x.first, d);
  }

  if (dict.find("identifier") == dict.end()) {
    Parsed_Parameter_Vector params;
    Parsed_Parameter *param = new Parsed_Parameter(loc);
    param->type = "string";
    param->name = "name";
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
  instance_uses[name].push_back(Instance_Scene_Entity(name,
                                                      loc,
                                                      graphics_state.current_material_name,
                                                      std::move(dict) /* RenderMan*/,
                                                      render_from_instance));
}

void Ri::Option(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
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

void Ri::Opacity(Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Opacity is unimplemented" << std::endl;
}

void Ri::Orientation(const std::string &orientation, File_Loc loc)
{
  std::cout << "Orientation is unimplemented" << std::endl;
}

void Ri::Paraboloid(float rmax,
                    float zmin,
                    float zmax,
                    float thetamax,
                    Parsed_Parameter_Vector params,
                    File_Loc loc)
{
  std::cout << "Paraboloid is unimplemented" << std::endl;
}

void Ri::Patch(const std::string &type, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Patch is unimplemented" << std::endl;
}

void Ri::PatchMesh(const std::string &type,
                   int nu,
                   const std::string &uwrap,
                   int nv,
                   const std::string &vwrap,
                   Parsed_Parameter_Vector params,
                   File_Loc loc)
{
  std::cout << "PatchMesh is unimplemented" << std::endl;
}

void Ri::Pattern(const std::string &name,
                 const std::string &handle,
                 Parsed_Parameter_Vector params,
                 File_Loc loc)
{
  VERIFY_WORLD("Pattern");

  Parsed_Parameter *param = new Parsed_Parameter(loc);
  param->type = "string";
  param->name = "shader_type";
  param->add_string("pattern");
  param->add_string(name);
  param->add_string(handle);
  params.push_back(param);

  Parameter_Dictionary dict(std::move(params));
  osl_parameters.push_back(dict);
}

void Ri::Perspective(float fov, File_Loc loc)
{
  std::cout << "Perspective is unimplemented" << std::endl;
}

void Ri::PixelFilter(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "PixelFilter is unimplemented." << std::endl;
}

void Ri::PixelSampleImager(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "PixelSampleImager is unimplemented" << std::endl;
}

void Ri::PixelSamples(float xsamples, float ysamples, File_Loc loc)
{
  std::cout << "PixelSamples is unimplemented" << std::endl;
}

void Ri::PixelVariance(float variance, File_Loc loc)
{
  std::cout << "PixelVariance is unimplemented" << std::endl;
}

void Ri::Points(int npoints, Parsed_Parameter_Vector params, File_Loc loc)
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
  for (int i = 0; single_only && i < n_loops.size(); i++)
    single_only = ((int)n_loops[i] == 1);

  if (single_only)
    PointsPolygons(n_vertices, vertices, params, loc);
}

void Ri::PointsPolygons(std::vector<int> n_vertices,
                        std::vector<int> vertices,
                        Parsed_Parameter_Vector params,
                        File_Loc loc)
{
  Parsed_Parameter *param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "vertices";
  for (int i = 0; i < vertices.size(); ++i)
    param->add_int(vertices[i]);
  params.push_back(param);

  int nfaces = n_vertices.size();
  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nvertices";
  for (int i = 0; i < n_vertices.size(); ++i)
    param->add_int(n_vertices[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nfaces";
  param->add_int(nfaces);
  params.push_back(param);

  // Fix attributes whose storage class couldn't be determined at parsing
  // time
  for (auto &p : params) {
    if ((p->type == "point" || p->type == "normal" || p->type == "color") &&
        p->storage == Container_Type::Constant && p->floats.size() > 3) {
      int num_vals = p->floats.size() / 3;
      if (num_vals == n_vertices.size())
        p->storage = Container_Type::Uniform;
      else if (num_vals == vertices.size())
        p->storage = Container_Type::Vertex;
    }
  }

  bool has_varying_normals = false;
  for (auto &p : params) {
    if (p->type == "normal") {
      if (p->storage == Container_Type::FaceVarying || p->storage == Container_Type::Varying ||
          p->storage == Container_Type::Vertex)
        has_varying_normals = true;
      break;
    }
  }

  if (has_varying_normals) {
    param = new Parsed_Parameter(loc);
    param->type = "bool";
    param->name = "smooth";
    param->add_bool(true);
    params.push_back(param);
  }
  Shape("mesh", params, loc);
}

void Ri::Polygon(int nvertices, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Polygon is unimplemented" << std::endl;
}

void Ri::ProcDelayedReadArchive(Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "ProcDelayedReadArchive is unimplemented" << std::endl;
}

void Ri::ProcDynamicLoad(Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "ProcDynamicLoad is unimplemented" << std::endl;
}

void Ri::Procedural(const std::string &proc_name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Procedural is unimplemented" << std::endl;
}

void Ri::Procedural2(const std::string &proc_name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Procedural2 is unimplemented" << std::endl;
}

void Ri::ProcFree(Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "ProcFree is unimplemented" << std::endl;
}

void Ri::ProcRunProgram(Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "ProcRunProgram is unimplemented" << std::endl;
}

void Ri::Projection(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
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

void Ri::Quantize(
    const std::string &type, int one, int min, int max, float ditheramplitude, File_Loc loc)
{
  std::cout << "Quantize is unimplemented" << std::endl;
}

void Ri::ReadArchive(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "ReadArchive is unimplemented" << std::endl;
}

void Ri::RelativeDetail(float relativedetail, File_Loc loc)
{
  std::cout << "RelativeDetail is unimplemented" << std::endl;
}

void Ri::Resource(const std::string &handle,
                  const std::string &type,
                  Parsed_Parameter_Vector params,
                  File_Loc loc)
{
  std::cout << "Resource is unimplemented" << std::endl;
}

void Ri::ResourceBegin(File_Loc loc)
{
  std::cout << "ResourceBegin is unimplemented" << std::endl;
}

void Ri::ResourceEnd(File_Loc loc)
{
  std::cout << "ResourceEnd is unimplemented" << std::endl;
}

void Ri::ReverseOrientation(File_Loc loc)
{
  VERIFY_WORLD("ReverseOrientation");
  graphics_state.reverse_orientation = !graphics_state.reverse_orientation;
}

void Ri::Rotate(float angle, float ax, float ay, float az, File_Loc loc)
{
  graphics_state.for_active_transforms(
      [=](auto t) { return t * transform_rotate(DEG2RADF(angle), make_float3(ax, ay, az)); });
}

void Ri::Scale(float sx, float sy, float sz, File_Loc loc)
{
  graphics_state.for_active_transforms(
      [=](auto t) { return t * transform_scale(make_float3(sx, sy, sz)); });
}

void Ri::ScopedCoordinateSystem(const std::string &, File_Loc loc)
{
  std::cout << "ScopedCoordinateSystem is unimplemented" << std::endl;
}

void Ri::ScreenWindow(float left, float right, float bottom, float top, File_Loc loc)
{
  std::cout << "ScreenWindow is unimplemented" << std::endl;
}

void Ri::ShadingInterpolation(const std::string &type, File_Loc loc)
{
  std::cout << "ShadingInterpolation is unimplemented" << std::endl;
}

void Ri::ShadingRate(float size, File_Loc loc)
{
  std::cout << "ShadingRate is unimplemented" << std::endl;
}

void Ri::Shutter(float opentime, float closetime, File_Loc loc)
{
  std::cout << "Shutter is unimplemented" << std::endl;
}

void Ri::Sides(int nsides, File_Loc loc)
{
  Parsed_Parameter *param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "sides";
  param->add_float(nsides);
  graphics_state.shape_attributes.push_back(param);
}

void Ri::Skew(
    float angle, float dx1, float dy1, float dz1, float dx2, float dy2, float dz2, File_Loc loc)
{
  std::cout << "Skew is unimplemented" << std::endl;
}

void Ri::SolidBegin(const std::string &type, File_Loc loc)
{
  std::cout << "SolidBegin is unimplemented" << std::endl;
}

void Ri::SolidEnd(File_Loc loc)
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
  int n_stacks = 25;

  // Create top vertices
  float phi = start_phi;
  for (int j = 0; j < n_slices; ++j) {
    float theta = thetamax * float(j) / float(n_slices);
    float3 dir = spherical_to_direction(-theta, phi);
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
      float3 dir = spherical_to_direction(-theta, phi);
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
    float3 dir = spherical_to_direction(-theta, phi);
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
  for (int j = 0; j < n_stacks - 1; j++) {
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

  for (int j = 0; j < n_stacks - 1; j++) {
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

  Parsed_Parameter *param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "vertices";
  for (int i = 0; i < polys.size(); ++i)
    param->add_int(polys[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nvertices";
  for (int i = 0; i < poly_counts.size(); ++i) {
    param->add_int(poly_counts[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->storage = Container_Type::Vertex;
  param->type = "point";
  param->name = "P";
  param->elem_per_item = 3;
  for (int i = 0; i < pts.size(); ++i) {
    param->add_float(pts[i]);
  }
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->storage = Container_Type::FaceVarying;
  param->type = "float";
  param->name = "uv";
  param->elem_per_item = 2;
  for (int i = 0; i < uvs.size(); ++i) {
    param->add_float(uvs[i].x);
    param->add_float(uvs[i].y);
  }
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nfaces";
  param->add_int(poly_count);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "bool";
  param->name = "smooth";
  param->add_bool(true);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->storage = Container_Type::Varying;
  param->type = "normal";
  param->name = "N";
  param->elem_per_item = 3;
  for (int i = 0; i < norms.size(); ++i)
    param->add_float(norms[i]);
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
  Parsed_Parameter *param = new Parsed_Parameter(loc);
  param->type = "string";
  param->name = "scheme";
  param->add_string(scheme);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nfaces";
  param->add_int(nfaces);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "vertices";
  for (int i = 0; i < vertices.size(); ++i)
    param->add_int(vertices[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nvertices";
  for (int i = 0; i < n_vertices.size(); ++i)
    param->add_int(n_vertices[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "string";
  param->name = "tags";
  for (int i = 0; i < tags.size(); ++i)
    param->add_string(tags[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "nargs";
  for (int i = 0; i < nargs.size(); ++i)
    param->add_int(nargs[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "int";
  param->name = "intargs";
  for (int i = 0; i < intargs.size(); ++i)
    param->add_int(intargs[i]);
  params.push_back(param);

  param = new Parsed_Parameter(loc);
  param->type = "float";
  param->name = "floatargs";
  for (int i = 0; i < floatargs.size(); ++i)
    param->add_float(floatargs[i]);
  params.push_back(param);

  Shape("subdivision_mesh", params, loc);
}

void Ri::Surface(const std::string &name, Parsed_Parameter_Vector params, File_Loc loc)
{
  std::cout << "Surface is unimplemented" << std::endl;
}

void Ri::System(const std::string &cmd, File_Loc loc)
{
  std::cout << "System is unimplemented" << std::endl;
}

void Ri::Texture(const std::string &name,
                 const std::string &type,
                 const std::string &texname,
                 Parsed_Parameter_Vector params,
                 File_Loc loc)
{
  std::cout << "Texture is unimplemented." << std::endl;
}

void Ri::TextureCoordinates(
    float s1, float t1, float s2, float t2, float s3, float t3, float s4, float t4, File_Loc loc)
{
  std::cout << "TextureCoordinates is unimplemented" << std::endl;
}

void Ri::Torus(float majorrad,
               float minorrad,
               float phimin,
               float phimax,
               float thetamax,
               Parsed_Parameter_Vector params,
               File_Loc loc)
{
  std::cout << "Torus is unimplemented" << std::endl;
}

void Ri::transform(float transform[16], File_Loc loc)
{
  graphics_state.for_active_transforms([=](auto t) {
    // Stomp the current transform
    ProjectionTransform projection = *(ProjectionTransform *)&transform[0];
    return projection_transpose(projection);
  });
}

void Ri::TransformBegin(File_Loc loc)
{
  pushed_graphics_states.push_back(graphics_state);
  push_stack.push_back(std::make_pair('t', loc));
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
  else
    CHECK_EQ(push_stack.back().first, 't');
  push_stack.pop_back();
}

void Ri::Translate(float dx, float dy, float dz, File_Loc loc)
{
  graphics_state.for_active_transforms(
      [=](auto t) { return t * transform_translate(make_float3(dx, dy, dz)); });
}

void Ri::TrimCurve(int nloops,
                   int ncurves[],
                   int order[],
                   float knot[],
                   float min[],
                   float max[],
                   int n[],
                   float u[],
                   float v[],
                   float w[],
                   File_Loc loc)
{
  std::cout << "TrimCurve is unimplemented" << std::endl;
}

void Ri::WorldBegin(File_Loc loc)
{
  VERIFY_OPTIONS("WorldBegin");
  // Reset graphics state for _WorldBegin_
  current_block = Block_State::World_Block;
  for (int i = 0; i < Max_Transforms; ++i)
    graphics_state.ctm[i] = projection_identity();
  named_coordinate_systems["world"] = graphics_state.ctm;

  Parsed_Parameter_Vector params;
  Parsed_Parameter *param;

  auto options = _rib_state.options["Ri"];
  auto *opt_param = options["PixelFilterName"];
  if (opt_param != nullptr) {
    std::string name = opt_param->strings[0];
    float xradius = 2, yradius = 2;

    auto *opt_width = options["PixelFilterWidth"];
    if (opt_width != nullptr) {
      xradius = opt_width->floats[0];
      yradius = opt_width->floats[1];
    }

    param = new Parsed_Parameter(loc);
    param->type = "float";
    param->name = "xradius";
    param->add_float(xradius);
    params.push_back(param);

    param = new Parsed_Parameter(loc);
    param->type = "float";
    param->name = "yradius";
    param->add_float(yradius);
    params.push_back(param);

    Parameter_Dictionary dict(std::move(params));
    filter = Scene_Entity(name, std::move(dict), loc);
  }

  opt_param = options["PixelVariance"];
  if (opt_param != nullptr) {
    float pixel_var = opt_param->floats[0];

    param = new Parsed_Parameter(loc);
    param->type = "float";
    param->name = "pixel_variance";
    param->add_float(pixel_var);
    sampler.parameters.push_back(param);
  }

  options = _rib_state.options["searchpath"];
  opt_param = options["shader"];
  std::string paths = opt_param->strings[0] + ":";

  options = _rib_state.options["hider"];
  opt_param = options["minsamples"];
  if (opt_param != nullptr) {
    int min_samples = opt_param->ints[0];
    opt_param = options["maxsamples"];
    int max_samples = opt_param->ints[0];

    param = new Parsed_Parameter(loc);
    param->type = "int";
    param->name = "min_samples";
    param->add_int(min_samples);
    sampler.parameters.push_back(param);

    param = new Parsed_Parameter(loc);
    param->type = "int";
    param->name = "max_samples";
    param->add_int(max_samples);
    sampler.parameters.push_back(param);
  }
}

void Ri::WorldEnd(File_Loc loc)
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

  param = new Parsed_Parameter(File_Loc());
  param->type = "string";
  param->name = "shader_default";
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
  if (active_instance_definition)
    active_instance_definition->entity.shapes.push_back(std::move(entity));
  else
    shapes.push_back(std::move(entity));
}

CCL_NAMESPACE_END
