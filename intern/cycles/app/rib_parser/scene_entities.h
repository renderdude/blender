#ifndef SCENE_ENTITIES_H
#define SCENE_ENTITIES_H

#include <algorithm>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "error.h"
#include "param_dict.h"
#include "util/log.h"
#include "util/projection.h"
#include "util/transform.h"

CCL_NAMESPACE_BEGIN

// Max_Transforms Definition
constexpr int Max_Transforms = 2;

// Transform_Set Definition
struct Transform_Set {
  // Transform_Set Public Methods
  ProjectionTransform &operator[](int i)
  {
    CHECK_GE(i, 0);
    CHECK_LT(i, Max_Transforms);
    return t[i];
  }
  const ProjectionTransform &operator[](int i) const
  {
    CHECK_GE(i, 0);
    CHECK_LT(i, Max_Transforms);
    return t[i];
  }
  friend Transform_Set inverse(const Transform_Set &ts)
  {
    Transform_Set tInv;
    for (int i = 0; i < Max_Transforms; ++i)
      tInv.t[i] = projection_inverse(ts.t[i]);
    return tInv;
  }
  bool is_animated() const
  {
    for (int i = 0; i < Max_Transforms - 1; ++i)
      if (t[i] != t[i + 1])
        return true;
    return false;
  }

 private:
  ProjectionTransform t[Max_Transforms];
};

static constexpr int Start_Transform_Bits = 1 << 0;
static constexpr int End_Transform_Bits = 1 << 1;
static constexpr int All_Transforms_Bits = (1 << Max_Transforms) - 1;

// Graphics_State Definition
struct Graphics_State {
  template<typename F> void for_active_transforms(F func)
  {
    for (int i = 0; i < Max_Transforms; ++i)
      if (active_transform_bits & (1 << i))
        ctm[i] = func(ctm[i]);
  }

  // Graphics_State Public Members
  std::string current_inside_medium, current_outside_medium;

  int current_material_index = 0;
  std::string current_material_name;

  std::string area_light_name;
  Parameter_Dictionary area_light_params;
  File_Loc area_light_loc;

  Parsed_Parameter_Vector shape_attributes;
  Parsed_Parameter_Vector light_attributes;
  Parsed_Parameter_Vector material_attributes;
  Parsed_Parameter_Vector medium_attributes;
  Parsed_Parameter_Vector texture_attributes;
  // RenderMan
  std::unordered_map<std::string, Parsed_Parameter_Vector> rib_attributes;

  bool reverse_orientation = false;
  Transform_Set ctm;
  uint32_t active_transform_bits = All_Transforms_Bits;
  float transform_start_time = 0, transform_end_time = 1;
};

class Camera_Transform {
 public:
  /// @name Initialization
  ///@{
  Camera_Transform() = default;
  explicit Camera_Transform(const ProjectionTransform &world_from_camera, float3 world_offset)
  {
    /*
     case Rendering_Coordinate_System::Camera: {
       // Compute __world_from_render_ for camera-space rendering
       Real tMid = (world_from_camera.start_time + world_from_camera.end_time) / 2;
       _world_from_render = world_from_camera.interpolate(tMid);
       break;
     }
        case Rendering_Coordinate_System::Camera_World: {
          // Compute __world_from_render_ for camera-world space rendering
          Real tMid = (world_from_camera.start_time + world_from_camera.end_time) / 2;
          Point3f pCamera = world_from_camera(Point3f(0, 0, 0), tMid);
          _world_from_render = translate(Vector3f(pCamera));
          break;
        }
        case Rendering_Coordinate_System::World: {
          // Compute __world_from_render_ for world-space rendering
        */
    _world_from_render = ProjectionTransform(transform_translate(world_offset));
    //         break;
    //       }
    // Compute __render_from_camera_ transformation
    ProjectionTransform render_from_world = projection_inverse(_world_from_render);
    _render_from_camera = render_from_world * world_from_camera;
  }

  ///@}

  /// @name Access
  ///@{
  const ProjectionTransform &render_from_camera() const
  {
    return _render_from_camera;
  }

  const ProjectionTransform &world_from_render() const
  {
    return _world_from_render;
  }

  ProjectionTransform render_from_world() const
  {
    return projection_inverse(_world_from_render);
  }

  ProjectionTransform camera_from_render(float time = 0.f) const
  {
    return projection_inverse(_render_from_camera);
  }

  ProjectionTransform camera_from_world(float time = 0.f) const
  {
    return projection_inverse(_world_from_render * _render_from_camera);
  }

 private:
  ProjectionTransform _render_from_camera;
  ProjectionTransform _world_from_render;

};  // end of class Camera_Transform

struct Scene_Entity {
  // Scene_Entity Public Methods
  Scene_Entity() = default;
  Scene_Entity(const std::string &name, Parameter_Dictionary parameters, File_Loc loc)
      : name(name), loc(loc), parameters(parameters)
  {
  }

  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[ Scene_Entity name: " << name;
    ss << " parameters: " << parameters.to_string();
    ss << " loc: " << loc.to_string() << " ]";
    return ss.str();
  }

  // Scene_Entity Public Members
  std::string name;
  File_Loc loc;
  Parameter_Dictionary parameters;
};

struct Transformed_Scene_Entity : public Scene_Entity {
  Transformed_Scene_Entity() = default;
  Transformed_Scene_Entity(const std::string &name,
                           Parameter_Dictionary parameters,
                           File_Loc loc,
                           const ProjectionTransform &render_from_object)
      : Scene_Entity(name, parameters, loc), render_from_object(render_from_object)
  {
  }

  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[ Transformed_Scene_Entity name: " << name;
    ss << " parameters: " << parameters.to_string();
    ss << " loc: " << loc.to_string();
    //         ss << " render_from_object: " << render_from_object << " ]";
    return ss.str();
  }

  ProjectionTransform render_from_object;
};

// Camera_Scene_Entity Definition
struct Camera_Scene_Entity : public Scene_Entity {
  // Camera_Scene_Entity Public Methods
  Camera_Scene_Entity() = default;
  Camera_Scene_Entity(const std::string &name,
                      Parameter_Dictionary parameters,
                      File_Loc loc,
                      const Camera_Transform &camera_transform,
                      const std::string &medium)
      : Scene_Entity(name, parameters, loc), camera_transform(camera_transform), medium(medium)
  {
  }

  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[ Camera_Scene_Entity name: " << name;
    ss << " parameters: " << parameters.to_string();
    ss << " loc: " << loc.to_string();
    ss << " camera_transform: ";  //<< camera_transform;
    ss << " medium: " << medium << " ]";
    return ss.str();
  }

  Camera_Transform camera_transform;
  std::vector<float> crop_window;
  std::string medium;
};

// Display_Scene_Entity Definition
struct Display_Scene_Entity : public Scene_Entity {
  // Display_Scene_Entity Public Methods
  Display_Scene_Entity() = default;
  Display_Scene_Entity(const std::string &name,
                       Parameter_Dictionary parameters,
                       File_Loc loc,
                       const std::string &camera_name)
      : Scene_Entity(name, parameters, loc), camera_name(camera_name)
  {
  }

  std::string camera_name;
};

struct Shape_Scene_Entity : public Scene_Entity {
  Shape_Scene_Entity() = default;
  Shape_Scene_Entity(const std::string &name,
                     Parameter_Dictionary parameters,
                     File_Loc loc,
                     const ProjectionTransform *render_from_object,
                     const ProjectionTransform *object_from_render,
                     int light_index,
                     Graphics_State gstate)
      : Scene_Entity(name, parameters, loc),
        render_from_object(render_from_object),
        object_from_render(object_from_render),
        light_index(light_index),
        graphics_state(gstate)
  {
  }

  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[ Shape_Scene_Entity name: " << name;
    ss << " parameters: " << parameters.to_string();
    ss << " loc: " << loc.to_string();
    ss << " render_from_object: ";  // << render_from_object;
    ss << " object_from_render: ";  // << object_from_render;
    ss << " light_index: " << light_index << "]";
    return ss.str();
  }

  const ProjectionTransform *render_from_object = nullptr, *object_from_render = nullptr;
  int light_index = -1;
  Graphics_State graphics_state;
};

struct Animated_Shape_Scene_Entity : public Transformed_Scene_Entity {
  Animated_Shape_Scene_Entity() = default;
  Animated_Shape_Scene_Entity(const std::string &name,
                              Parameter_Dictionary parameters,
                              File_Loc loc,
                              const ProjectionTransform &render_from_object,
                              const ProjectionTransform *identity,
                              bool reverse_orientation,
                              int material_index,
                              const std::string &material_name,
                              int light_index,
                              const std::string &inside_medium,
                              const std::string &outside_medium)
      : Transformed_Scene_Entity(name, parameters, loc, render_from_object),
        identity(identity),
        reverse_orientation(reverse_orientation),
        material_index(material_index),
        material_name(material_name),
        light_index(light_index),
        inside_medium(inside_medium),
        outside_medium(outside_medium)
  {
  }

  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[ Animated_Shape_Scene_Entity name: " << name;
    ss << " parameters: " << parameters.to_string();
    ss << " loc: " << loc.to_string();
    ss << " render_from_object: ";  // << render_from_object;
    ss << " object_from_render: ";  // << object_from_render;
    ss << " reverse_orientation: " << reverse_orientation;
    ss << " material_index: " << material_index;
    ss << " material_name: " << material_name;
    ss << " light_index: " << light_index;
    ss << " inside_medium: " << inside_medium;
    ss << " outside_medium: " << outside_medium << "]";
    return ss.str();
  }

  const ProjectionTransform *identity = nullptr;
  bool reverse_orientation = false;
  int material_index;  // one of these two...  std::variant?
  std::string material_name;
  int light_index = -1;
  std::string inside_medium, outside_medium;
};

struct Light_Scene_Entity : public Scene_Entity {
  Light_Scene_Entity() = default;
  Light_Scene_Entity(const std::string &name,
                     const std::string &light_type,
                     Parameter_Dictionary parameters,
                     File_Loc loc,
                     const ProjectionTransform *render_from_light)
      : Scene_Entity(name, parameters, loc),
        render_from_light(render_from_light),
        light_type(light_type)
  {
  }

  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[ Light_Scene_Entity name: " << name;
    ss << " parameters: " << parameters.to_string();
    ss << " loc: " << loc.to_string();
    ss << " render_from_object: ";  // << render_from_object;
    ss << " light_type: " << light_type;
    return ss.str();
  }

  const ProjectionTransform *render_from_light = nullptr;
  std::string light_type;
};

struct Instance_Definition_Scene_Entity {
  Instance_Definition_Scene_Entity() = default;
  Instance_Definition_Scene_Entity(const std::string &name, File_Loc loc) : name(name), loc(loc)
  {
  }

  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[ Instance_Definition_Scene_Entity name: " << name;
    ss << " loc: " << loc.to_string();
    ss << " lights: ";           // << lights;
    ss << " shapes: ";           // << shapes;
    ss << " animated_shapes: ";  // << animated_shapes;
    ss << " ]";
    return ss.str();
  }

  std::string name;
  File_Loc loc;
  std::vector<Light_Scene_Entity> lights;
  std::vector<Shape_Scene_Entity> shapes;
  std::vector<Animated_Shape_Scene_Entity> animated_shapes;
};

using Medium_Scene_Entity = Transformed_Scene_Entity;
using Texture_Scene_Entity = Transformed_Scene_Entity;

struct Instance_Scene_Entity {
  Instance_Scene_Entity() = default;
  Instance_Scene_Entity(const std::string &n,
                        File_Loc loc,
                        const std::string &material_name,
                        const ProjectionTransform &render_from_instance_anim)
      : name(n),
        loc(loc),
        material_name(material_name),
        render_from_instance_anim(new ProjectionTransform(render_from_instance_anim))
  {
  }

  Instance_Scene_Entity(const std::string &n,
                        File_Loc loc,
                        const std::string &material_name,
                        const ProjectionTransform *render_from_instance)
      : name(n), loc(loc), material_name(material_name), render_from_instance(render_from_instance)
  {
  }

  // For RenderMan
  Instance_Scene_Entity(const std::string &n,
                        File_Loc loc,
                        const std::string &material_name,
                        Mapped_Parameter_Dictionary params,
                        const ProjectionTransform &render_from_instance_anim)
      : name(n),
        loc(loc),
        material_name(material_name),
        parameters(params),
        render_from_instance_anim(new ProjectionTransform(render_from_instance_anim))
  {
  }

  Instance_Scene_Entity(const std::string &n,
                        File_Loc loc,
                        const std::string &material_name,
                        Mapped_Parameter_Dictionary params,
                        const ProjectionTransform *render_from_instance)
      : name(n),
        loc(loc),
        material_name(material_name),
        parameters(params),
        render_from_instance(render_from_instance)
  {
  }

  std::string to_string() const
  {
    std::stringstream ss;
    ss << "[ Instance_Scene_Entity name: " << name;
    ss << " parameters: ";
    std::for_each(parameters.begin(),
                  parameters.end(),
                  [&](std::pair<std::string, Parameter_Dictionary> parameter) {
                    ss << "[ " << parameter.first << ": " << parameter.second.to_string() << " ] ";
                  });
    ss << " loc: " << loc.to_string();
    ss << " render_from_instance_anim: ";  // << render_from_instance_anim;
    ss << " render_from_instance: ";       // << render_from_instance;
    return ss.str();
  }

  std::string name;
  File_Loc loc;
  std::string material_name;
  Mapped_Parameter_Dictionary parameters;
  ProjectionTransform *render_from_instance_anim = nullptr;
  const ProjectionTransform *render_from_instance = nullptr;
};

CCL_NAMESPACE_END

#endif  // SCENE_ENTITIES_H
