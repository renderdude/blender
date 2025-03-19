#include "app/rib_parser/parsed_parameter.h"
#include "app/rib_parser/scene_entities.h"
#include "kernel/types.h"
#include "scene/camera.h"
#include "scene/mesh.h"
#include "scene/object.h"
#include "scene/shader_graph.h"
#include "scene/shader_nodes.h"

#include "app/rib_parser/exporters/attribute.h"
#include "app/rib_parser/exporters/geometry.h"
#include "util/log.h"
#include "util/vector.h"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <opensubdiv/vtr/types.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "app/rib_parser/util/mikktspace/mikktspace.hh"

CCL_NAMESPACE_BEGIN

/* Tangent Space */

template<bool is_subd> struct MikkMeshWrapper {
  MikkMeshWrapper(const char *layer_name, const Mesh *mesh, float3 *tangent, float *tangent_sign)
      : mesh(mesh), tangent(tangent), tangent_sign(tangent_sign)
  {
    const AttributeSet &attributes = is_subd ? mesh->subd_attributes : mesh->attributes;

    Attribute *attr_vN = attributes.find(ATTR_STD_VERTEX_NORMAL);
    vertex_normal = attr_vN->data_float3();

    Attribute *attr_uv = attributes.find(ustring(layer_name));
    if (attr_uv != nullptr) {
      texface = attr_uv->data_float2();
    }
  }

  int GetNumFaces()
  {
    if constexpr (is_subd) {
      return mesh->get_num_subd_faces();
    }
    else {
      return mesh->num_triangles();
    }
  }

  int GetNumVerticesOfFace(const int face_num)
  {
    if constexpr (is_subd) {
      return mesh->get_subd_num_corners()[face_num];
    }
    else {
      return 3;
    }
  }

  int CornerIndex(const int face_num, const int vert_num)
  {
    if constexpr (is_subd) {
      const Mesh::SubdFace &face = mesh->get_subd_face(face_num);
      return face.start_corner + vert_num;
    }
    else {
      return face_num * 3 + vert_num;
    }
  }

  int VertexIndex(const int face_num, const int vert_num)
  {
    int corner = CornerIndex(face_num, vert_num);
    if constexpr (is_subd) {
      return mesh->get_subd_face_corners()[corner];
    }
    else {
      return mesh->get_triangles()[corner];
    }
  }

  mikk::float3 GetPosition(const int face_num, const int vert_num)
  {
    const float3 vP = mesh->get_verts()[VertexIndex(face_num, vert_num)];
    return mikk::float3(vP.x, vP.y, vP.z);
  }

  mikk::float3 GetTexCoord(const int face_num, const int vert_num)
  {
    /* TODO: Check whether introducing a template boolean in order to
     * turn this into a constexpr is worth it. */
    if (texface != NULL) {
      const int corner_index = CornerIndex(face_num, vert_num);
      float2 tfuv = texface[corner_index];
      return mikk::float3(tfuv.x, tfuv.y, 1.0f);
    }
    if (orco != NULL) {
      const int vertex_index = VertexIndex(face_num, vert_num);
      const float2 uv = map_to_sphere((orco[vertex_index] + orco_loc) * inv_orco_size);
      return mikk::float3(uv.x, uv.y, 1.0f);
    }
    return mikk::float3(0.0f, 0.0f, 1.0f);
  }

  mikk::float3 GetNormal(const int face_num, const int vert_num)
  {
    float3 vN;
    if (is_subd) {
      const Mesh::SubdFace &face = mesh->get_subd_face(face_num);
      if (face.smooth) {
        const int vertex_index = VertexIndex(face_num, vert_num);
        vN = vertex_normal[vertex_index];
      }
      else {
        vN = face.normal(mesh);
      }
    }
    else {
      if (mesh->get_smooth()[face_num]) {
        const int vertex_index = VertexIndex(face_num, vert_num);
        vN = vertex_normal[vertex_index];
      }
      else {
        const Mesh::Triangle tri = mesh->get_triangle(face_num);
        vN = tri.compute_normal(mesh->get_verts().data());
      }
    }
    return mikk::float3(vN.x, vN.y, vN.z);
  }

  void SetTangentSpace(const int face_num, const int vert_num, mikk::float3 T, bool orientation)
  {
    const int corner_index = CornerIndex(face_num, vert_num);
    tangent[corner_index] = make_float3(T.x, T.y, T.z);
    if (tangent_sign != NULL) {
      tangent_sign[corner_index] = orientation ? 1.0f : -1.0f;
    }
  }

  const Mesh *mesh;
  int num_faces;

  float3 *vertex_normal;
  float2 *texface;
  float3 *orco;
  float3 orco_loc, inv_orco_size;

  float3 *tangent;
  float *tangent_sign;
};

static void mikk_compute_tangents(const char *layer_name, Mesh *mesh, bool need_sign)
{
  /* Create tangent attributes. */
  const bool is_subd = mesh->get_num_subd_faces();
  AttributeSet &attributes = is_subd ? mesh->subd_attributes : mesh->attributes;
  Attribute *attr;
  ustring name = ustring((string(layer_name) + ".tangent").c_str());
  attr = attributes.add(ATTR_STD_UV_TANGENT, name);
  float3 *tangent = attr->data_float3();

  /* Create bitangent sign attribute. */
  float *tangent_sign = nullptr;
  if (need_sign) {
    Attribute *attr_sign;
    ustring name_sign = ustring((string(layer_name) + ".tangent_sign").c_str());

    attr_sign = attributes.add(ATTR_STD_UV_TANGENT_SIGN, name_sign);
    tangent_sign = attr_sign->data_float();
  }

  /* Setup userdata. */
  if (is_subd) {
    MikkMeshWrapper<true> userdata(layer_name, mesh, tangent, tangent_sign);
    /* Compute tangents. */
    mikk::Mikktspace(userdata).genTangSpace();
  }
  else {
    MikkMeshWrapper<false> userdata(layer_name, mesh, tangent, tangent_sign);
    /* Compute tangents. */
    mikk::Mikktspace(userdata).genTangSpace();
  }
}

void RIBCyclesMesh::build_instance_definition(
    Instance_Definition_Scene_Entity const *inst_def)
{
  _shape = inst_def->shapes[0];

  auto shade_pp = _shape.graphics_state.rib_attributes.find("shade");
  if (shade_pp != _shape.graphics_state.rib_attributes.end()) {
    int nfaces = _shape.parameters.get_int_array("nfaces")[0];
    // Check if we're using the entire set of prims for shading
    auto &param = shade_pp->second;
    if (param[0]->elem_per_item < nfaces) {
      _shape = reduce_geometry_by_faceset(_shape, param[0]->ints());
    }
  }

  if (inst_def->shapes.size() > 1) {
    fprintf(stderr,
            "An instance definition, %s, contains more than one shape.\n",
            inst_def->name.c_str());
    fprintf(stderr, "Only using the first shape found.\n");
  }

  initialize(inst_def->name);

  populate();

  _geomTransform = *_shape.render_from_object;
}

void RIBCyclesMesh::build_instance(Instance_Scene_Entity &inst)
{
  array<Node *> usedShaders(1);
  usedShaders[0] = _scene->default_surface;

  for (auto *shader : _scene->shaders) {
    if (!shader->name.compare(inst.material_name)) {
      usedShaders[0] = shader;
      break;
    }
  }

  needs_emission_normalization = false;

  if (inst.parameters.find("Ri") != inst.parameters.end()) {
    if (inst.parameters.at("Ri").get_one_bool("areaNormalize", false)) {
      for (Node *shader : usedShaders) {
        static_cast<Shader *>(shader)->tag_used(_scene);
        for (ShaderNode *snode : static_cast<Shader *>(shader)->graph->nodes) {
          needs_emission_normalization |= snode->has_surface_emission();
        }
      }
    }
  }

  _geom->set_used_shaders(usedShaders);

  // Need to check if we need to update the subd transform
  if (_geom->get_subdivision_type() != Mesh::SUBDIVISION_NONE) {
    const Transform tfm = projection_to_transform(*inst.render_from_instance) *
                          _geom->get_subd_objecttoworld();

    _geom->set_subd_objecttoworld(tfm);
  }

  // Can't populate the primvars until the material is set, otherwise it won't be able to tell
  // that it needs uv's or other attributes
  if (_uv_param) {
    create_uv_map(_uv_param);
  }

  // Must happen after material ID update, so that attribute decisions can be made
  // based on it (e.g. check whether an attribute is actually needed)
  bool rebuild = (_geom->triangles_is_modified()) || (_geom->subd_start_corner_is_modified()) ||
                 (_geom->subd_num_corners_is_modified()) || (_geom->subd_shader_is_modified()) ||
                 (_geom->subd_smooth_is_modified()) || (_geom->subd_ptex_offset_is_modified()) ||
                 (_geom->subd_face_corners_is_modified());

  if (_geom->is_modified() || rebuild) {
    _scene->mutex.lock();
    _geom->tag_update(_scene, rebuild);
    _geom->compute_bounds();
    _scene->mutex.unlock();
  }

  _scene->mutex.lock();
  _instance = _scene->create_node<Object>();
  _scene->mutex.unlock();
  std::string instance_id = initialize_instance(inst);

  // Make sure the first object attribute is the instanceId
  assert(!_instance->attributes.empty() && _instance->attributes.front().name() == instance_id);

  // Default to a single instance with an identity transform
  _instance->attributes.front() = ParamValue(instance_id, -1.0f);

  // Update transform
  const float metersPerUnit = 1.;

  const Transform tfm = transform_scale(make_float3(metersPerUnit)) *
                        projection_to_transform(*(inst.render_from_instance) * _geomTransform);
  _instance->set_tfm(tfm);

  uint visibility = PATH_RAY_ALL_VISIBILITY;
  int rib_vis = inst.parameters.at("visibility").get_one_int("camera", 1);
  if (rib_vis == 0) {
    visibility &= ~PATH_RAY_CAMERA;
  }
  rib_vis = inst.parameters.at("visibility").get_one_int("indirect", 1);
  if (rib_vis == 0) {
    visibility &= ~PATH_RAY_REFLECT;
  }
  rib_vis = inst.parameters.at("visibility").get_one_int("transmission", 1);
  if (rib_vis == 0) {
    visibility &= ~PATH_RAY_SHADOW;
  }

  _instance->set_visibility(visibility);
  if (needs_emission_normalization) {
    normalize_emission(_instance);
  }

  _scene->mutex.lock();
  _instance->tag_update(_scene);
  _instance->compute_bounds(_instance->use_motion());
  _scene->mutex.unlock();

  _bounds.grow(_instance->bounds);
}

void RIBCyclesMesh::initialize(std::string name)
{
  // Create geometry
  _scene->mutex.lock();
  _geom = _scene->create_node<Mesh>();
  _scene->mutex.unlock();
  _geom->name = name;

  vector<int3> tmp_vec;
  triangles.swap(tmp_vec);
}

std::string RIBCyclesMesh::initialize_instance(Instance_Scene_Entity &inst)
{
  _instance->set_geometry(_geom);

  std::string id = inst.parameters.at("identifier").get_one_string("name", "");
  _instance->attributes.emplace_back(id, -1.0f);
  _instance->set_color(make_float3(0.8f, 0.8f, 0.8f));
  _instance->set_random_id(hash_string(id.c_str()));

  return id;
}

void RIBCyclesMesh::populate()
{
  separate_face_varying_normals();

  populate_topology();
  populate_points();

  // Must happen after topology update, so that normals attribute size can be calculated
  populate_normals();

  populate_primvars();
}

void RIBCyclesMesh::separate_face_varying_normals()
{
  Parsed_Parameter const *param = _shape.parameters.get_parameter("N");
  if (param != nullptr && param->storage == Container_Type::FaceVarying) {
    std::unordered_map<int, vector<float3>> normal_map;
    vector<float3> normals = _shape.parameters.get_normal_array("N");
    Parsed_Parameter *points = _shape.parameters.get_parameter("P");

    // Extract the parameters associated with the points
    vector<Parsed_Parameter *> varying_primvars;
    for (auto *pp : _shape.parameters.get_parameter_vector()) {
      if (pp->storage == Container_Type::Varying || pp->storage == Container_Type::Vertex) {
        varying_primvars.push_back(pp);
      }
    }

    vector<int> vertIndx = _shape.parameters.get_int_array("vertices");
    const vector<int> vertCounts = _shape.parameters.get_int_array("nvertices");
    int index_offset = 0;

    for (size_t i = 0; i < vertCounts.size(); i++) {
      for (int j = 0; j < vertCounts[i]; j++) {
        int v0 = vertIndx[index_offset + j];
        float3 N = normals[index_offset + j];
        if (normal_map[v0].empty()) {
          normal_map[v0].push_back(N);
        }
        else {
          bool separate_face = false;
          for (auto &n : normal_map[v0]) {
            if (std::fabs(dot(n, N)) < 0.26) {  // angle > 75 degrees
              separate_face = true;
              for (auto *vp : varying_primvars) {
                if (vp->type == Parameter_Type::Real || vp->type == Parameter_Type::Point2 ||
                    vp->type == Parameter_Type::Point3)
                {
                  for (int i = 0; i < vp->elem_per_item; i++) {
                    vp->floats().push_back(vp->floats()[v0 * vp->elem_per_item + i]);
                  }
                }
                else {
                  std::cerr << "Missed primvar type for " << vp->name << std::endl;
                }
              }
              vertIndx[index_offset + j] = points->floats().size() / 3 - 1;
              break;
            }
          }
          if (!separate_face) {
            normal_map[v0].push_back(N);
          }
        }
      }

      index_offset += vertCounts[i];
    }
    _shape.parameters.get_parameter("vertices")->ints().swap(vertIndx);
  }
}

void RIBCyclesMesh::populate_normals()
{
  _geom->attributes.remove(ATTR_STD_VERTEX_NORMAL);

  // Authored normals should only exist on triangle meshes
  if (_geom->get_subdivision_type() != Mesh::SUBDIVISION_NONE) {
    return;
  }

  Parsed_Parameter const *param = _shape.parameters.get_parameter("N");

  vector<float3> normals;
  Container_Type interpolation;

  // If no normals exist, create vertex normals
  if (param == nullptr) {
    array<float3> &v = _geom->get_verts();
    vector<float3> N(v.size(), make_float3(0, 0, 0));

    for (auto i = 0; i < _geom->num_triangles(); i++) {
      Mesh::Triangle tri = _geom->get_triangle(i);
      float3 n = tri.compute_normal(v.data());
      N[tri.v[0]] += n;
      N[tri.v[1]] += n;
      N[tri.v[2]] += n;
    }

    for (auto i = 0; i < N.size(); i++) {
      auto n = normalize(N[i]);
      N[i] = n;
    }

    normals.swap(N);
    interpolation = Container_Type::Vertex;
  }
  else {
    normals = _shape.parameters.get_normal_array("N");
    interpolation = param->storage;
  }

  float orientation = _shape.graphics_state.reverse_orientation ? -1.0f : 1.0f;

  if (interpolation == Container_Type::Constant) {
    const float3 constantNormal = normals[0];

    float3 *const N = _geom->attributes.add(ATTR_STD_VERTEX_NORMAL)->data_float3();
    for (size_t i = 0; i < _geom->get_verts().size(); ++i) {
      N[i] = orientation * constantNormal;
    }
  }
  else if (interpolation == Container_Type::Uniform) {
    float3 *const N = _geom->attributes.add(ATTR_STD_VERTEX_NORMAL)->data_float3();
    const vector<int> vertIndx = _shape.parameters.get_int_array("vertices");
    const vector<int> vertCounts = _shape.parameters.get_int_array("nvertices");
    int index_offset = 0;

    for (size_t i = 0; i < vertCounts.size(); i++) {
      for (int j = 0; j < vertCounts[i]; j++) {
        int v0 = vertIndx[index_offset + j];
        N[v0] += orientation * normals[index_offset];
      }

      index_offset += vertCounts[i];
    }

    // Now normalize
    for (size_t i = 0; i < _geom->get_verts().size(); ++i) {
      N[i] = normalize(N[i]);
    }
  }
  else if (interpolation == Container_Type::Vertex || interpolation == Container_Type::Varying) {
    float3 *const N = _geom->attributes.add(ATTR_STD_VERTEX_NORMAL)->data_float3();
    for (size_t i = 0; i < _geom->get_verts().size(); ++i) {
      N[i] = orientation * normals[i];
    }
  }
  else if (interpolation == Container_Type::FaceVarying) {
    // Cycles has no standard attribute for face-varying normals, so this is a lossy transformation
    float3 *const N = _geom->attributes.add(ATTR_STD_VERTEX_NORMAL)->data_float3();
    const vector<int> vertIndx = _shape.parameters.get_int_array("vertices");
    const vector<int> vertCounts = _shape.parameters.get_int_array("nvertices");
    int index_offset = 0;

    for (size_t i = 0; i < vertCounts.size(); i++) {
      for (int j = 0; j < vertCounts[i]; j++) {
        int v0 = vertIndx[index_offset + j];
        N[v0] += orientation * normals[index_offset + j];
      }

      index_offset += vertCounts[i];
    }

    // Now normalize
    for (size_t i = 0; i < _geom->get_verts().size(); ++i) {
      N[i] = normalize(N[i]);
    }
  }
}

static std::unordered_map<Container_Type, AttributeElement> interpolations = {
    {std::make_pair(Container_Type::FaceVarying, ATTR_ELEMENT_CORNER)},
    {std::make_pair(Container_Type::Uniform, ATTR_ELEMENT_FACE)},
    {std::make_pair(Container_Type::Vertex, ATTR_ELEMENT_VERTEX)},
    {std::make_pair(Container_Type::Varying, ATTR_ELEMENT_VERTEX)},
    {std::make_pair(Container_Type::Constant, ATTR_ELEMENT_OBJECT)},
};

void RIBCyclesMesh::populate_primvars()
{
  Scene *const scene = (Scene *)_geom->get_owner();

  const bool subdivision = _geom->get_subdivision_type() != Mesh::SUBDIVISION_NONE;
  AttributeSet &attributes = subdivision ? _geom->subd_attributes : _geom->attributes;

  Parsed_Parameter_Vector const &paramv = _shape.parameters.get_parameter_vector();

  for (auto *const param : paramv) {
    // Skip special primvars that are handled separately
    if (param->name == "P" || param->name == "N" || param->name == "nfaces" ||
        param->name == "nvertices" || param->name == "vertices")
    {
      continue;
    }

    const ustring name(param->name);
    AttributeStandard std = ATTR_STD_NONE;
    if (param->name == "st" || param->name == "uv") {
      param->name = "uv";
      _uv_param = param;
      continue;
    }
    if (param->storage == Container_Type::Vertex) {
      if (param->name == "color") {
        std = ATTR_STD_VERTEX_COLOR;
      }
      else if (param->name == "N") {
        std = ATTR_STD_VERTEX_NORMAL;
      }
    }
    /*
    else if (desc.name == HdTokens->displayColor &&
             interpolation.first == HdInterpolationConstant) {
      if (value.IsHolding<VtVec3fArray>() && value.GetArraySize() == 1) {
        const GfVec3f color = value.UncheckedGet<VtVec3fArray>()[0];
        _instances[0]->set_color(make_float3(color[0], color[1], color[2]));
      }
    }
    */

    Parsed_Parameter *result = param;
    // Skip attributes that are not needed
    if ((std != ATTR_STD_NONE && _geom->need_attribute(scene, std)) ||
        _geom->need_attribute(scene, name))
    {

      if (!subdivision) {
        // Adjust attributes for polygons that were triangulated
        if (param->storage == Container_Type::Uniform) {
          result = compute_triangulated_uniform_primvar(param);
          if (!result) {
            continue;
          }
        }
        else if (param->storage == Container_Type::FaceVarying) {
          result = compute_triangulated_face_varying_primvar(param);
          if (!result) {
            continue;
          }
        }
      }

      apply_primvars(attributes, name, result, interpolations[param->storage], std);
    }
  }
}

void RIBCyclesMesh::create_uv_map(Parsed_Parameter *param)
{
  const bool subdivision = _geom->get_subdivision_type() != Mesh::SUBDIVISION_NONE;
  AttributeSet &attributes = subdivision ? _geom->subd_attributes : _geom->attributes;

  AttributeStandard uv_std = ATTR_STD_UV;
  ustring uv_name = ustring("uv");
  AttributeStandard tangent_std = ATTR_STD_UV_TANGENT;
  ustring tangent_name = ustring("uv.tangent");

  /* Denotes whether UV map was requested directly. */
  const bool need_uv = _geom->need_attribute(_scene, uv_name) ||
                       _geom->need_attribute(_scene, uv_std);
  /* Denotes whether tangent was requested directly. */
  const bool need_tangent = _geom->need_attribute(_scene, tangent_name) ||
                            _geom->need_attribute(_scene, tangent_std);

  Parsed_Parameter *result = param;
  if (need_uv || need_tangent) {
    if (!subdivision) {
      // Adjust attributes for polygons that were triangulated
      if (param->storage == Container_Type::Uniform) {
        result = compute_triangulated_uniform_primvar(param);
      }
      else if (param->storage == Container_Type::FaceVarying) {
        result = compute_triangulated_face_varying_primvar(param);
      }
    }

    apply_primvars(attributes, uv_name, result, interpolations[param->storage], uv_std);
  }

  /* UV tangent */
  if (need_tangent) {
    AttributeStandard sign_std = ATTR_STD_UV_TANGENT_SIGN;
    ustring sign_name = ustring("uv.tangent_sign");
    bool need_sign = (_geom->need_attribute(_scene, sign_name) ||
                      _geom->need_attribute(_scene, sign_std));
    Parsed_Parameter const *param = _shape.parameters.get_parameter("N");
    bool need_normals = (param == nullptr || param->storage == Container_Type::FaceVarying) &&
                        subdivision;
    if (need_normals) {
      _geom->add_vertex_normals();
    }
    mikk_compute_tangents("uv", _geom, need_sign);
    if (need_normals) {
      _geom->attributes.remove(ATTR_STD_VERTEX_NORMAL);
    }
  }
}

void RIBCyclesMesh::populate_points()
{
  auto points = _shape.parameters.get_point3_array("P");
  array<float3> P_array;
  P_array = points;

  _geom->set_verts(P_array);
}

void RIBCyclesMesh::populate_topology()
{
  // Clear geometry before populating it again with updated topology
  _geom->clear(true);

  /* Get RIB refinement level
    const HdDisplayStyle displayStyle = GetDisplayStyle(sceneDelegate);
    _topology = HdMeshTopology(GetMeshTopology(sceneDelegate), displayStyle.refineLevel);
  */

  std::string subdivScheme = _shape.parameters.get_one_string("scheme", "");
  if (subdivScheme == "bilinear") {
    _geom->set_subdivision_type(Mesh::SUBDIVISION_LINEAR);
  }
  else if (subdivScheme == "catmull-clark") {
    _geom->set_subdivision_type(Mesh::SUBDIVISION_CATMULL_CLARK);
  }
  else {
    _geom->set_subdivision_type(Mesh::SUBDIVISION_NONE);
  }

  const bool smooth = _shape.parameters.get_one_bool("smooth", false);
  const bool subdivision = _geom->get_subdivision_type() != Mesh::SUBDIVISION_NONE;

  // Initialize lookup table from polygon face to material shader index
  std::vector<int> faceShaders(_shape.parameters.get_one_int("nfaces", 0), 0);

  int shader = 0;
  const vector<int> vertIndx = _shape.parameters.get_int_array("vertices");
  const vector<int> vertCounts = _shape.parameters.get_int_array("nvertices");

  if (!subdivision) {
    compute_triangle_indices(vertIndx, vertCounts, triangles);

    auto points = _shape.parameters.get_point3_array("P");
    _geom->reserve_mesh(points.size(), triangles.size());

    for (size_t i = 0; i < triangles.size(); ++i) {
      const int3 triangle = triangles[i];
      _geom->add_triangle(triangle[0], triangle[1], triangle[2], shader, smooth);
    }
  }
  else {
    /* TODO: subdivision meshes

    PxOsdSubdivTags subdivTags = GetSubdivTags(sceneDelegate);
    _topology.SetSubdivTags(subdivTags);
  */
    size_t numNgons = 0;
    size_t numCorners = 0;
    for (int vertCount : vertCounts) {
      numNgons += (vertCount == 4) ? 0 : 1;
      numCorners += vertCount;
    }

    _geom->reserve_subd_faces(_shape.parameters.get_one_int("nfaces", 0), numCorners);

    // TODO: Handle hole indices
    size_t faceIndex = 0;
    size_t indexOffset = 0;
    for (int vertCount : vertCounts) {
      _geom->add_subd_face(&vertIndx[indexOffset], vertCount, faceShaders[faceIndex], smooth);

      faceIndex++;
      indexOffset += vertCount;
    }

    /*
    const VtIntArray creaseLengths = subdivTags.GetCreaseLengths();
    if (!creaseLengths.empty()) {
      size_t numCreases = 0;
      for (int creaseLength : creaseLengths) {
        numCreases += creaseLength - 1;
      }

      _geom->reserve_subd_creases(numCreases);

      const VtIntArray creaseIndices = subdivTags.GetCreaseIndices();
      const VtFloatArray creaseWeights = subdivTags.GetCreaseWeights();

      indexOffset = 0;
      size_t creaseLengthOffset = 0;
      size_t createWeightOffset = 0;
      for (int creaseLength : creaseLengths) {
        for (int j = 0; j < creaseLength - 1; ++j, ++createWeightOffset) {
          const int v0 = creaseIndices[indexOffset + j];
          const int v1 = creaseIndices[indexOffset + j + 1];

          float weight = creaseWeights.size() == creaseLengths.size() ?
                             creaseWeights[creaseLengthOffset] :
                             creaseWeights[createWeightOffset];

          _geom->add_edge_crease(v0, v1, weight);
        }

        indexOffset += creaseLength;
        creaseLengthOffset++;
      }

      const VtIntArray cornerIndices = subdivTags.GetCornerIndices();
      const VtFloatArray cornerWeights = subdivTags.GetCornerWeights();

      for (size_t i = 0; i < cornerIndices.size(); ++i) {
        _geom->add_vertex_crease(cornerIndices[i], cornerWeights[i]);
      }
    }
  */

    const float metersPerUnit = 1.;

    const Transform tfm = transform_scale(make_float3(metersPerUnit)) *
                          projection_to_transform((*_shape.render_from_object));

    _geom->set_subd_dicing_rate(1.0f);
    _geom->set_subd_max_level(16);
    _geom->set_subd_objecttoworld(tfm);
  }
}

void RIBCyclesMesh::populate_shader_graph([[maybe_unused]] bool initializing) {}

void RIBCyclesMesh::compute_triangle_indices(const vector<int> &vertices,
                                             const vector<int> &nvertices,
                                             vector<int3> &indices)
{
  int index_offset = 0;

  for (size_t i = 0; i < nvertices.size(); i++) {
    for (int j = 0; j < nvertices[i] - 2; j++) {
      int v0 = vertices[index_offset];
      int v1 = vertices[index_offset + j + 1];
      int v2 = vertices[index_offset + j + 2];

      // Reverse orientation for cycles
      indices.push_back(make_int3(v0, v1, v2));
    }

    index_offset += nvertices[i];
  }
}

Parsed_Parameter *RIBCyclesMesh::compute_triangulated_uniform_primvar(
    const Parsed_Parameter *param)
{
  Parsed_Parameter *result = new Parsed_Parameter(*param);
  vector<float> tmp;
  result->floats().swap(tmp);

  const vector<int> nvertices = _shape.parameters.get_int_array("nvertices");

  int index_offset = 0;
  for (size_t i = 0; i < nvertices.size(); i++) {
    for (int j = 0; j < nvertices[i] - 2; j++) {
      result->floats().push_back(param->floats()[index_offset]);
    }
    index_offset++;
  }

  return result;
}

Parsed_Parameter *RIBCyclesMesh::compute_triangulated_face_varying_primvar(
    const Parsed_Parameter *param)
{
  const vector<int> nvertices = _shape.parameters.get_int_array("nvertices");
  auto per_facevarying_size = std::reduce(nvertices.begin(), nvertices.end());
  per_facevarying_size -= nvertices.size() * 2;
  int elem_per_item = param->elem_per_item;

  Parsed_Parameter *result = new Parsed_Parameter(*param);
  vector<float> tmp;
  tmp.reserve(per_facevarying_size * 3 * elem_per_item);
  result->floats().swap(tmp);

  int index_offset = 0;

  for (size_t i = 0; i < nvertices.size(); i++) {
    for (int j = 0; j < nvertices[i] - 2; j++) {
      int ind = index_offset;
      for (int k = 0; k < elem_per_item; ++k) {
        result->floats().push_back(param->floats()[elem_per_item * ind + k]);
      }
      ind = index_offset + j + 1;
      for (int k = 0; k < elem_per_item; ++k) {
        result->floats().push_back(param->floats()[elem_per_item * ind + k]);
      }
      ind = index_offset + j + 2;
      for (int k = 0; k < elem_per_item; ++k) {
        result->floats().push_back(param->floats()[elem_per_item * ind + k]);
      }
    }

    index_offset += nvertices[i];
  }

  return result;
}

Shape_Scene_Entity RIBCyclesMesh::reduce_geometry_by_faceset(Shape_Scene_Entity const &shape,
                                                             vector<int> const &faceset)
{
  Shape_Scene_Entity new_shape(shape);
  std::unordered_map<int, std::unordered_map<Container_Type, vector<int>>> reindexed;
  const vector<int> nvertices = shape.parameters.get_int_array("nvertices");
  const vector<int> vertIndx = shape.parameters.get_int_array("vertices");
  vector<int> new_vertices;
  // Step 1: Extract the vertices that are in the faceset
  int index_offset = 0, face_index = 0;
  for (size_t i = 0; i < nvertices.size(); i++) {
    if (i < faceset[face_index]) {
      index_offset += nvertices[i];
    }
    else {
      for (size_t j = 0; j < nvertices[i]; j++) {
        new_vertices.push_back(vertIndx[index_offset + j]);
      }
      index_offset += nvertices[i];
      face_index++;
      if (face_index >= faceset.size()) {
        break;
      }
    }
  }
  // Step 2: Copy over the points remapping the indices
  vector<int> unique_vertices(new_vertices);
  std::sort(unique_vertices.begin(), unique_vertices.end());
  unique_vertices.erase(unique(unique_vertices.begin(), unique_vertices.end()),
                        unique_vertices.end());
  std::map<int, int> index_map;
  index_offset = 0;
  for (auto uv : unique_vertices) {
    index_map[uv] = index_offset++;
  }

  // Step 3: Subset!
  for (auto *pp : new_shape.parameters.get_parameter_vector()) {
    vector<float> floats;
    vector<int> ints;
    vector<std::string> strings;
    vector<uint8_t> bools;
    switch (pp->storage) {
      case Container_Type::Constant:
      case Container_Type::Reference:
      case Container_Type::Uniform: {
        // Simply copy into the temporaries to make step 4 simpler
        if (pp->has_floats()) {
          floats = pp->floats();
        }
        if (pp->has_ints()) {
          ints = pp->ints();
        }
        if (pp->has_strings()) {
          strings = pp->strings();
        }
        if (pp->has_bools()) {
          bools = pp->bools();
        }
        break;
      }
      case Container_Type::Varying:
      case Container_Type::Vertex: {
        int vert_index = 0;
        for (size_t i = 0; i <= unique_vertices.back(); i++) {
          if (i == unique_vertices[vert_index]) {
            int index = i * pp->elem_per_item;
            if (pp->has_floats() && pp->floats().size() > index) {
              for (size_t j = 0; j < pp->elem_per_item; j++) {
                floats.push_back(pp->floats()[index + j]);
              }
            }
            if (pp->has_ints() && pp->ints().size() > index) {
              for (size_t j = 0; j < pp->elem_per_item; j++) {
                ints.push_back(pp->ints()[index + j]);
              }
            }
            if (pp->has_strings() && pp->strings().size() > index) {
              for (size_t j = 0; j < pp->elem_per_item; j++) {
                strings.push_back(pp->strings()[index + j]);
              }
            }
            if (pp->has_bools() && pp->bools().size() > index) {
              for (size_t j = 0; j < pp->elem_per_item; j++) {
                bools.push_back(pp->bools()[index + j]);
              }
            }
            vert_index++;
          }
        }
        break;
      }
      case Container_Type::FaceVarying: {
        int index_offset = 0, face_index = 0;
        for (size_t i = 0; i < nvertices.size(); i++) {
          int num_elems = nvertices[i] * pp->elem_per_item;
          if (face_index < faceset.size() && i == faceset[face_index]) {
            if (pp->has_floats() && pp->floats().size() > index_offset) {
              for (size_t j = 0; j < num_elems; j++) {
                floats.push_back(pp->floats()[index_offset + j]);
              }
            }
            if (pp->has_ints() && pp->ints().size() > index_offset) {
              for (size_t j = 0; j < num_elems; j++) {
                ints.push_back(pp->ints()[index_offset + j]);
              }
            }
            if (pp->has_strings() && pp->strings().size() > index_offset) {
              for (size_t j = 0; j < num_elems; j++) {
                strings.push_back(pp->strings()[index_offset + j]);
              }
            }
            if (pp->has_bools() && pp->bools().size() > index_offset) {
              for (size_t j = 0; j < num_elems; j++) {
                bools.push_back(pp->bools()[index_offset + j]);
              }
            }
            face_index++;
          }
          index_offset += num_elems;
        }
        break;
      }
    }

    // Step 4: Swap in the subsetted parameters
    if (pp->has_floats()) {
      pp->floats() = floats;
    }
    if (pp->has_ints()) {
      pp->ints() = ints;
    }
    if (pp->has_strings()) {
      pp->strings() = strings;
    }
    if (pp->has_bools()) {
      pp->bools() = bools;
    }
  }

  // Step 5: Fix remaining parameters not fixed in the loop above
  auto *pp = new_shape.parameters.get_parameter("nfaces");
  pp->ints()[0] = faceset.size();
  pp = new_shape.parameters.get_parameter("nvertices");
  vector<int> ints;
  for (int idx = 0; idx < faceset.size(); idx++) {
    ints.push_back(nvertices[faceset[idx]]);
  }
  pp->ints() = ints;
  pp = new_shape.parameters.get_parameter("vertices");
  ints.clear();
  for (int idx = 0; idx < new_vertices.size(); idx++) {
    ints.push_back(index_map[new_vertices[idx]]);
  }
  pp->ints() = ints;

  return new_shape;
}

// This will likely produce weird results if 2+ objects of different sizes share
// the same shader with normalization on
void RIBCyclesMesh::normalize_emission(Object *instance)
{
  Mesh *mesh = static_cast<Mesh *>(instance->get_geometry());
  bool transform_applied = mesh->transform_applied;
  Transform tfm = instance->get_tfm();
  std::unordered_map<int, float> totarea;

  size_t mesh_num_triangles = mesh->num_triangles();
  for (size_t i = 0; i < mesh_num_triangles; i++) {
    int shader_index = mesh->get_shader()[i];

    if (totarea.find(shader_index) == totarea.end()) {
      totarea[shader_index] = 0.f;
    }
    Mesh::Triangle t = mesh->get_triangle(i);
    float3 p1 = mesh->get_verts()[t.v[0]];
    float3 p2 = mesh->get_verts()[t.v[1]];
    float3 p3 = mesh->get_verts()[t.v[2]];

    if (!transform_applied) {
      p1 = transform_point(&tfm, p1);
      p2 = transform_point(&tfm, p2);
      p3 = transform_point(&tfm, p3);
    }

    totarea[shader_index] += triangle_area(p1, p2, p3);
  }

  for (auto &[key, value] : totarea) {
    Shader *shader = static_cast<Shader *>(mesh->get_used_shaders()[key]);
    for (auto *node : shader->graph->nodes) {
      if (node->is_a(PrincipledBsdfNode::node_type)) {
        PrincipledBsdfNode *bsdf = (PrincipledBsdfNode *)node;
        float strength = bsdf->get_emission_strength();
        bsdf->set_emission_strength(strength / value);
      }
      else if (node->is_a(EmissionNode::node_type)) {
        EmissionNode *bsdf = (EmissionNode *)node;
        float strength = bsdf->get_strength();
        bsdf->set_strength(strength / value);
      }
    }
  }
}
CCL_NAMESPACE_END
