#include "app/rib_parser/error.h"
#include "app/rib_parser/parsed_parameter.h"
#include "app/rib_parser/scene_entities.h"
#include "scene/camera.h"
#include "scene/mesh.h"
#include "scene/object.h"
#include "scene/shader_graph.h"
#include "scene/shader_nodes.h"

#include "app/rib_parser/exporters/attribute.h"
#include "app/rib_parser/exporters/geometry.h"
#include "util/vector.h"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <opensubdiv/vtr/types.h>
#include <string>
#include <unordered_map>

#include "app/rib_parser/util/mikktspace/mikktspace.hh"

CCL_NAMESPACE_BEGIN

/* Tangent Space */

template<bool is_subd> struct MikkMeshWrapper {
  MikkMeshWrapper(const char *layer_name, const Mesh *mesh, float3 *tangent, float *tangent_sign)
      : mesh(mesh), texface(NULL), orco(NULL), tangent(tangent), tangent_sign(tangent_sign)
  {
    const AttributeSet &attributes = is_subd ? mesh->subd_attributes : mesh->attributes;

    Attribute *attr_vN = attributes.find(ATTR_STD_VERTEX_NORMAL);
    vertex_normal = attr_vN->data_float3();

    Attribute *attr_uv = attributes.find(ustring(layer_name));
    if (attr_uv != NULL) {
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
    else if (orco != NULL) {
      const int vertex_index = VertexIndex(face_num, vert_num);
      const float2 uv = map_to_sphere((orco[vertex_index] + orco_loc) * inv_orco_size);
      return mikk::float3(uv.x, uv.y, 1.0f);
    }
    else {
      return mikk::float3(0.0f, 0.0f, 1.0f);
    }
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
        vN = tri.compute_normal(&mesh->get_verts()[0]);
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
  float *tangent_sign = NULL;
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

void RIBCyclesMesh::export_geometry()
{
  _shape = _inst_def->shapes[0];

  auto shade_pp = _shape.graphics_state.rib_attributes.find("shade");
  if (shade_pp != _shape.graphics_state.rib_attributes.end()) {
    int nfaces = _shape.parameters.get_int_array("nfaces")[0];
    // Check if we're using the entire set of prims for shading
    auto &param = shade_pp->second;
    if (param[0]->elem_per_item < nfaces)
      _shape = reduce_geometry_by_faceset(_shape, param[0]->ints);
  }

  if (_inst_def->shapes.size() > 1) {
    fprintf(stderr,
            "An instance definition, %s, contains more than one shape.\n",
            _inst_def->name.c_str());
    fprintf(stderr, "Only using the first shape found.\n");
  }

  _instances.resize(_inst_v.size());

  for (size_t i = 0; i < _inst_v.size(); ++i) {
    std::string material_id = _inst_v[i].material_name;
    if (_instanced_geom.find(material_id) != _instanced_geom.end())
      _geom = _instanced_geom[material_id];
    else {
      initialize(material_id);

      array<Node *> usedShaders(1);
      usedShaders[0] = _scene->default_surface;

      for (auto shader : _scene->shaders) {
        if (!shader->name.compare(material_id)) {
          usedShaders[0] = shader;
          break;
        }
      }

      for (Node *shader : usedShaders) {
        static_cast<Shader *>(shader)->tag_used(_scene);
      }

      _geom->set_used_shaders(usedShaders);

      // Must happen after material ID update, so that attribute decisions can be made
      // based on it (e.g. check whether an attribute is actually needed)
      bool rebuild = false;
      populate(rebuild);

      if (_geom->is_modified() || rebuild) {
        _geom->tag_update(_scene, rebuild);
        _geom->compute_bounds();
      }
    }

    _instances[i] = _scene->create_node<Object>();
    initialize_instance(static_cast<int>(i));

    if (i == 0) {
      std::string instance_id = _inst_v[0].parameters.at("identifier").get_one_string("name", "");
      // Make sure the first object attribute is the instanceId
      assert(_instances[0]->attributes.size() >= 1 &&
             _instances[0]->attributes.front().name() == instance_id);

      if (_inst_v.size() > 1) {
        _instances[0]->attributes.front() = ParamValue(instance_id, +0.0f);
      }
      else {
        // Default to a single instance with an identity transform
        _instances[0]->attributes.front() = ParamValue(instance_id, -1.0f);
      }
    }
    // Update transform
    const float metersPerUnit = 1.;

    const Transform tfm = transform_scale(make_float3(metersPerUnit)) *
                          projection_to_transform((*_inst_v[i].render_from_instance) *
                                                  (*_shape.render_from_object));
    _instances[i]->set_tfm(tfm);

    /* Not sure where to pull visibility from a RIB file
    if (HdChangeTracker::IsVisibilityDirty(*dirtyBits, id)) {
      for (Object *instance : _instances) {
        instance->set_visibility(Base::IsVisible() ? ~0 : 0);
      }
    }
  */
  }
  for (Object *instance : _instances) {
    instance->tag_update(_scene);
    instance->compute_bounds(instance->use_motion());
    _bounds.grow(instance->bounds);
  }
}

void RIBCyclesMesh::initialize(std::string material_name)
{
  // Create geometry
  _geom = _scene->create_node<Mesh>();
  _geom->name = _inst_def->name;

  vector<int3> tmp_vec;
  triangles.swap(tmp_vec);

  _instanced_geom[material_name] = _geom;
}

void RIBCyclesMesh::initialize_instance(int index)
{
  Object *instance = _instances[index];
  instance->set_geometry(_geom);

  std::string id = _inst_v[index].parameters.at("identifier").get_one_string("name", "");
  instance->attributes.emplace_back(id,
                                    _instances.size() == 1 ? -1.0f : static_cast<float>(index));
  instance->set_color(make_float3(0.8f, 0.8f, 0.8f));
  instance->set_random_id(hash_uint2(hash_string(_geom->name.c_str()), index));
}

void RIBCyclesMesh::populate(bool &rebuild)
{
  separate_face_varying_normals();

  populate_topology();
  populate_points();

  // Must happen after topology update, so that normals attribute size can be calculated
  populate_normals();

  // Must happen after topology update, so that appropriate attribute set can be selected
  populate_primvars();

  rebuild = (_geom->triangles_is_modified()) || (_geom->subd_start_corner_is_modified()) ||
            (_geom->subd_num_corners_is_modified()) || (_geom->subd_shader_is_modified()) ||
            (_geom->subd_smooth_is_modified()) || (_geom->subd_ptex_offset_is_modified()) ||
            (_geom->subd_face_corners_is_modified());
}

void RIBCyclesMesh::separate_face_varying_normals()
{
  Parsed_Parameter const *param = _shape.parameters.get_parameter("N");
  if (param != nullptr && param->storage == Container_Type::FaceVarying) {
    std::unordered_map<int, vector<float3>> normal_map;
    vector<float3> normals = _shape.parameters.get_normal_array("N");
    Parsed_Parameter* points = _shape.parameters.get_parameter("P");

    // Extract the parameters associated with the points
    vector<Parsed_Parameter*> varying_primvars;
    for (auto pp: _shape.parameters.get_parameter_vector())
      if (pp->storage == Container_Type::Varying || pp->storage == Container_Type::Vertex)
        varying_primvars.push_back(pp);

    vector<int> vertIndx = _shape.parameters.get_int_array("vertices");
    const vector<int> vertCounts = _shape.parameters.get_int_array("nvertices");
    int index_offset = 0;

    for (size_t i = 0; i < vertCounts.size(); i++) {
      for (int j = 0; j < vertCounts[i]; j++) {
        int v0 = vertIndx[index_offset + j];
        float3 N = normals[index_offset + j];
        if (normal_map[v0].size() == 0)
          normal_map[v0].push_back(N);
        else {
          bool seperate_face = false;
          for (auto& n: normal_map[v0]) {
            if (std::fabs(dot(n , N)) < 0.26) { // angle > 75 degrees
              seperate_face = true;
              for (auto vp: varying_primvars) {
                if (vp->type == "float" || vp->type == "point") {
                  for (int i = 0; i < vp->elem_per_item; i++)
                    vp->floats.push_back(vp->floats[v0 * vp->elem_per_item + i]);
                }
                else
                  std::cerr << "Missed primvar type: " << vp->type << std::endl;
              }
              vertIndx[index_offset + j] = points->floats.size()/3 - 1;
              break;
            }
          }
          if (!seperate_face) normal_map[v0].push_back(N);
        }
      }

      index_offset += vertCounts[i];
    }
    _shape.parameters.get_parameter("vertices")->ints.swap(vertIndx);
  }
}

void RIBCyclesMesh::populate_normals()
{
  _geom->attributes.remove(ATTR_STD_FACE_NORMAL);
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
    float3 *const N = _geom->attributes.add(ATTR_STD_FACE_NORMAL)->data_float3();
    for (size_t i = 0; i < _geom->num_triangles(); ++i) {
      N[i] = orientation * normals[i / 2];
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
#if 1
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
    for (size_t i = 0; i < _geom->get_verts().size(); ++i)
      N[i] = normalize(N[i]);
#else
    float3 *const N = _geom->attributes.add(ATTR_STD_FACE_NORMAL)->data_float3();
    int index_offset = 0;
    size_t tri_index = 0;
    const vector<int> vertCounts = _shape.parameters.get_int_array("nvertices");
    for (size_t i = 0; i < vertCounts.size(); i++) {
      for (int j = 0; j < vertCounts[i] - 2; j++) {
        int v0 = index_offset;
        int v1 = index_offset + j + 1;
        int v2 = index_offset + j + 2;

        float3 average_normal = normals[v0] + normals[v1] + normals[v2];
        N[tri_index++] = orientation * normalize(average_normal);
      }

      index_offset += vertCounts[i];
    }
#endif
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

  for (const auto param : paramv) {
    // Skip special primvars that are handled separately
    if (param->name == "P" || param->name == "N" || param->name == "nfaces" ||
        param->name == "nvertices" || param->name == "vertices") {
      continue;
    }

    const ustring name(param->name);
    AttributeStandard std = ATTR_STD_NONE;
    if (param->name == "st" || param->name == "uv") {
      param->name = "uv";
      create_uv_map(param);
      continue;
    }
    else if (param->storage == Container_Type::Vertex) {
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
        _geom->need_attribute(scene, name)) {

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
    mikk_compute_tangents("uv", _geom, need_sign);
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

    size_t numNgons = 0;
    size_t numCorners = 0;
    for (int vertCount : vertCounts) {
      numNgons += (vertCount == 4) ? 0 : 1;
      numCorners += vertCount;
    }

    _geom->reserve_subd_faces(_topology.GetNumFaces(), numNgons, numCorners);

    // TODO: Handle hole indices
    size_t faceIndex = 0;
    size_t indexOffset = 0;
    for (int vertCount : vertCounts) {
      _geom->add_subd_face(&vertIndx[indexOffset], vertCount, faceShaders[faceIndex], smooth);

      faceIndex++;
      indexOffset += vertCount;
    }

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

    _geom->set_subd_dicing_rate(1.0f);
    _geom->set_subd_max_level(_topology.GetRefineLevel());
    _geom->set_subd_objecttoworld(_instances[0]->get_tfm());
    */
  }
}

void RIBCyclesMesh::populate_shader_graph(bool initializing)
{
}

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
  result->floats.swap(tmp);

  const vector<int> nvertices = _shape.parameters.get_int_array("nvertices");

  int index_offset = 0;
  for (size_t i = 0; i < nvertices.size(); i++) {
    for (int j = 0; j < nvertices[i] - 2; j++) {
      result->floats.push_back(param->floats[index_offset]);
    }
    index_offset++;
  }

  return result;
}

Parsed_Parameter *RIBCyclesMesh::compute_triangulated_face_varying_primvar(
    const Parsed_Parameter *param)
{
  Parsed_Parameter *result = new Parsed_Parameter(*param);
  vector<float> tmp;
  result->floats.swap(tmp);

  const vector<int> nvertices = _shape.parameters.get_int_array("nvertices");

  int index_offset = 0;
  int elem_per_item = param->elem_per_item;

  for (size_t i = 0; i < nvertices.size(); i++) {
    for (int j = 0; j < nvertices[i] - 2; j++) {
      int ind = index_offset;
      for (int k = 0; k < elem_per_item; ++k)
        result->floats.push_back(param->floats[elem_per_item * ind + k]);
      ind = index_offset + j + 1;
      for (int k = 0; k < elem_per_item; ++k)
        result->floats.push_back(param->floats[elem_per_item * ind + k]);
      ind = index_offset + j + 2;
      for (int k = 0; k < elem_per_item; ++k)
        result->floats.push_back(param->floats[elem_per_item * ind + k]);
    }

    index_offset += nvertices[i];
  }

  return result;
}

Shape_Scene_Entity RIBCyclesMesh::reduce_geometry_by_faceset(Shape_Scene_Entity const &shape,
                                                             vector<int> &faceset)
{
  Shape_Scene_Entity new_shape(shape);
  std::unordered_map<int, std::unordered_map<Container_Type, vector<int>>> reindexed;
  const vector<int> nvertices = shape.parameters.get_int_array("nvertices");
  const vector<int> vertIndx = shape.parameters.get_int_array("vertices");
  vector<int> new_vertices;
  // Step 1: Extract the vertices that are in the faceset
  int index_offset = 0, face_index = 0;
  for (size_t i = 0; i < nvertices.size(); i++) {
    if (i < faceset[face_index])
      index_offset += nvertices[i];
    else {
      for (size_t j = 0; j < nvertices[i]; j++)
        new_vertices.push_back(vertIndx[index_offset + j]);
      index_offset += nvertices[i];
      face_index++;
      if (face_index >= faceset.size())
        break;
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
  for (auto pp : new_shape.parameters.get_parameter_vector()) {
    vector<float> floats;
    vector<int> ints;
    vector<std::string> strings;
    vector<uint8_t> bools;
    switch (pp->storage) {
      case Container_Type::Constant:
      case Container_Type::Reference:
      case Container_Type::Uniform: {
        // Simply copy into the temporaries to make step 4 simpler
        floats = pp->floats;
        ints = pp->ints;
        strings = pp->strings;
        bools = pp->bools;
        break;
      }
      case Container_Type::Varying:
      case Container_Type::Vertex: {
        int vert_index = 0;
        for (size_t i = 0; i <= unique_vertices.back(); i++) {
          if (i == unique_vertices[vert_index]) {
            int index = i * pp->elem_per_item;
            if (pp->floats.size() > index)
              for (size_t j = 0; j < pp->elem_per_item; j++)
                floats.push_back(pp->floats[index + j]);
            if (pp->ints.size() > index)
              for (size_t j = 0; j < pp->elem_per_item; j++)
                ints.push_back(pp->ints[index + j]);
            if (pp->strings.size() > index)
              for (size_t j = 0; j < pp->elem_per_item; j++)
                strings.push_back(pp->strings[index + j]);
            if (pp->bools.size() > index)
              for (size_t j = 0; j < pp->elem_per_item; j++)
                bools.push_back(pp->bools[index + j]);
            vert_index++;
          }
        }
        break;
      }
      case Container_Type::FaceVarying: {
        int index_offset = 0, face_index = 0;
        for (size_t i = 0; i < nvertices.size(); i++) {
          int num_elems = nvertices[i] * pp->elem_per_item;
          if (i == faceset[face_index]) {
            if (pp->floats.size() > index_offset)
              for (size_t j = 0; j < num_elems; j++)
                floats.push_back(pp->floats[index_offset + j]);
            if (pp->ints.size() > index_offset)
              for (size_t j = 0; j < num_elems; j++)
                ints.push_back(pp->ints[index_offset + j]);
            if (pp->strings.size() > index_offset)
              for (size_t j = 0; j < num_elems; j++)
                strings.push_back(pp->strings[index_offset + j]);
            if (pp->bools.size() > index_offset)
              for (size_t j = 0; j < num_elems; j++)
                bools.push_back(pp->bools[index_offset + j]);
            face_index++;
          }
          index_offset += num_elems;
        }
        break;
      }
    }

    // Step 4: Swap in the subsetted parameters
    pp->floats = floats;
    pp->ints = ints;
    pp->strings = strings;
    pp->bools = bools;
  }

  // Step 5: Fix remaining parameters not fixed in the loop above
  auto pp = new_shape.parameters.get_parameter("nfaces");
  pp->ints[0] = faceset.size();
  pp = new_shape.parameters.get_parameter("nvertices");
  vector<int> ints;
  for (int idx = 0; idx < faceset.size(); idx++)
    ints.push_back(nvertices[faceset[idx]]);
  pp->ints = ints;
  pp = new_shape.parameters.get_parameter("vertices");
  ints.clear();
  for (int idx = 0; idx < new_vertices.size(); idx++)
    ints.push_back(index_map[new_vertices[idx]]);
  pp->ints = ints;

  return new_shape;
}

CCL_NAMESPACE_END
