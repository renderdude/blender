#include "curves.h"
#include "app/rib_parser/exporters/attribute.h"
#include "app/rib_parser/parsed_parameter.h"
#include "scene/hair.h"
#include "scene/object.h"

CCL_NAMESPACE_BEGIN

void RIBCyclesCurves::export_curves()
{
  _shape = _inst_def->shapes[0];

  if (_inst_def->shapes.size() > 1) {
    fprintf(stderr,
            "An instance definition, %s, contains more than one shape.\n",
            _inst_def->name.c_str());
    fprintf(stderr, "Only using the first shape found.\n");
  }

  _instances.resize(_inst_v.size());

  for (size_t i = 0; i < _inst_v.size(); ++i) {
    std::string material_id = _inst_v[i].material_name;
    if (_instanced_geom.find(material_id) != _instanced_geom.end()) {
      _geom = _instanced_geom[material_id];
    }
    else {
      initialize(material_id);

      array<Node *> usedShaders(1);
      usedShaders[0] = _scene->default_surface;

      for (auto *shader : _scene->shaders) {
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

void RIBCyclesCurves::initialize(std::string material_name)
{
  // Create geometry
  _geom = _scene->create_node<Hair>();
  _geom->name = _inst_def->name;

  _instanced_geom[material_name] = _geom;
}

void RIBCyclesCurves::initialize_instance(int index)
{
  Object *instance = _instances[index];
  instance->set_geometry(_geom);

  std::string id = _inst_v[index].parameters.at("identifier").get_one_string("name", "");
  instance->attributes.emplace_back(id,
                                    _instances.size() == 1 ? -1.0f : static_cast<float>(index));
  instance->set_color(make_float3(0.8f, 0.8f, 0.8f));
  instance->set_random_id(hash_uint2(hash_string(_geom->name.c_str()), index));
}

void RIBCyclesCurves::populate(bool &rebuild)
{
  populate_topology();
  populate_points();
  populate_widths();
  populate_primvars();

  rebuild = (_geom->curve_keys_is_modified()) || (_geom->curve_radius_is_modified());
}

void RIBCyclesCurves::populate_widths()
{
  auto width_param = _shape.parameters.get_parameter("width");
  auto widths = _shape.parameters.get_float_array("width");
  array<float> radiusDataCycles;
  radiusDataCycles.reserve(widths.size());

  if (width_param->storage == Container_Type::Constant) {
    const float constantRadius = widths[0] * 0.5f;

    for (size_t i = 0; i < _geom->num_keys(); ++i) {
      radiusDataCycles.push_back_reserved(constantRadius);
    }
  }
  else if (width_param->storage == Container_Type::Vertex) {
    CHECK(widths.size() == _geom->num_keys());

    for (size_t i = 0; i < _geom->num_keys(); ++i) {
      radiusDataCycles.push_back_reserved(widths[i] * 0.5f);
    }
  }

  _geom->set_curve_radius(radiusDataCycles);
}

static std::unordered_map<Container_Type, AttributeElement> interpolations = {
    {std::make_pair(Container_Type::Uniform, ATTR_ELEMENT_CURVE)},
    {std::make_pair(Container_Type::Vertex, ATTR_ELEMENT_CURVE_KEY)},
    {std::make_pair(Container_Type::Varying, ATTR_ELEMENT_CURVE_KEY)},
    {std::make_pair(Container_Type::Constant, ATTR_ELEMENT_OBJECT)},
};

void RIBCyclesCurves::populate_primvars()
{
  Scene *const scene = (Scene *)_geom->get_owner();

  AttributeSet &attributes = _geom->attributes;

  Parsed_Parameter_Vector const &paramv = _shape.parameters.get_parameter_vector();

  for (const auto param : paramv) {
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
      std = ATTR_STD_UV;
    }
    else if (param->name == "color" && param->storage == Container_Type::Constant) {
      auto color = param->floats();
      _instances[0]->set_color(make_float3(color[0], color[1], color[2]));
    }

    Parsed_Parameter *result = param;
    // Skip attributes that are not needed
    if ((std != ATTR_STD_NONE && _geom->need_attribute(scene, std)) ||
        _geom->need_attribute(scene, name))
    {
      apply_primvars(attributes, name, result, interpolations[param->storage], std);
    }
  }
}

void RIBCyclesCurves::populate_points()
{
  auto points = _shape.parameters.get_point3_array("P");

  array<float3> pointsDataCycles;
  pointsDataCycles.reserve(points.size());

  for (const auto &point : points) {
    pointsDataCycles.push_back_reserved(point);
  }

  _geom->set_curve_keys(pointsDataCycles);
}

void RIBCyclesCurves::populate_topology()
{
  // Clear geometry before populating it again with updated topology
  _geom->clear(true);

  std::string type = _shape.parameters.get_one_string("type", "");
  auto nvertices = _shape.parameters.get_int_array("nvertices");
  int count = 0;
  for (int nv : nvertices) {
    count += nv;
  }

  _geom->reserve_curves(nvertices.size(), count);

  for (int curve = 0, key = 0; curve < nvertices.size(); ++curve) {
    // Always reference shader at index zero, which is the primitive material
    _geom->add_curve(key, 0);

    key += nvertices[curve];
  }
}

CCL_NAMESPACE_END
