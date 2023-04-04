#include "param_dict.h"

#include <sstream>

#include "util/log.h"

CCL_NAMESPACE_BEGIN

template<> struct Parameter_Type_Traits<Parameter_Type::Boolean> {
  static constexpr char typeName[] = "bool";
  static constexpr int nPerItem = 1;
  using Return_Type = uint8_t;
  static bool convert(const uint8_t *v, const File_Loc *loc)
  {
    return *v;
  }
  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.bools;
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::Boolean>::typeName[];

// Bxdf Parameter_Type_Traits Definition
template<> struct Parameter_Type_Traits<Parameter_Type::Bxdf> {
  static constexpr char typeName[] = "bxdf";
  static constexpr int nPerItem = 1;
  using Return_Type = std::string;
  static std::string convert(const std::string *s, const File_Loc *loc)
  {
    return *s;
  }
  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.strings;
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::Bxdf>::typeName[];

// Color Parameter_Type_Traits Definition
template<> struct Parameter_Type_Traits<Parameter_Type::Color> {
  // Parameter_Type::Point3 Type Traits
  using Return_Type = float3;

  static constexpr char typeName[] = "color";

  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.floats;
  }

  static constexpr int nPerItem = 3;

  static float3 convert(const float *v, const File_Loc *loc)
  {
    return make_float3(v[0], v[1], v[2]);
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::Color>::typeName[];

template<> struct Parameter_Type_Traits<Parameter_Type::Real> {
  static constexpr char typeName[] = "float";
  static constexpr int nPerItem = 1;
  using Return_Type = float;
  static float convert(const float *v, const File_Loc *loc)
  {
    return *v;
  }
  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.floats;
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::Real>::typeName[];

template<> struct Parameter_Type_Traits<Parameter_Type::Integer> {
  static constexpr char typeName[] = "int";
  static constexpr int nPerItem = 1;
  using Return_Type = int;
  static int convert(const int *i, const File_Loc *loc)
  {
    return *i;
  }
  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.ints;
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::Integer>::typeName[];

template<> struct Parameter_Type_Traits<Parameter_Type::Point2> {
  static constexpr char typeName[] = "point2";
  static constexpr int nPerItem = 2;
  using Return_Type = float2;
  static float2 convert(const float *v, const File_Loc *loc)
  {
    return make_float2(v[0], v[1]);
  }
  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.floats;
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::Point2>::typeName[];

template<> struct Parameter_Type_Traits<Parameter_Type::Vector2> {
  static constexpr char typeName[] = "vector2";
  static constexpr int nPerItem = 2;
  using Return_Type = float2;
  static float2 convert(const float *v, const File_Loc *loc)
  {
    return make_float2(v[0], v[1]);
  }
  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.floats;
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::Vector2>::typeName[];

// Point3 Parameter_Type_Traits Definition
template<> struct Parameter_Type_Traits<Parameter_Type::Point3> {
  // Parameter_Type::Point3 Type Traits
  using Return_Type = float3;

  static constexpr char typeName[] = "point";

  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.floats;
  }

  static constexpr int nPerItem = 3;

  static float3 convert(const float *v, const File_Loc *loc)
  {
    return make_float3(v[0], v[1], v[2]);
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::Point3>::typeName[];

template<> struct Parameter_Type_Traits<Parameter_Type::Vector3> {
  static constexpr char typeName[] = "vector";
  static constexpr int nPerItem = 3;
  using Return_Type = float3;
  static float3 convert(const float *v, const File_Loc *loc)
  {
    return make_float3(v[0], v[1], v[2]);
  }
  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.floats;
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::Vector3>::typeName[];

template<> struct Parameter_Type_Traits<Parameter_Type::Normal> {
  static constexpr char typeName[] = "normal";
  static constexpr int nPerItem = 3;
  using Return_Type = float3;
  static float3 convert(const float *v, const File_Loc *loc)
  {
    return make_float3(v[0], v[1], v[2]);
  }
  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.floats;
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::Normal>::typeName[];

template<> struct Parameter_Type_Traits<Parameter_Type::String> {
  static constexpr char typeName[] = "string";
  static constexpr int nPerItem = 1;
  using Return_Type = std::string;
  static std::string convert(const std::string *s, const File_Loc *loc)
  {
    return *s;
  }
  static const auto &get_values(const Parsed_Parameter &param)
  {
    return param.strings;
  }
};

constexpr char Parameter_Type_Traits<Parameter_Type::String>::typeName[];

///////////////////////////////////////////////////////////////////////////
// Parameter_Dictionary

Parameter_Dictionary::Parameter_Dictionary(Parsed_Parameter_Vector p) : params(std::move(p))
{
  _num_owned_params = params.size();
  std::reverse(params.begin(), params.end());
  check_parameter_types();
}

Parameter_Dictionary::Parameter_Dictionary(Parsed_Parameter_Vector p0,
                                           const Parsed_Parameter_Vector &params1)
    : params(std::move(p0))
{
  _num_owned_params = params.size();
  std::reverse(params.begin(), params.end());
  params.insert(params.end(), params1.rbegin(), params1.rend());
  check_parameter_types();
}

void Parameter_Dictionary::check_parameter_types()
{
  for (const Parsed_Parameter *p : params) {
    if (p->type == Parameter_Type_Traits<Parameter_Type::Boolean>::typeName) {
      if (p->bools.empty()) {
        std::stringstream ss;
        ss << "\"" << p->name << "\":";
        ss << "non-Boolean values provided for Boolean-valued parameter";
        error_exit(&p->loc, ss.str());
      }
    }
    else if (p->type == Parameter_Type_Traits<Parameter_Type::Boolean>::typeName ||
             p->type == Parameter_Type_Traits<Parameter_Type::Bxdf>::typeName ||
             p->type == Parameter_Type_Traits<Parameter_Type::Color>::typeName ||
             p->type == Parameter_Type_Traits<Parameter_Type::Real>::typeName ||
             p->type == Parameter_Type_Traits<Parameter_Type::Integer>::typeName ||
             p->type == Parameter_Type_Traits<Parameter_Type::Point2>::typeName ||
             p->type == Parameter_Type_Traits<Parameter_Type::Vector2>::typeName ||
             p->type == Parameter_Type_Traits<Parameter_Type::Point3>::typeName ||
             p->type == Parameter_Type_Traits<Parameter_Type::Vector3>::typeName ||
             p->type == Parameter_Type_Traits<Parameter_Type::Normal>::typeName) {
      if ((p->ints.empty() && p->floats.empty()) &&
          (p->storage == Container_Type::Reference && p->strings.empty())) {
        std::stringstream ss;
        ss << "\"" << p->name << "\":";
        ss << "non-numeric values provided for numeric-valued parameter";
        error_exit(&p->loc, ss.str());
      }
    }
    else if (p->type == Parameter_Type_Traits<Parameter_Type::String>::typeName ||
             p->type == "texture") {
      if (p->strings.empty()) {
        std::stringstream ss;
        ss << "\"" << p->name << "\":";
        ss << "non-string values provided for string-valued parameter";
        error_exit(&p->loc, ss.str());
      }
    }
    else {
      std::stringstream ss;
      ss << "\"" << p->type << "\":";
      ss << "unknown parameter type";
      error_exit(&p->loc, ss.str());
    }
  }
}

// Parameter_Dictionary Method Definitions
template<Parameter_Type PT>
typename Parameter_Type_Traits<PT>::Return_Type Parameter_Dictionary::lookup_single(
    const std::string &name, typename Parameter_Type_Traits<PT>::Return_Type defaultValue) const
{
  // Search _params_ for parameter _name_
  using traits = Parameter_Type_Traits<PT>;
  for (const Parsed_Parameter *p : params) {
    if (p->name != name || p->type != traits::typeName)
      continue;
    // Extract parameter values from _p_
    const auto &values = traits::get_values(*p);

    // Issue error if an incorrect number of parameter values were provided
    if (values.empty()) {
      std::stringstream ss;
      ss << "No values provided for parameter ";
      ss << "\"" << name << "\".";
      error_exit(&p->loc, ss.str());
    }
    if (values.size() > traits::nPerItem) {
      std::stringstream ss;
      ss << "Expected " << traits::nPerItem << " values for parameter ";
      ss << "\"" << name << "\".";
      error_exit(&p->loc, ss.str());
    }

    // Return parameter values as _Return_Type_
    p->looked_up = true;
    return traits::convert(values.data(), &p->loc);
  }

  return defaultValue;
}

void Parameter_Dictionary::free_parameters()
{
  for (int i = 0; i < _num_owned_params; ++i)
    delete params[i];
  params.clear();
}

float Parameter_Dictionary::get_one_float(const std::string &name, float def) const
{
  return lookup_single<Parameter_Type::Real>(name, def);
}

int Parameter_Dictionary::get_one_int(const std::string &name, int def) const
{
  return lookup_single<Parameter_Type::Integer>(name, def);
}

bool Parameter_Dictionary::get_one_bool(const std::string &name, bool def) const
{
  return lookup_single<Parameter_Type::Boolean>(name, def);
}

float2 Parameter_Dictionary::get_one_point2(const std::string &name, float2 def) const
{
  return lookup_single<Parameter_Type::Point2>(name, def);
}

float2 Parameter_Dictionary::get_one_vector2(const std::string &name, float2 def) const
{
  return lookup_single<Parameter_Type::Vector2>(name, def);
}

float3 Parameter_Dictionary::get_one_color(const std::string &name, float3 def) const
{
  return lookup_single<Parameter_Type::Color>(name, def);
}

float3 Parameter_Dictionary::get_one_point3(const std::string &name, float3 def) const
{
  return lookup_single<Parameter_Type::Point3>(name, def);
}

float3 Parameter_Dictionary::get_one_vector3(const std::string &name, float3 def) const
{
  return lookup_single<Parameter_Type::Vector3>(name, def);
}

float3 Parameter_Dictionary::get_one_normal(const std::string &name, float3 def) const
{
  return lookup_single<Parameter_Type::Normal>(name, def);
}

std::string Parameter_Dictionary::get_one_string(const std::string &name,
                                                 const std::string &def) const
{
  return lookup_single<Parameter_Type::String>(name, def);
}

Parsed_Parameter const *Parameter_Dictionary::get_parameter(const std::string &name) const
{
  for (const Parsed_Parameter *p : params) {
    if (p->name != name)
      continue;

    return p;
  }

  return nullptr;
}

Parsed_Parameter *Parameter_Dictionary::get_parameter(const std::string &name)
{
  for (Parsed_Parameter *p : params) {
    if (p->name != name)
      continue;

    return p;
  }

  return nullptr;
}

template<typename Return_Type, typename ValuesType, typename C>
static vector<Return_Type> return_array(const ValuesType &values,
                                        const Parsed_Parameter &param,
                                        int nPerItem,
                                        C convert)
{
  if (values.empty()) {
    std::stringstream ss;
    ss << "No values provided for ";
    ss << "\"" << param.name << "\".";
    error_exit(&param.loc, ss.str());
  }
  if (values.size() % nPerItem) {
    std::stringstream ss;
    ss << "Number of values provided for \"" << param.name << "\" ";
    ss << "not a multiple of " << nPerItem;
    error_exit(&param.loc, ss.str());
  }

  param.looked_up = true;
  size_t n = values.size() / nPerItem;
  vector<Return_Type> v(n);
  for (size_t i = 0; i < n; ++i)
    v[i] = convert(&values[nPerItem * i], &param.loc);
  return v;
}

template<typename Return_Type, typename G, typename C>
vector<Return_Type> Parameter_Dictionary::lookup_array(const std::string &name,
                                                       Parameter_Type type,
                                                       const char *typeName,
                                                       int nPerItem,
                                                       G get_values,
                                                       C convert) const
{
  for (const Parsed_Parameter *p : params)
    if (p->name == name && p->type == typeName)
      return return_array<Return_Type>(get_values(*p), *p, nPerItem, convert);

  return {};
}

template<Parameter_Type PT>
vector<typename Parameter_Type_Traits<PT>::Return_Type> Parameter_Dictionary::lookup_array(
    const std::string &name) const
{
  using traits = Parameter_Type_Traits<PT>;
  return lookup_array<typename traits::Return_Type>(
      name, PT, traits::typeName, traits::nPerItem, traits::get_values, traits::convert);
}

vector<float> Parameter_Dictionary::get_float_array(const std::string &name) const
{
  return lookup_array<Parameter_Type::Real>(name);
}

vector<int> Parameter_Dictionary::get_int_array(const std::string &name) const
{
  return lookup_array<Parameter_Type::Integer>(name);
}

vector<uint8_t> Parameter_Dictionary::get_bool_array(const std::string &name) const
{
  return lookup_array<Parameter_Type::Boolean>(name);
}

vector<float2> Parameter_Dictionary::get_point2_array(const std::string &name) const
{
  return lookup_array<Parameter_Type::Point2>(name);
}

vector<float2> Parameter_Dictionary::get_vector2_array(const std::string &name) const
{
  return lookup_array<Parameter_Type::Vector2>(name);
}

vector<float3> Parameter_Dictionary::get_color_array(const std::string &name) const
{
  return lookup_array<Parameter_Type::Color>(name);
}

vector<float3> Parameter_Dictionary::get_point3_array(const std::string &name) const
{
  return lookup_array<Parameter_Type::Point3>(name);
}

vector<float3> Parameter_Dictionary::get_vector3_array(const std::string &name) const
{
  return lookup_array<Parameter_Type::Vector3>(name);
}

vector<float3> Parameter_Dictionary::get_normal_array(const std::string &name) const
{
  return lookup_array<Parameter_Type::Normal>(name);
}

vector<std::string> Parameter_Dictionary::get_string_array(const std::string &name) const
{
  return lookup_array<Parameter_Type::String>(name);
}

std::string Parameter_Dictionary::get_texture(const std::string &name) const
{
  for (const Parsed_Parameter *p : params) {
    if (p->name != name || p->type != "texture")
      continue;

    if (p->strings.empty()) {
      std::stringstream ss;
      ss << "No string values provided for parameter ";
      ss << "\"" << name << "\".";
      error_exit(&p->loc, ss.str());
    }
    if (p->strings.size() > 1) {
      std::stringstream ss;
      ss << "More than one value provided for parameter ";
      ss << "\"" << name << "\".";
      error_exit(&p->loc, ss.str());
    }
    p->looked_up = true;
    return p->strings[0];
  }

  return "";
}

void Parameter_Dictionary::remove(const std::string &name, const char *typeName)
{
  for (auto iter = params.begin(); iter != params.end(); ++iter)
    if ((*iter)->name == name && (*iter)->type == typeName) {
      params.erase(iter);
      return;
    }
}

void Parameter_Dictionary::remove_float(const std::string &name)
{
  remove(name, Parameter_Type_Traits<Parameter_Type::Real>::typeName);
}

void Parameter_Dictionary::remove_int(const std::string &name)
{
  remove(name, Parameter_Type_Traits<Parameter_Type::Integer>::typeName);
}

void Parameter_Dictionary::remove_bool(const std::string &name)
{
  remove(name, Parameter_Type_Traits<Parameter_Type::Boolean>::typeName);
}

void Parameter_Dictionary::remove_point2(const std::string &name)
{
  remove(name, Parameter_Type_Traits<Parameter_Type::Point2>::typeName);
}

void Parameter_Dictionary::remove_vector2(const std::string &name)
{
  remove(name, Parameter_Type_Traits<Parameter_Type::Vector2>::typeName);
}

void Parameter_Dictionary::remove_point3(const std::string &name)
{
  remove(name, Parameter_Type_Traits<Parameter_Type::Point3>::typeName);
}

void Parameter_Dictionary::remove_vector3(const std::string &name)
{
  remove(name, Parameter_Type_Traits<Parameter_Type::Vector3>::typeName);
}

void Parameter_Dictionary::remove_normal(const std::string &name)
{
  remove(name, Parameter_Type_Traits<Parameter_Type::Normal>::typeName);
}

void Parameter_Dictionary::remove_string(const std::string &name)
{
  remove(name, Parameter_Type_Traits<Parameter_Type::String>::typeName);
}

void Parameter_Dictionary::remove_texture(const std::string &name)
{
  remove(name, "texture");
}

void Parameter_Dictionary::remove_spectrum(const std::string &name)
{
  remove(name, "spectrum");
  remove(name, "rgb");
  remove(name, "blackbody");
}

void Parameter_Dictionary::remove_bxdf(const std::string &name)
{
  remove(name, Parameter_Type_Traits<Parameter_Type::Bxdf>::typeName);
}

void Parameter_Dictionary::rename_parameter(const std::string &before, const std::string &after)
{
  for (Parsed_Parameter *p : params)
    if (p->name == before)
      p->name = after;
}

void Parameter_Dictionary::rename_used_textures(const std::map<std::string, std::string> &m)
{
  for (Parsed_Parameter *p : params) {
    if (p->type != "texture")
      continue;

    CHECK_EQ(1, p->strings.size());
    auto iter = m.find(p->strings[0]);
    if (iter != m.end())
      p->strings[0] = iter->second;
  }
}

void Parameter_Dictionary::report_unused() const
{
  // type / name
  vector<std::pair<const std::string *, const std::string *>> seen;

  for (const Parsed_Parameter *p : params) {
    if (p->may_be_unused)
      continue;

    bool haveSeen = std::find_if(seen.begin(),
                                 seen.end(),
                                 [&p](std::pair<const std::string *, const std::string *> p2) {
                                   return *p2.first == p->type && *p2.second == p->name;
                                 }) != seen.end();
    if (p->looked_up) {
      // A parameter may be used when creating an initial Material, say,
      // but then an override from a Shape may shadow it such that its
      // name is already in the seen array.
      if (!haveSeen)
        seen.push_back(std::make_pair(&p->type, &p->name));
    }
    else if (haveSeen) {
      // It's shadowed by another parameter; that's fine.
    }
    else {
      std::stringstream ss;
      ss << "\"" << p->name << "\":";
      ss << " unused parameter.";
      error_exit(&p->loc, ss.str());
    }
  }
}

std::string Parameter_Dictionary::to_string() const
{
  std::string s = "[ Parameter_Dictionary params: ";
  for (const Parsed_Parameter *p : params) {
    s += "[ " + p->to_string() + "] ";
  }
  s += "]";
  return s;
}

const File_Loc *Parameter_Dictionary::loc(const std::string &name) const
{
  for (const Parsed_Parameter *p : params)
    if (p->name == name)
      return &p->loc;
  return nullptr;
}

CCL_NAMESPACE_END
