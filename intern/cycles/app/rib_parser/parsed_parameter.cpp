#include "parsed_parameter.h"
#include <cassert>
#include <sstream>

CCL_NAMESPACE_BEGIN

Parsed_Parameter::Parsed_Parameter(Parsed_Parameter const &pp)
{
  type = pp.type;
  name = pp.name;
  elem_per_item = pp.elem_per_item;
  storage = pp.storage;
  loc = pp.loc;
  if (std::holds_alternative<vector<std::string>>(pp.payload)) {
    payload = vector<std::string>();
    for (auto s : std::get<vector<std::string>>(pp.payload)) {
      std::get<vector<std::string>>(payload).push_back(s);
    }
  }
  else {
    payload = pp.payload;
  }
  looked_up = pp.looked_up;
  may_be_unused = pp.may_be_unused;
}

void Parsed_Parameter::add_bool(bool v)
{
  assert(has_bools());
  bools().push_back(v);
}

void Parsed_Parameter::add_float(float v)
{
  assert(has_floats());
  floats().push_back(v);
}

void Parsed_Parameter::add_int(int i)
{
  assert(has_ints());
  ints().push_back(i);
}

void Parsed_Parameter::add_pointer(void* v)
{
  assert(has_pointerss());
  pointers().push_back(v);
}

void Parsed_Parameter::add_string(std::string_view str)
{
  // assert(has_strings());
  strings().push_back({str.begin(), str.end()});
}

std::string Parsed_Parameter::to_string() const
{
  std::stringstream ss;
  ss << "\"";
  switch (type) {
    case Parameter_Type::Boolean:
      ss << "bool";
      break;
    case Parameter_Type::Color:
    case Parameter_Type::Vector2:
    case Parameter_Type::Vector3:
    case Parameter_Type::Normal:
    case Parameter_Type::Point2:
    case Parameter_Type::Point3:
    case Parameter_Type::Real:
      ss << "float";
      break;
    case Parameter_Type::Integer:
      ss << "int";
      break;
    case Parameter_Type::Bxdf:
    case Parameter_Type::Parameter:
    case Parameter_Type::String:
    case Parameter_Type::Texture:
      ss << "string";
      break;
    case Parameter_Type::Pointer:
      ss << "pointer";
      break;
    case Parameter_Type::Unknown:
      break;
  }
  ss << " " << name << "\" [ ";
  if (std::holds_alternative<vector<uint8_t>>(payload)) {
    for (bool b : bools()) {
      ss << (b ? "true " : "false ");
    }
  }
  else if (std::holds_alternative<vector<float>>(payload)) {
    for (float d : floats()) {
      ss << d << " ";
    }
  }
  else if (std::holds_alternative<vector<int>>(payload)) {
    for (int i : ints()) {
      ss << i << " ";
    }
  }
  else if (std::holds_alternative<vector<std::string>>(payload)) {
    for (const auto &s : strings()) {
      ss << "\"" << s << "\" ";
    }
  }
  else if (std::holds_alternative<vector<void*>>(payload)) {
    for (const auto &s : pointers()) {
      ss << "\"" << s << "\" ";
    }
  }
  ss << "] ";

  return ss.str();
}

CCL_NAMESPACE_END
