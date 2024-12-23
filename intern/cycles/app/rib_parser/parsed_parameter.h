#ifndef PARSED_PARAMETER_H
#define PARSED_PARAMETER_H

#include <string>
#include <string_view>
#include <variant>

#include "error.h"
#include "util/vector.h"

CCL_NAMESPACE_BEGIN
/** @brief Parsed_Parameter.
 * @details
 */

enum class Container_Type { Constant, FaceVarying, Reference, Uniform, Varying, Vertex };

enum class Parameter_Type {
  Boolean,
  Bxdf,
  Color,
  Integer,
  Normal,
  Parameter,
  Pointer,
  Point2,
  Point3,
  Real,
  String,
  Texture,
  Vector2,
  Vector3,
  Unknown
};

using Payload_Type =
    std::variant<vector<uint8_t>, vector<float>, vector<int>, vector<std::string>, vector<void *>>;

class Parsed_Parameter {
 public:
  /// @name Initialization
  ///@{
  Parsed_Parameter() = default;
  Parsed_Parameter(File_Loc loc) : loc(loc) {}
  Parsed_Parameter(Parameter_Type type, std::string name, File_Loc loc)
      : type(type), name(name), loc(loc)
  {
    switch (type) {
      case Parameter_Type::Boolean:
        payload = vector<uint8_t>();
        break;
      case Parameter_Type::Color:
      case Parameter_Type::Vector2:
      case Parameter_Type::Vector3:
      case Parameter_Type::Normal:
      case Parameter_Type::Point2:
      case Parameter_Type::Point3:
      case Parameter_Type::Real:
        payload = vector<float>();
        break;
      case Parameter_Type::Integer:
        payload = vector<int>();
        break;
      case Parameter_Type::Bxdf:
      case Parameter_Type::Parameter:
      case Parameter_Type::String:
      case Parameter_Type::Texture:
        payload = vector<std::string>();
        break;
      case Parameter_Type::Pointer:
        payload = vector<void *>();
        break;
      case Parameter_Type::Unknown:
        break;
    }
  }
  Parsed_Parameter(Parsed_Parameter const &pp);
  ///@}

  /// @name Access
  ///@{
  vector<float> const &floats() const
  {
    return std::get<vector<float>>(payload);
  }
  vector<float> &floats()
  {
    return std::get<vector<float>>(payload);
  }
  vector<int> const &ints() const
  {
    return std::get<vector<int>>(payload);
  }
  vector<int> &ints()
  {
    return std::get<vector<int>>(payload);
  }
  vector<std::string> const &strings() const
  {
    return std::get<vector<std::string>>(payload);
  }
  vector<std::string> &strings()
  {
    return std::get<vector<std::string>>(payload);
  }
  vector<uint8_t> const &bools() const
  {
    return std::get<vector<uint8_t>>(payload);
  }
  vector<uint8_t> &bools()
  {
    return std::get<vector<uint8_t>>(payload);
  }
  vector<void*> const &pointers() const
  {
    return std::get<vector<void*>>(payload);
  }
  vector<void*> &pointers()
  {
    return std::get<vector<void*>>(payload);
  }

  Parameter_Type type = Parameter_Type::Unknown;
  std::string name;
  int elem_per_item = 1;
  Container_Type storage = Container_Type::Constant;
  File_Loc loc;
  Payload_Type payload;
  mutable bool looked_up = false;
  bool may_be_unused = false;
  ///@}
  /// @name Status Report
  ///@{
  bool floats_are(float value) const
  {
    bool result = true;
    for (auto it = floats().begin(); it != floats().end() && result; ++it) {
      result &= (*it == value);
    }

    return result;
  }

  bool has_bools() const
  {
    return std::holds_alternative<vector<uint8_t>>(payload);
  }
  bool has_ints() const
  {
    return std::holds_alternative<vector<int>>(payload);
  }
  bool has_floats() const
  {
    return std::holds_alternative<vector<float>>(payload);
  }
  bool has_strings() const
  {
    return std::holds_alternative<vector<std::string>>(payload);
  }
  bool has_pointerss() const
  {
    return std::holds_alternative<vector<void*>>(payload);
  }
  ///@}
  /// @name Conversion
  ///@{
  std::string to_string() const;
  ///@}
  /// @name Basic operations
  ///@{
  void add_float(float v);
  void add_int(int i);
  void add_string(std::string_view str);
  void add_bool(bool v);
  void add_pointer(void* v);
  ///@}
};  // end of class Parsed_Parameter

using Parsed_Parameter_Vector = vector<Parsed_Parameter *>;

CCL_NAMESPACE_END

#endif  // PARSED_PARAMETER_H
