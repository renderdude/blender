/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2022 NVIDIA Corporation
 * Copyright 2022 Blender Foundation */

#include "app/rib_parser/exporters/node_util.h"
#include "app/rib_parser/parsed_parameter.h"
#include "util/transform.h"
#include <cstdlib>
#include <sys/_types/_size_t.h>

CCL_NAMESPACE_BEGIN

template<typename DstType> DstType convert_to_cycles(Parsed_Parameter const* param)
{
  fprintf(stderr, "Could not convert Parsed_Parameter to Cycles type\n");
  exit(-1);
  return DstType(0);
}

template<> bool convert_to_cycles<bool>(Parsed_Parameter const* param)
{
  return (bool)(param->ints[0]);
}

template<> int convert_to_cycles<int>(Parsed_Parameter const* param)
{
  return param->ints[0];
}

template<> float convert_to_cycles<float>(Parsed_Parameter const* param)
{
  return param->floats[0];
}

template<> float2 convert_to_cycles<float2>(Parsed_Parameter const* param)
{
  return make_float2(param->floats[0], param->floats[1]);
}

template<> float3 convert_to_cycles<float3>(Parsed_Parameter const* param)
{
  return make_float3(param->floats[0], param->floats[1], param->floats[2]);
}

template<> ustring convert_to_cycles<ustring>(Parsed_Parameter const* param)
{
  return ustring(param->strings[0]);
}

Transform convert_to_transform(const vector<float> &values, size_t index = 0)
{
  return make_transform(values[index + 0],
                        values[index + 1],
                        values[index + 2],
                        values[index + 3],
                        values[index + 4],
                        values[index + 5],
                        values[index + 6],
                        values[index + 7],
                        values[index + 8],
                        values[index + 9],
                        values[index + 10],
                        values[index + 11]);
}

template<> Transform convert_to_cycles<Transform>(Parsed_Parameter const* param)
{
  return convert_to_transform(param->floats);
}

template<typename DstType> array<DstType> convert_to_cycles_array(Parsed_Parameter const* param)
{
  fprintf(stderr, "Could not convert Parsed_Parameter to Cycles type\n");
  exit(-1);
}

template<> array<float3> convert_to_cycles_array<float3>(Parsed_Parameter const* param)
{
  array<float3> cyclesArray;
  cyclesArray.reserve(param->floats.size() / 3);
  for (size_t i = 0; i < param->floats.size(); i += 3) {
    cyclesArray.push_back_reserved(
        make_float3(param->floats[i], param->floats[i + 1], param->floats[i + 2]));
  }
  return cyclesArray;
}

template<> array<ustring> convert_to_cycles_array<ustring>(Parsed_Parameter const* param)
{
  array<ustring> cyclesArray;
  cyclesArray.reserve(param->strings.size());
  for (const auto &element : param->strings) {
    cyclesArray.push_back_reserved(ustring(element));
  }
  return cyclesArray;
}

template<> array<Transform> convert_to_cycles_array<Transform>(Parsed_Parameter const* param)
{
  array<Transform> cyclesArray;
  cyclesArray.reserve(param->floats.size() / 16);
  for (size_t i = 0; i < param->floats.size(); i += 16) {
    cyclesArray.push_back_reserved(convert_to_transform(param->floats, i));
  }
  return cyclesArray;
}

void set_node_value(Node *node, const SocketType &socket, Parsed_Parameter const* param)
{
  switch (socket.type) {
    default:
    case SocketType::UNDEFINED:
      fprintf(stderr, "Unexpected conversion: SocketType::UNDEFINED\n");
      break;

    case SocketType::BOOLEAN:
      node->set(socket, convert_to_cycles<bool>(param));
      break;
    case SocketType::FLOAT:
      node->set(socket, convert_to_cycles<float>(param));
      break;
    case SocketType::INT:
      node->set(socket, convert_to_cycles<int>(param));
      break;
    case SocketType::UINT:
      node->set(socket, convert_to_cycles<unsigned int>(param));
      break;
    case SocketType::COLOR:
    case SocketType::VECTOR:
    case SocketType::POINT:
    case SocketType::NORMAL:
      node->set(socket, convert_to_cycles<float3>(param));
      break;
    case SocketType::POINT2:
      node->set(socket, convert_to_cycles<float2>(param));
      break;
    case SocketType::CLOSURE:
      // Handled by node connections
      break;
    case SocketType::STRING:
      node->set(socket, convert_to_cycles<ustring>(param));
      break;
    case SocketType::ENUM:
      // Enum's can accept a string or an int
      if (param->strings.size() > 0) {
        node->set(socket, convert_to_cycles<ustring>(param));
      }
      else {
        node->set(socket, convert_to_cycles<int>(param));
      }
      break;
    case SocketType::TRANSFORM:
      node->set(socket, convert_to_cycles<Transform>(param));
      break;
    case SocketType::NODE:
      // TODO: renderIndex->GetRprim()->cycles_node ?
      fprintf(stderr, "Unimplemented conversion: SocketType::NODE\n");
      break;

    case SocketType::BOOLEAN_ARRAY: {
      auto cyclesArray = convert_to_cycles_array<bool>(param);
      node->set(socket, cyclesArray);
      break;
    }
    case SocketType::FLOAT_ARRAY: {
      auto cyclesArray = convert_to_cycles_array<float>(param);
      node->set(socket, cyclesArray);
      break;
    }
    case SocketType::INT_ARRAY: {
      auto cyclesArray = convert_to_cycles_array<int>(param);
      node->set(socket, cyclesArray);
      break;
    }
    case SocketType::COLOR_ARRAY:
    case SocketType::VECTOR_ARRAY:
    case SocketType::POINT_ARRAY:
    case SocketType::NORMAL_ARRAY: {
      auto cyclesArray = convert_to_cycles_array<float3>(param);
      node->set(socket, cyclesArray);
      break;
    }
    case SocketType::POINT2_ARRAY: {
      auto cyclesArray = convert_to_cycles_array<float2>(param);
      node->set(socket, cyclesArray);
      break;
    }
    case SocketType::STRING_ARRAY: {
      auto cyclesArray = convert_to_cycles_array<ustring>(param);
      node->set(socket, cyclesArray);
      break;
    }
    case SocketType::TRANSFORM_ARRAY: {
      auto cyclesArray = convert_to_cycles_array<Transform>(param);
      node->set(socket, cyclesArray);
      break;
    }
    case SocketType::NODE_ARRAY: {
      // TODO: renderIndex->GetRprim()->cycles_node ?
      fprintf(stderr, "Unimplemented conversion: SocketType::NODE_ARRAY\n");
      break;
    }
  }
}

CCL_NAMESPACE_END
