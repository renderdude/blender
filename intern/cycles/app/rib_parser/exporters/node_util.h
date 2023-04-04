/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2022 NVIDIA Corporation
 * Copyright 2022 Blender Foundation */

#pragma once

#include "app/rib_parser/parsed_parameter.h"
#include "graph/node.h"


CCL_NAMESPACE_BEGIN

void set_node_value(Node *node, const SocketType &socket, Parsed_Parameter const* param);

CCL_NAMESPACE_END
