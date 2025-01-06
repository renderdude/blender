#pragma once

#include "scene/attribute.h"
#include "app/rib_parser/parsed_parameter.h"

CCL_NAMESPACE_BEGIN

void apply_primvars(AttributeSet &attributes,
                   const ustring &name,
                   Parsed_Parameter* value,
                   AttributeElement elem,
                   AttributeStandard std);

CCL_NAMESPACE_END
