#include "ri_api_distributed.h"

CCL_NAMESPACE_BEGIN

void Ri_Distributed::Pattern(const std::string &name,
             const std::string &handle,
             Parsed_Parameter_Vector params,
             File_Loc loc) {
	Ri::Pattern(name, handle, params, loc);
}

CCL_NAMESPACE_END
