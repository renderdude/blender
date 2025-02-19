#pragma once

#include "scene/scene.h"

#include "app/rib_parser/scene_entities.h"

CCL_NAMESPACE_BEGIN

void export_lights(Scene *scene,
                   Instance_Scene_Entity &inst,
                   Instance_Definition_Scene_Entity *inst_def);

CCL_NAMESPACE_END
