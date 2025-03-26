#pragma once

#include <vector>

#include "scene/scene.h"
#include "session/buffers.h"
#include "session/session.h"

CCL_NAMESPACE_BEGIN

struct Options {
  unique_ptr<Session> session;
  Scene *scene;
  BufferParams* buffer_params;
  string filepath;
  int width, height;
  SceneParams scene_params;
  SessionParams session_params;
  bool quiet;
  bool show_help, interactive, pause;
  string output_filepath;
  string output_pass;
  string display_type;
  string display_server = "127.0.0.1:14158";
  std::vector<float> crop_window = {0., 1., 0., 1.};
  bool crop_window_set = false;
#ifdef WITH_CYCLES_DISTRIBUTED
  bool reverse_connect = false;
  bool is_distributed = false;
  string directory = "/tmp";
#endif
};

CCL_NAMESPACE_END
