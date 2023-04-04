/*
 * Copyright 2011-2022 Blender Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include <atomic>
#include <string>
#include <vector>

#include "session/display_driver.h"

#include "util/function.h"
#include "util/unique_ptr.h"

CCL_NAMESPACE_BEGIN

class Display_Item;

struct Display_Info {
  /* Indicates whether texture creation was attempted and succeeded.
   * Used to avoid multiple attempts of texture creation on GPU issues or GPU context
   * misconfiguration. */
  bool creation_attempted = false;
  bool is_created = false;

  /* Is true when new data was written to the PBO, meaning, the texture might need to be resized
   * and new data is to be uploaded to the GPU. */
  bool need_update = false;

  /* Content of the texture is to be filled with zeroes. */
  bool need_clear = true;

  /* Dimensions of the texture in pixels. */
  int width = 0;
  int height = 0;
  int full_width = 0;
  int full_height = 0;

  /* Dimensions of the underlying PBO. */
  int buffer_width = 0;
  int buffer_height = 0;

  std::vector<half4> pixels;
};

using Display_Pair = std::pair<DisplayDriver::Params, Display_Info>;

class TEVDisplayDriver : public DisplayDriver {
 public:
  /* Callbacks for enabling and disabling the OpenGL context. Must be provided to support enabling
   * the context on the Cycles render thread independent of the main thread. */
  TEVDisplayDriver(std::string display_server);
  ~TEVDisplayDriver();

  virtual void clear() override;
  void update();

 protected:
  void connect_to_display_server(const std::string &host);
  void disconnect_from_display_server();

  virtual void next_tile_begin() override;

  virtual bool update_begin(const Params &params, int texture_width, int texture_height) override;
  virtual void update_end() override;

  virtual half4 *map_texture_buffer() override;
  virtual void unmap_texture_buffer() override;

  virtual void draw(const Params &params) override;

  /* Texture which contains pixels of the render result. */
  Display_Info texture_;
  Display_Pair _item;

  std::string _display_server;
  Display_Item *_current_item = nullptr;
};

CCL_NAMESPACE_END
