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

#include "util/texture.h"
#include <atomic>
#include <chrono>
#include <iostream>
#include <list>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <utility>

#ifdef IS_WINDOWS
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <Ws2tcpip.h>
#  include <winsock2.h>
#  undef NOMINMAX
using socket_t = SOCKET;
#else
using socket_t = int;
#  include <arpa/inet.h>
#  include <errno.h>
#  include <netdb.h>
#  include <netinet/in.h>
#  include <signal.h>
#  include <sys/socket.h>
#  include <sys/time.h>
#  include <unistd.h>
#  define SOCKET_ERROR (-1)
#  define INVALID_SOCKET (-1)
#endif

#include "app/tev/display_driver.h"

#include "util/log.h"
#include "util/murmurhash.h"
#include "util/span.h"
#include "util/string.h"

CCL_NAMESPACE_BEGIN

/* --------------------------------------------------------------------
 * Communication
 */

enum SocketError : int {
#ifdef IS_WINDOWS
  Again = EAGAIN,
  ConnRefused = WSAECONNREFUSED,
  WouldBlock = WSAEWOULDBLOCK,
#else
  Again = EAGAIN,
  ConnRefused = ECONNREFUSED,
  WouldBlock = EWOULDBLOCK,
#endif
};

static int close_socket(socket_t socket)
{
#ifdef IS_WINDOWS
  return closesocket(socket);
#else
  return close(socket);
#endif
}

static std::atomic<int> num_active_channels{0};

class IPC_Channel {
 public:
  IPC_Channel(const std::string &host);
  ~IPC_Channel();

  IPC_Channel(const IPC_Channel &) = delete;
  IPC_Channel &operator=(const IPC_Channel &) = delete;

  bool send(p_std::span<const uint8_t> message);

  bool connected() const
  {
    return socketFd != INVALID_SOCKET;
  }

 private:
  void connect();
  void disconnect();

  int numFailures = 0;
  std::string address, port;
  socket_t socketFd = INVALID_SOCKET;
};

IPC_Channel::IPC_Channel(const std::string &hostname)
{
  if (num_active_channels++ == 0) {
#ifdef IS_WINDOWS
    WSADATA wsaData;
    int err = WSAStartup(MAKEWORD(2, 2), &wsaData);
    if (err != NO_ERROR)
      LOG_FATAL("Unable to initialize WinSock: %s", error_string(err));
#else
    // We don't care about getting a SIGPIPE if the display server goes
    // away...
    signal(SIGPIPE, SIG_IGN);
#endif
  }

  size_t split = hostname.find_last_of(':');
  if (split == std::string::npos) {
    std::cerr << "Expected \"host:port\" for display server address. Given \"" << hostname << "\"."
              << std::endl;
  }
  else {
    address = std::string(hostname.begin(), hostname.begin() + split);
    port = std::string(hostname.begin() + split + 1, hostname.end());

    connect();
  }
}

void IPC_Channel::connect()
{
  CHECK_EQ(socketFd, INVALID_SOCKET);

  LOG(INFO) << "Trying to connect to display server";

  struct addrinfo hints = {}, *addrinfo;
  hints.ai_family = PF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;
  int err = getaddrinfo(address.c_str(), port.c_str(), &hints, &addrinfo);
  if (err)
    LOG(ERROR) << gai_strerror(err);

  socketFd = INVALID_SOCKET;
  for (struct addrinfo *ptr = addrinfo; ptr; ptr = ptr->ai_next) {
    socketFd = socket(ptr->ai_family, ptr->ai_socktype, ptr->ai_protocol);
    if (socketFd == INVALID_SOCKET) {
      LOG(INFO) << "socket() failed";
      continue;
    }

#ifdef IS_LINUX
    struct timeval timeout;
    timeout.tv_sec = 3;
    timeout.tv_usec = 0;
    if (setsockopt(socketFd, SOL_SOCKET, SO_SNDTIMEO, &timeout, sizeof(timeout)) == SOCKET_ERROR) {
      LOG(INFO) << "setsockopt() failed";
    }
#endif  // IS_LINUX

    if (::connect(socketFd, ptr->ai_addr, ptr->ai_addrlen) == SOCKET_ERROR) {
#ifdef IS_WINDOWS
      int err = WSAGetLastError();
#else
      int err = errno;
#endif
      if (err == SocketError::ConnRefused)
        LOG(INFO) << "Connection refused. Will try again...";
      else
        LOG(INFO) << "connect() failed";

      close_socket(socketFd);
      socketFd = INVALID_SOCKET;
      continue;
    }

    break;  // success
  }

  freeaddrinfo(addrinfo);
  if (socketFd != INVALID_SOCKET)
    LOG(INFO) << "connected to display server";
}

IPC_Channel::~IPC_Channel()
{
  if (connected())
    disconnect();

  if (--num_active_channels == 0) {
#ifdef IS_WINDOWS
    WSACleanup();
#endif
  }
}

void IPC_Channel::disconnect()
{
  CHECK(connected());

  close_socket(socketFd);
  socketFd = INVALID_SOCKET;
}

bool IPC_Channel::send(p_std::span<const uint8_t> message)
{
  if (!connected()) {
    connect();
    if (!connected())
      return false;
  }

  // Start with the length of the message.
  // FIXME: annoying coupling w/sending code's message buffer layout...
  int *startPtr = (int *)message.data();
  *startPtr = message.size();

  int bytesSent = ::send(socketFd, (const char *)message.data(), message.size(), 0 /* flags */);
  if (bytesSent == message.size())
    return true;

  LOG(ERROR) << "send to display server failed";
  disconnect();
  return false;
}

namespace {
enum Display_Directive : uint8_t {
  Open_Image = 0,
  Reload_Image = 1,
  Close_Image = 2,
  Update_Image = 3,
  Create_Image = 4,
};

void serialize(uint8_t **ptr, const std::string &s)
{
  for (size_t i = 0; i < s.size(); ++i, *ptr += 1)
    **ptr = s[i];
  **ptr = '\0';
  *ptr += 1;
}

template<typename T> void serialize(uint8_t **ptr, T value)
{
  memcpy(*ptr, &value, sizeof(T));
  *ptr += sizeof(T);
}

}  // namespace

class Display_Item {
 public:
  Display_Item(const std::string &title,
               int2 resolution,
               int2 tileSize,
               bool tiled,
               bool reuse_image,
               std::vector<std::string> channel_names,
               std::function<void(int4 b, Display_Pair &p, p_std::span<p_std::span<float>>)>
                   get_tile_values);

  bool Display(IPC_Channel &channel);
  std::queue<Display_Pair> params;

 private:
  bool send_open_image(IPC_Channel &channel);

  bool openedImage = false;
  bool tiled = false;
  std::string title;
  int2 resolution;
  std::function<void(int4 b, Display_Pair &p, p_std::span<p_std::span<float>>)> get_tile_values;
  std::vector<std::string> channel_names;
  int2 tileSize = make_int2(128, 128);

  struct Image_Channel_Buffer {
    Image_Channel_Buffer(const std::string &channelName,
                         int nTiles,
                         int2 tileSize,
                         const std::string &title);

    void set_tile_bounds(int x, int y, int width, int height);
    bool send_if_changed(IPC_Channel &channel, int tileIndex, int2 tileSize);

    std::vector<uint8_t> buffer;
    int tileBoundsOffset = 0, channelValuesOffset = 0;
    std::vector<uint64_t> tileHashes;

    int setCount;
  };
  std::vector<Image_Channel_Buffer> channel_buffers;
};

Display_Item::Display_Item(
    const std::string &baseTitle,
    int2 resolution,
    int2 tileSize,
    bool tiled,
    bool reuse_image,
    std::vector<std::string> channel_names,
    std::function<void(int4 b, Display_Pair &p, p_std::span<p_std::span<float>>)> get_tile_values)

    : openedImage(reuse_image),
      tiled(tiled),
      resolution(resolution),
      get_tile_values(get_tile_values),
      channel_names(channel_names),
      tileSize(tileSize)
{
  std::stringstream stream_buf;
  stream_buf << baseTitle << " (";
#ifdef IS_WINDOWS
  stream_buf << get_current_thread_id();
#else
  stream_buf << getpid();
#endif
  stream_buf << ")";
  title = stream_buf.str();

  int nTiles = ((resolution.x + tileSize.x - 1) / tileSize.x) *
               ((resolution.y + tileSize.y - 1) / tileSize.y);

  for (const std::string &channelName : channel_names)
    channel_buffers.push_back(Image_Channel_Buffer(channelName, nTiles, tileSize, title));
}

Display_Item::Image_Channel_Buffer::Image_Channel_Buffer(const std::string &channelName,
                                                         int nTiles,
                                                         int2 tileSize,
                                                         const std::string &title)
{
  int bufferAlloc = tileSize.x * tileSize.y * sizeof(float) + title.size() + 32;

  buffer.resize(bufferAlloc);

  uint8_t *ptr = buffer.data();
  serialize(&ptr, int(0));  // reserve space for message length
  serialize(&ptr, Display_Directive::Update_Image);
  uint8_t grabFocus = 0;
  serialize(&ptr, grabFocus);
  serialize(&ptr, title);
  serialize(&ptr, channelName);

  tileBoundsOffset = ptr - buffer.data();
  // Note: may not be float-aligned, but that's not a problem on x86...
  // TODO: fix this. The problem is that it breaks the whole idea of
  // passing a span<float> to the callback function...
  channelValuesOffset = tileBoundsOffset + 4 * sizeof(int);

  // Zero-initialize the buffer color contents before computing the hash
  // for a fully-zero tile (which corresponds to the initial state on the
  // viewer side.)
  memset(buffer.data() + channelValuesOffset, 0, tileSize.x * tileSize.y * sizeof(float));
  uint64_t zeroHash = util_murmur_hash3(
      buffer.data() + channelValuesOffset, tileSize.x * tileSize.y * sizeof(float), 0);
  tileHashes.assign(nTiles, zeroHash);
}

void Display_Item::Image_Channel_Buffer::set_tile_bounds(int x, int y, int width, int height)
{
  uint8_t *ptr = buffer.data() + tileBoundsOffset;

  serialize(&ptr, x);
  serialize(&ptr, y);
  serialize(&ptr, width);
  serialize(&ptr, height);

  setCount = width * height;
}

bool Display_Item::Image_Channel_Buffer::send_if_changed(IPC_Channel &ipc_channel,
                                                         int tileIndex,
                                                         int2 tileSize)
{
  int excess = setCount - tileSize.x * tileSize.y;
  if (excess > 0)
    memset(
        buffer.data() + channelValuesOffset + setCount * sizeof(float), 0, excess * sizeof(float));

  uint64_t hash = util_murmur_hash3(
      buffer.data() + channelValuesOffset, tileSize.x * tileSize.y * sizeof(float), 0);
  if (hash == tileHashes[tileIndex])
    return true;

  if (!ipc_channel.send(
          p_std::make_span(buffer.data(), channelValuesOffset + setCount * sizeof(float))))
    return false;

  tileHashes[tileIndex] = hash;
  return true;
}

bool Display_Item::Display(IPC_Channel &ipc_channel)
{
  if (!openedImage) {
    if (!send_open_image(ipc_channel))
      // maybe next time
      return false;
    openedImage = true;
  }

  while (!params.empty()) {
    std::vector<p_std::span<float>> display_values(channel_buffers.size());
    for (int c = 0; c < channel_buffers.size(); ++c) {
      float *ptr = (float *)(channel_buffers[c].buffer.data() +
                             channel_buffers[c].channelValuesOffset);
      display_values[c] = p_std::make_span(ptr, tileSize.x * tileSize.y);
    }

    Display_Pair &param = params.front();

    if (tiled) {
      int x = param.first.full_offset.x;
      int y = param.first.full_offset.y;
      int height = std::min(y + tileSize.y, resolution.y) - y;
      int width = std::min(x + tileSize.x, resolution.x) - x;

      for (int c = 0; c < channel_buffers.size(); ++c)
        channel_buffers[c].set_tile_bounds(
            x, (param.first.full_size.y - tileSize.y) - y, width, height);

      int4 b = make_int4(0, 0, width, height);
      get_tile_values(b, param, p_std::make_span(display_values));

      int tileIndex = (y / tileSize.y) * ((resolution.x + tileSize.x - 1) / tileSize.x) + (x / tileSize.x);

      // send the RGB buffers only if they're different than
      // the last version sent.
      for (int c = 0; c < channel_buffers.size(); ++c)
        if (!channel_buffers[c].send_if_changed(ipc_channel, tileIndex, tileSize)) {
          // Welp. Stop for now...
          openedImage = false;
          return false;
        }
    }
    else {
      int tileIndex = 0;
      for (int y = 0; y < resolution.y; y += tileSize.y)
        for (int x = 0; x < resolution.x; x += tileSize.x, ++tileIndex) {
          int height = std::min(y + tileSize.y, resolution.y) - y;
          int width = std::min(x + tileSize.x, resolution.x) - x;

          for (int c = 0; c < channel_buffers.size(); ++c)
            channel_buffers[c].set_tile_bounds(x, y, width, height);

          int4 b = make_int4(x, y, x + width, y + height);
          get_tile_values(b, param, p_std::make_span(display_values));

          // send the RGB buffers only if they're different than
          // the last version sent.
          for (int c = 0; c < channel_buffers.size(); ++c)
            if (!channel_buffers[c].send_if_changed(ipc_channel, tileIndex, tileSize)) {
              // Welp. Stop for now...
              openedImage = false;
              return false;
            }
        }
    }

    params.pop();
  }

  return true;
}

bool Display_Item::send_open_image(IPC_Channel &ipc_channel)
{
  // Initial "open the image" message
  uint8_t buffer[1024];
  uint8_t *ptr = buffer;

  serialize(&ptr, int(0));  // reserve space for message length
  serialize(&ptr, Display_Directive::Create_Image);
  uint8_t grabFocus = 1;
  serialize(&ptr, grabFocus);
  serialize(&ptr, title);

  int nChannels = channel_names.size();
  serialize(&ptr, resolution.x);
  serialize(&ptr, resolution.y);
  serialize(&ptr, nChannels);
  for (int c = 0; c < nChannels; ++c)
    serialize(&ptr, channel_names[c]);

  return ipc_channel.send(p_std::make_span(buffer, ptr - buffer));
}

static std::atomic<bool> exitThread{false};
static std::mutex mutex;
static std::thread updateThread;
static std::list<Display_Item *> dynamicItems;

static IPC_Channel *channel;

static void update_dynamic_items()
{
  int wait_time = 250;
  int retries = 0;

  while (!exitThread) {
    std::this_thread::sleep_for(std::chrono::milliseconds(wait_time));

    std::lock_guard<std::mutex> lock(mutex);
    for (auto *item : dynamicItems)
      item->Display(*channel);

    if (!channel->connected()) {
      retries++;
      if (retries == 3) {
        if (wait_time < 10000) {
          wait_time *= 4;
          retries = 0;
        }
      }
    }
    else if (wait_time > 250)
      wait_time = 250;
  }

  // One last time to get the last bits
  std::lock_guard<std::mutex> lock(mutex);
  for (auto *item : dynamicItems)
    item->Display(*channel);

  dynamicItems.clear();
  delete channel;
  channel = nullptr;
}

/* --------------------------------------------------------------------
 * TEVDisplayDriver.
 */

TEVDisplayDriver::TEVDisplayDriver(std::string display_server) : _display_server(display_server)
{
  connect_to_display_server(_display_server);
}

TEVDisplayDriver::~TEVDisplayDriver()
{
  disconnect_from_display_server();
}

/* --------------------------------------------------------------------
 * Update procedure.
 */

void TEVDisplayDriver::update()
{
}

void TEVDisplayDriver::next_tile_begin()
{
  /* Assuming no tiles used in interactive display. */
}

bool TEVDisplayDriver::update_begin(const Params &params, int texture_width, int texture_height)
{
  /* Note that it's the responsibility of TEVDisplayDriver to ensure updating and drawing
   * the texture does not happen at the same time. This is achieved indirectly.
   *
   * This locking is not performed on the Cycles side, because that would cause lock inversion. */

  std::lock_guard<std::mutex> lock(mutex);
  /* Update texture dimensions if needed. */
  if (texture_.full_width != params.size.x || texture_.full_height != params.size.y) {
    texture_.full_width = params.size.x;
    texture_.full_height = params.size.y;

    texture_.pixels.resize(texture_.full_width * texture_.full_height);

    int2 tile_size = make_int2(128, 128);
    if (params.size.x != params.full_size.x || params.size.y != params.full_size.y)
      tile_size = params.size;

    bool reuse_image = _current_item != nullptr;
    _current_item = new Display_Item(
        "Test",
        params.full_size,
        tile_size,
        (params.full_size.x != params.size.x || params.full_size.y != params.size.y),
        reuse_image,
        {"R", "G", "B"},
        [](int4 b, Display_Pair &p, p_std::span<p_std::span<float>> display_value) {
          DisplayDriver::Params &param = p.first;
          Display_Info &texture = p.second;
          int index = 0;
          int origin = (texture.height - 1) * texture.width;
          int x_stride = param.size.x / texture.width;
          int y_stride = param.size.y / texture.height;
          for (int y = b.y; y < b.w; ++y) {
            int yy = y / y_stride;
            for (int x = b.x; x < b.z; ++x) {
              int xx = x / x_stride;
              int offset = origin - yy * texture.width + xx;
              float4 val = half4_to_float4_image(texture.pixels[offset]);
              for (int j = 0; j < 3; ++j)
                display_value[j][index] = val[j];
              ++index;
            }
          }
        });
    dynamicItems.push_back(_current_item);
  }

  if (texture_.width != texture_width || texture_.height != texture_height) {
    texture_.width = texture_width;
    texture_.height = texture_height;

    /* Texture did change, and no pixel storage was provided. Tag for an explicit zeroing out to
     * avoid undefined content. */
    texture_.need_clear = true;
  }

  _item.first = params;

  /* New content will be provided to the texture in one way or another, so mark this in a
   * centralized place. */
  texture_.need_update = true;

  return true;
}

void TEVDisplayDriver::update_end()
{
  std::lock_guard<std::mutex> lock(mutex);
  _item.second = texture_;
  _current_item->params.push(_item);
}

/* --------------------------------------------------------------------
 * Texture buffer mapping.
 */

half4 *TEVDisplayDriver::map_texture_buffer()
{
  half4 *mapped_rgba_pixels = texture_.pixels.data();
  if (!mapped_rgba_pixels) {
    LOG(ERROR) << "Error mapping TEVDisplayDriver pixel buffer object.";
  }

  if (texture_.need_clear) {
    const int64_t texture_width = texture_.width;
    const int64_t texture_height = texture_.height;
    memset(reinterpret_cast<void *>(mapped_rgba_pixels),
           0,
           texture_width * texture_height * sizeof(half4));
    texture_.need_clear = false;
  }

  return mapped_rgba_pixels;
}

void TEVDisplayDriver::unmap_texture_buffer()
{
}

/* --------------------------------------------------------------------
 * Drawing.
 */

void TEVDisplayDriver::clear()
{
  texture_.need_clear = true;
}

void TEVDisplayDriver::draw(const Params &params)
{
  /* See do_update_begin() for why no locking is required here. */
  if (texture_.need_clear) {
    /* Texture is requested to be cleared and was not yet cleared.
     * Do early return which should be equivalent of drawing all-zero texture. */
    return;
  }
}

void TEVDisplayDriver::connect_to_display_server(const std::string &host)
{
  CHECK(channel == nullptr);
  channel = new IPC_Channel(host);

  updateThread = std::thread(update_dynamic_items);
}

void TEVDisplayDriver::disconnect_from_display_server()
{
  if (updateThread.get_id() != std::thread::id()) {
    exitThread = true;
    updateThread.join();
    updateThread = std::thread();
    exitThread = false;
  }
}

CCL_NAMESPACE_END
