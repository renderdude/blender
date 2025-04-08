#pragma once

#include <string>
#include <unistd.h>

#include "mpl/mpl.hpp"
#include "mpl/tag.hpp"

CCL_NAMESPACE_BEGIN

/** @brief distributed.
 * @details
 */

class Distributed {
 public:
  /// @name Initialization
  ///@{
  Distributed(bool reverse_connect);
  ///@}

  /// @name Duplication
  ///@{
  Distributed(Distributed const &) = default;
  Distributed &operator=(Distributed const &) = default;
  ///@}

  /// @name Move
  ///@{
  Distributed(Distributed &&) = default;
  Distributed &operator=(Distributed &&) = default;
  ///@}

  /// @name Destruction
  ///@{
  ~Distributed() = default;
  ///@}

  /// @name Access
  ///@{
  mpl::communicator comm_world, inter_comm_world;
  ///@}

  /// @name Status report
  ///@{
  bool is_render_server;
  ///@}
  /// @name MPI Tags
  ///@{
  enum class message_tags {
    base_message = 100,
    asset_checksum,
    chunk_count,
    chunk_size,
    chunk_data,
    file_count,
    file_name,
    request_asset_file,
    request_asset_checksum,
    request_rib_file,
    requested_asset_file,
    rib_checksum,
  };
  ///@}

  /// @name Conversion
  ///@{
  std::string to_string() const;
  ///@}

 protected:
  ///@{
  mpl::inter_communicator *_inter_comm;
  ///@}

 private:
  /* data */

};  // end of class Distributed
CCL_NAMESPACE_END
