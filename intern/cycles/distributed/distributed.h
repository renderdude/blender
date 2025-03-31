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
  static constexpr int request_checksum = 100;
  static constexpr int checksum = 101;
  static constexpr int request_file = 102;
  static constexpr int requested_file = 103;
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

}; // end of class Distributed
CCL_NAMESPACE_END
