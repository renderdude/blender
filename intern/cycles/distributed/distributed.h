#ifndef DISTRIBUTED_H
#define DISTRIBUTED_H

#include <string>
#include <unistd.h>

#include "mpl/mpl.hpp"

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
  /// @name Measurement
  ///@{
  ///@}
  /// @name Comparison
  ///@{
  ///@}
  /// @name Status report
  ///@{
  bool is_distributed;
  ///@}
  /// @name Status setting
  ///@{
  ///@}
  /// @name Cursor movement
  ///@{
  ///@}
  /// @name Element change
  ///@{
  ///@}
  /// @name Removal
  ///@{
  ///@}
  /// @name Resizing
  ///@{
  ///@}
  /// @name Transformation
  ///@{
  ///@}
  /// @name Conversion
  ///@{
  std::string to_string() const;
  ///@}
  /// @name Basic operations
  ///@{
  ///@}
  /// @name Miscellaneous
  ///@{
  ///@}
  /// @name Obsolete
  ///@{
  ///@}
  /// @name Inapplicable
  ///@{
  ///@}

protected:
  ///@{
  ///@}

private:
  /* data */

}; // end of class Distributed
CCL_NAMESPACE_END

#endif // DISTRIBUTED_H
