#pragma once

#include "ri_api.h"

/** @brief Ri_distributed.
 * @details
 * Distributed version of the Ri API. It's main function is to override the calls
 * of the API that deal with external resources, e.g., Archives, Shaders, Textures, ...
 * and ensure the correct version is local in the working directory on the remote server
 */
CCL_NAMESPACE_BEGIN

class Ri_Distributed : public Ri {
 public:
  /// @name Initialization
  ///@{
  Ri_Distributed(Options &options) : Ri(options) {}
  ///@}

  /// @name Destruction
  ///@{
  ~Ri_Distributed() = default;
  ///@}

  /// @name API Overrides
  ///@{
  ///@}

 protected:
  ///@{
  ///@}

 private:
  /* data */

};  // end of class Ri_Distributed

CCL_NAMESPACE_END
