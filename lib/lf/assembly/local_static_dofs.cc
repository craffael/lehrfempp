/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementation of description of location of local shape functions
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#include "local_static_dofs.h"

namespace lf::assemble {

size_type LocalStaticDOFs2D::NoLocDofs(lf::base::RefEl refel) const {
  switch (refel) {
    // internal indexing must match convention used in constructor
    case lf::base::RefEl::kPoint(): {
      return no_loc_dofs_per_refel_[0];
    }
    case lf::base::RefEl::kSegment(): {
      return no_loc_dofs_per_refel_[1];
    }
    case lf::base::RefEl::kTria(): {
      return no_loc_dofs_per_refel_[2];
    }
    case lf::base::RefEl::kQuad(): {
      return no_loc_dofs_per_refel_[3];
    }
    default: {
      LF_ASSERT_MSG(false, "Illegal refel type " << refel);
      break;
    }
  }
  return 0;
}  // end NoLocDofs

}  // namespace lf::assemble
