/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementation of local computations for Lagrange FE for
 * 2nd-order linear elliptic boundary value problems.
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "loc_comp_ellbvp.h"

namespace lf::uscalfe {

std::shared_ptr<spdlog::logger>&
ReactionDiffusionElementMatrixProviderLogger() {
  static auto logger = base::InitLogger(
      "lf::uscalfe::ReactionDiffusionElementMatrixProviderLogger");
  return logger;
}

std::shared_ptr<spdlog::logger>& MassEdgeMatrixProviderLogger() {
  static auto logger =
      base::InitLogger("lf::uscalfe::MassEdgeMatrixProviderLogger");
  return logger;
}

std::shared_ptr<spdlog::logger>& ScalarLoadElementVectorProviderLogger() {
  static auto logger =
      base::InitLogger("lf::uscalfe::ScalarLoadElementVectorProviderLogger");
  return logger;
}

std::shared_ptr<spdlog::logger>& ScalarLoadEdgeVectorProviderLogger() {
  static auto logger = base::InitLogger("ScalarLoadEdgeVectorProviderLogger");
  return logger;
}

}  // end namespace lf::uscalfe
