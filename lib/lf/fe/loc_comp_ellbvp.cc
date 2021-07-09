/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementation of local computations for FE for
 * 2nd-order linear elliptic boundary value problems.
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "loc_comp_ellbvp.h"

namespace lf::fe {

std::shared_ptr<spdlog::logger>& DiffusionElementMatrixProviderLogger() {
  static auto logger =
      base::InitLogger("lf::fe::DiffusionElementMatrixProviderLogger");
  return logger;
}

std::shared_ptr<spdlog::logger>& MassElementMatrixProviderLogger() {
  static auto logger =
      base::InitLogger("lf::fe::MassElementMatrixProviderLogger");
  return logger;
}

std::shared_ptr<spdlog::logger>& MassEdgeMatrixProviderLogger() {
  static auto logger = base::InitLogger("lf::fe::MassEdgeMatrixProviderLogger");
  return logger;
}

std::shared_ptr<spdlog::logger>& ScalarLoadElementVectorProviderLogger() {
  static auto logger =
      base::InitLogger("lf::fe::ScalarLoadElementVectorProviderLogger");
  return logger;
}

std::shared_ptr<spdlog::logger>& ScalarLoadEdgeVectorProviderLogger() {
  static auto logger =
      base::InitLogger("lf::fe::ScalarLoadEdgeVectorProviderLogger");
  return logger;
}

}  // end namespace lf::fe
