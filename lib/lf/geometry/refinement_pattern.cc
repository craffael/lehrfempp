/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementation of parallelogram test for lattice polygons
 * @author Ralf Hiptmair
 * @date Ocotber 2018
 * @copyright MIT License
 */

#include "refinement_pattern.h"

namespace lf::geometry {
bool isParallelogram(
    const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &polygon) {
  // A parallelogram must have four vertices
  if (polygon.cols() != 4) {
    return false;
  }
  return ((polygon * (Eigen::Vector4i() << 1, -1, -1, 1).finished())
              .squaredNorm() == 0);
}
}  // namespace lf::geometry
