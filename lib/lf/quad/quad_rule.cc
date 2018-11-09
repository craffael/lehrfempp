/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief No implementation here
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "quad_rule.h"

namespace lf::quad {

CONTROLDECLARECOMMENT(QuadRule, out_ctrl_, "out_ctrl_",
                      "Output control for QuadRule");

std::ostream& operator<<(std::ostream& stream,
                         const lf::quad::QuadRule& quadrule) {
  quadrule.PrintInfo(stream);
  return stream;
}

}  // namespace lf::quad
