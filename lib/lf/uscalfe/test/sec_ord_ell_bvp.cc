/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementation of data structures and functions relevant for
 *        second-order elliptic boundary value problems.
 * 2nd-order linear elliptic boundary value problems.
 * @author Ralf Hiptmair
 * @date November 2018
 * @copyright MIT License
 */

#include "sec_ord_ell_bvp.h"

namespace lf::uscalfe::test {
// Definition of output control variable
unsigned int LFELinSys_ctrl = 0;
}  // namespace lf::uscalfe::test
