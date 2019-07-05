/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementation of local assembly
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#include "assembler.h"

namespace lf::assemble {
// Debugging output control variable, switched off by default
ADDOPTION(ass_mat_dbg_ctrl, Assembly_ctrl,
          "Debugging output control for AssembleMatrixLocally()");

}  // namespace lf::assemble
