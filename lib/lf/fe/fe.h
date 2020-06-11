#ifndef LF_FE_H
#define LF_FE_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Finite elements: local definition and assembly for elliptic BVPs
 * @author Tobias
 * @date May 2020
 * @copyright MIT License
 */

#include "fe_space_hierarchic.h"
#include "fe_tools.h"
#include "loc_comp_ellbvp.h"
#include "mesh_function_fe.h"
#include "mesh_function_grad_fe.h"

#include <lf/mesh/utils/utils.h>

/**
 * @brief Collects data structures and algorithms designed for scalar finite
 * element methods primarily meant for second-order elliptic boundary value
 * problems.
 *
 * TODO: ADD DOCUMENTATION
 */
namespace lf::fe {}  // namespace lf::fe

#endif
