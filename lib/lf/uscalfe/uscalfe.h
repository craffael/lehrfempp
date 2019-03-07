#ifndef LF_FE_H
#define LF_FE_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Lagrangian finite elements: local definition and assembly
 * for elliptic BVPs
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "fe_space_lagrange_o1.h"
#include "fe_tools.h"
#include "lin_fe.h"
#include "loc_comp_norms.h"
#include "mesh_function_constant.h"
#include "mesh_function_global.h"
#include "mesh_function_traits.h"
#include "uniform_scalar_fe_space.h"

/**
 * @brief Collects data structures and algorithms designed for scalar finite
 * element methods primarily meant for second-order elliptic boundary value
 * problems.
 */
namespace lf::uscalfe {}  // namespace lf::uscalfe

#endif
