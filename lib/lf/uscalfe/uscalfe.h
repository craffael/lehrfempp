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
#include "mesh_function_binary.h"
#include "mesh_function_constant.h"
#include "mesh_function_fe.h"
#include "mesh_function_global.h"
#include "mesh_function_grad_fe.h"
#include "mesh_function_traits.h"
#include "mesh_function_unary.h"
#include "scalar_uniform_fe_space.h"

namespace lf::uscalfe {}  // namespace lf::uscalfe

#endif
