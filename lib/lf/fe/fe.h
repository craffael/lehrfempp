#ifndef LF_FE_H
#define LF_FE_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Finite elements: local definition and assembly
 * elliptic BVPs
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include "fe_space_lagrange_o1.h"
#include "fe_space_lagrange_uniform.h"
#include "fe_tools.h"
#include "lin_fe.h"
#include "loc_comp_norms.h"
#include "mesh_function_constant.h"
#include "mesh_function_global.h"
#include "mesh_function_traits.h"

namespace lf::fe {}  // namespace lf::fe

#endif
