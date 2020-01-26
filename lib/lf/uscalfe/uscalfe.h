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
#include "fe_space_lagrange_o2.h"
#include "fe_space_lagrange_o3.h"
#include "fe_tools.h"
#include "lin_fe.h"
#include "loc_comp_ellbvp.h"
#include "loc_comp_norms.h"
#include "mesh_function_fe.h"
#include "mesh_function_grad_fe.h"
#include "uniform_scalar_fe_space.h"

#include <lf/mesh/utils/utils.h>

/**
 * @brief Collects data structures and algorithms designed for scalar finite
 * element methods primarily meant for second-order elliptic boundary value
 * problems.
 *
 * This namespace contains a number of classes/functions which
 * can be used to solve boundary value problems with uniform,
 * scalar Finite elements:
 * - Uniform means that the approximation space has uniform order of
 *   approximation over the whole mesh. Or in other words: The shape functions
 *   of a given approximation space depend only on the underlying reference
 *   element of a mesh entity.
 * - Scalar means that the approximation space is always scalar valued, the
 *   shape functions are scalar valued.
 *
 * Examples of approximation spaces that the methods/classes in this namespace
 * can represent/handle are:
 * - n-th order Lagrangian Shape functions
 * - 1st order Crouzeix-Raviart approximation space
 * - Broken spaces (e.g. for Discontinuous Galerkin Approximations)
 *
 * Here are examples of use cases not supported by lf::uscalfe:
 * - Approximation spaces where the polynomial degree depends on the mesh,
 *   respectively hierarchic approximation spaces (e.g. for hp-fem)
 * - Non-Scalar valued approximation spaces such as EdgeFunctions/Nedelec
 *   elements or Thomas-Raviart elements.
 */
namespace lf::uscalfe {
// Import operators/free functions from lf::mesh::utils so we can apply them
// also to mesh functions defined in lf::uscalfe (Argument Dependent Lookup)
using mesh::utils::operator*;
using mesh::utils::operator+;
using mesh::utils::operator-;
using mesh::utils::squaredNorm;
using mesh::utils::transpose;
}  // namespace lf::uscalfe

#endif
