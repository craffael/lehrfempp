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
 * @author Tobias Rohner
 * @date May 2020
 * @copyright MIT License
 */

#include "fe_point.h"
#include "fe_tools.h"
#include "hierarchic_scalar_fe_space.h"
#include "loc_comp_ellbvp.h"
#include "mesh_function_fe.h"
#include "mesh_function_grad_fe.h"
#include "prolongation.h"
#include "scalar_fe_space.h"
#include "scalar_reference_finite_element.h"

/**
 * @brief Collects data structures and algorithms designed for scalar finite
 * element methods primarily meant for second-order elliptic boundary value
 * problems.
 *
 * This namespace contains a number of classes/functions which
 * can be used to solve boundary vlaue problems with scalar finite elements.
 * This means that the shape functions must be scalar valued, but the shape
 * functions of a given approximation space may depend on the location in
 * the mesh instead of only the corresponding reference element. The
 * `lf::uscalfe` namespace is a specialization of this namespace to uniform
 * scalar finite elements. (Mainly Lagrangian FE).
 *
 * Examples of approximation spaces that the methodsclasses in this namespace
 * can represent/handle are:
 * - Lagrangian FE of any order
 * - Hierarchical FE Spaces
 * - Approximation spaces with local p-refinement
 * - Broken spaces (e.g. for Discontinuous Galerkin Approximations)
 */
namespace lf::fe {

/** Type for indices into global matrices/vectors */
using gdof_idx_t = lf::assemble::gdof_idx_t;
/** Type for indices referring to entity matrices/vectors */
using ldof_idx_t = lf::assemble::ldof_idx_t;
/** Type for vector length/matrix sizes */
using size_type = lf::assemble::size_type;
/** Type for (co-)dimensions */
using dim_t = lf::assemble::dim_t;
/** Type for global index of entities */
using glb_idx_t = lf::assemble::glb_idx_t;
/** Type for indexing sub-entities */
using sub_idx_t = lf::base::sub_idx_t;

// Import operators/free functions from lf::mesh::utils so we can apply them
// also to mesh functions defined in lf::fe (Argument Dependent Lookup)
using mesh::utils::operator*;
using mesh::utils::operator+;
using mesh::utils::operator-;
using mesh::utils::adjoint;
using mesh::utils::conjugate;
using mesh::utils::squaredNorm;
using mesh::utils::transpose;
}  // namespace lf::fe

#endif
