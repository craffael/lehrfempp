#ifndef LF_FETEST_H
#define LF_FETEST_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Functions for testing finite element facilities
 * @author Ralf Hiptmair
 * @date November 2018
 * @copyright MIT License
 */

#include <lf/refinement/mesh_hierarchy.h>
#include "fe_tools.h"

namespace lf::fe {
/**
 * @brief record interpolation errors in L2 norm and H1 norm on a sequence of
 * 2D hybrid meshes
 *
 * @tparam FFUNC functor type providing the scalar-valued function to be
 *               interpolated
 * @tparam GRADFUNC functor type for objects returning the gradient
 *
 * @param mesh_ptrs array of pointers to hybrid 2D meshes
 * @param f object encoding the function
 * @param grad_f gradient of f
 * @param rfs_tria_p pointer to description of local FE space on triangles
 * @param rfs_quad_p pointer to description of local FE space on quadrilaterals
 *
 * ### type requirement
 *
 * - FFUNC must provide an evaluation operator taking a single Eigen::VectorXd
 * argument and returning a scalar
 * - GRADFUNC must have an evaluation operator taking a single Eigen::VectorXd
 * argument and returning a vector.
 */
template <typename FFUNC, typename GRADFUNC>
std::vector<std::pair<double, double>> InterpolationErrors(
    std::vector<std::shared_ptr<const mesh::Mesh>> mesh_ptrs, FFUNC f,
    GRADFUNC grad_f,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p) {
  // Vector of error norms
  std::vector<std::pair<double, double>> err_norms{};

  // Helper class for L2 error computation
  LocalL2NormDifference lc_L2(rfs_tria_p, rfs_quad_p, f);
  // Helper class for H1 semi norm
  LocL2GradientFEDifference lc_H1(rfs_tria_p, rfs_quad_p, grad_f);

  // Loop over all meshes
  for (auto mesh_p : mesh_ptrs) {
    // Build finite element space and set up local-to-global index map
    UniformScalarFiniteElementSpace fe_space{mesh_p, rfs_tria_p, rfs_quad_p};
    const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
    // Perform (nodal) projection of the passed function onto the finite element
    // space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(fe_space, f, DefaultEntitySelector())};
    // Compute norms of interpolation error by means of numerical quadrature
    // whose order is controlled by the polynomials degree of the FE space
    double L2err = NormOfDifference(dofh, lc_L2, coeff_vec);
    double H1serr = NormOfDifference(dofh, lc_H1, coeff_vec);
    err_norms.emplace_back(L2err, H1serr);
  }
  return err_norms;
}

template <typename FFUNC, typename GRADFUNC>
inline std::vector<std::pair<double, double>> InterpolationErrors(
    lf::refinement::MeshHierarchy &multi_mesh, FFUNC f, GRADFUNC grad_f,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p) {
  return InterpolationErrors(multi_mesh.getMeshes(), f, grad_f, rfs_tria_p,
                             rfs_quad_p);
}

/**
 * @brief Track energy of FE-interpolated scalarr-valued function
 *         on a sequence of meshes
 *
 * @tparam SCALAR scalar type for the computations
 * @tparam FFUNC functor type providing the scalar-valued function to be
 *               interpolated ("reference function")
 * @tparam DIFF_COEFF functor type for diffusion coefficient
 * @tparam REAC_COEFF functor type for reaction coefficient
 *
 * @param mesh_ptrs array of pointers to hybrid 2D meshes
 * @param f object encoding the reference function
 * @param alpha object providing the diffusion coefficient
 * @param gamma object supplying the reaction coeffcient
 * @param rfs_tria_p pointer to description of local FE space on triangles
 * @param rfs_quad_p pointer to description of local FE space on quadrilaterals
 * @return vector of approximate energies on all meshes
 *
 * Given a scalar valued function, we compute its finite element projection
 * every mesh of the provided sequence. On every mesh we assemble the 
 * finite element Galerkin matrix for the diffusion coefficient `alpha`
 * and the reaction coefficient `gamma`, and use it to compute the 
 * energy norms of the interpolant. 
 *
 * This function can be used for debugging of finite element implementations
 * by comparing the energy norms of the interpolant with the exact energy
 * of the passed function.
 *
 * ### type requirements
 *
 * - FFUNC must provide an evaluation operator taking a single Eigen::VectorXd
 * argument and returning a scalar
 * - DIFF_COEFF must comply with `std::function<T(Eigen::Vector2d)>`, where `T`
 * is either a SCALAR compatible type of a matrix type like
 * `Eigen::Matrix<SCALAR,...>`.
 * - REAC_COEFF must behave like `std::function<SCALAR(Eigen::Vector2d)>`.
 *
 */
template <typename SCALAR, typename FFUNC, typename DIFF_COEFF,
          typename REAC_COEFF>
std::vector<double> EnergiesOfInterpolants(
    std::vector<std::shared_ptr<const mesh::Mesh>> mesh_ptrs, FFUNC f,
    DIFF_COEFF alpha, REAC_COEFF gamma,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p) {
  // Vector for returning the energies
  std::vector<double> energies{};

  // Loop over all meshes
  for (auto mesh_p : mesh_ptrs) {
    // Build finite element space and set up local-to-global index map
    UniformScalarFiniteElementSpace fe_space{mesh_p, rfs_tria_p, rfs_quad_p};
    const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

    // I: Perform (nodal) projection of the passed function onto the finite
    // element space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(fe_space, f, DefaultEntitySelector())};

    // II: Assemble finite element Galerkin matrix
    // Dimension of finite element space`
    const lf::assemble::size_type N_dofs(dofh.NoDofs());
    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
    // Actual assembly
    lf::fe::SecOrdBVPLagrFEFullInteriorGalMat(fe_space, alpha, gamma, A);

    // Computation of energy norm of interpolant
    double energy = coeff_vec.dot(A.MatVecMult(1.0, coeff_vec));
    energies.push_back(energy);
  }
  return energies;
}

template <typename SCALAR, typename FFUNC, typename DIFF_COEFF,
          typename REAC_COEFF>
std::vector<double> EnergiesOfInterpolants(
    lf::refinement::MeshHierarchy &multi_mesh, FFUNC f, DIFF_COEFF alpha,
    REAC_COEFF gamma,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p) {
  return EnergiesOfInterpolants<SCALAR>(multi_mesh.getMeshes(), f, alpha, gamma,
                                rfs_tria_p, rfs_quad_p);
}

}  // namespace lf::fe

#endif
