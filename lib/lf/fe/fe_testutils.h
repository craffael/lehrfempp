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

#include <lf/mesh/utils/utils.h>
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
  MeshFunctionL2NormDifference lc_L2(rfs_tria_p, rfs_quad_p, f);
  // Helper class for H1 semi norm
  MeshFunctionL2GradientDifference lc_H1(rfs_tria_p, rfs_quad_p, grad_f);

  // Loop over all meshes
  for (auto mesh_p : mesh_ptrs) {
    // Build finite element space and set up local-to-global index map
    UniformScalarFiniteElementSpace fe_space{mesh_p, rfs_tria_p, rfs_quad_p};
    const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
    // Perform (nodal) projection of the passed function onto the finite element
    // space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(fe_space, f, base::PredicateTrue{})};
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
std::vector<SCALAR> EnergiesOfInterpolants(
    std::vector<std::shared_ptr<const mesh::Mesh>> mesh_ptrs, FFUNC f,
    DIFF_COEFF alpha, REAC_COEFF gamma,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p) {
  static_assert(isMeshFunction<DIFF_COEFF>);
  static_assert(isMeshFunction<REAC_COEFF>);
  static_assert(isMeshFunction<FFUNC>);
  // Vector for returning the energies
  std::vector<SCALAR> energies{};

  // Loop over all meshes
  for (auto mesh_p : mesh_ptrs) {
    // Build finite element space and set up local-to-global index map
    UniformScalarFiniteElementSpace fe_space{mesh_p, rfs_tria_p, rfs_quad_p};
    const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

    // I: Perform (nodal) projection of the passed function onto the finite
    // element space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(fe_space, f, base::PredicateTrue{})};

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
std::vector<SCALAR> EnergiesOfInterpolants(
    lf::refinement::MeshHierarchy &multi_mesh, FFUNC f, DIFF_COEFF alpha,
    REAC_COEFF gamma,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p) {
  return EnergiesOfInterpolants<SCALAR>(multi_mesh.getMeshes(), f, alpha, gamma,
                                        rfs_tria_p, rfs_quad_p);
}

/**
 * @brief Evaluate the contribution of boundary terms to the energy
 *        for finite element interpolants of a given function.
 *
 * @tparam SCALAR scalar type for the computations
 * @tparam FFUNC functor type providing the scalar-valued function to be
 *               interpolated ("reference function")
 * @tparam IMP_COEFF functor type for coefficient in impedance boundary
 condition
 * @tparam EDGESELECTOR predicate for the selection of active edges
 *
 * @param mesh_ptrs array of pointers to hybrid 2D meshes
 * @param f object encoding the reference function
 * @param eta object providing the impedance coefficient
 * @param rfs_tria_p pointer to description of local FE space on triangles
 * @param rfs_quad_p pointer to description of local FE space on quadrilaterals
 * @param rfs_edge_p pointer to description of local FE space on segments
 * @param edge_sel selector predicate object for relevant edges on the boundary
 * @return vector of approximate energies on all meshes
 *
 * Given a scalar valued function, we compute its finite element projection on
 * every mesh of the provided sequence. On every mesh we assemble the
 * finite element Galerkin matrix corresponding to the bilinear form
 * @f[
    (u,v) \mapsto \int\limits_{\Gamma}
 \eta(\mathbf{x})u\,v\,\mathrm{dS}(\mathbf{x})
 * @f]
 * where @f$\Gamma@f$ is the union of all edges of a mesh selected by the
 * predicate `edge_sel`.
 *
 * This function can be used for debugging of finite element implementations
 * by comparing the energy norms of the interpolant with the exact "boundary
 energy
 * contribution" of the passed function.
 *
 * ### type requirements
 *
 * - FFUNC must provide an evaluation operator taking a single Eigen::VectorXd
 * argument and returning a scalar
 * - IMP_COEFF must behave like `std::function<SCALAR(Eigen::Vector2d)>`.
 * - EDGESELECTOR needs an evaluation operator
 *                `std::function<bool(const Entity &)>`
 */
template <typename SCALAR, typename FFUNC, typename IMP_COEFF,
          typename EDGESELECTOR>
std::vector<SCALAR> BoundaryEnergiesOfInterpolants(
    std::vector<std::shared_ptr<const mesh::Mesh>> mesh_ptrs, FFUNC f,
    IMP_COEFF eta,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_edge_p,
    EDGESELECTOR edge_sel) {
  // Vector for returning the energies
  std::vector<SCALAR> energies{};

  // Loop over all meshes
  for (auto mesh_p : mesh_ptrs) {
    // Build finite element space and set up local-to-global index map
    UniformScalarFiniteElementSpace fe_space{mesh_p, rfs_tria_p, rfs_quad_p,
                                             rfs_edge_p};
    const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

    // I: Collect flags for edges on the boundary
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space.Mesh(), 1)};
    auto bd_edge_sel = [&bd_flags,
                        &edge_sel](const lf::mesh::Entity &edge) -> bool {
      return (bd_flags(edge) && edge_sel(edge));
    };

    // II: Perform (nodal) projection of the passed function onto the finite
    // element space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(fe_space, f, base::PredicateTrue{})};

    // III: Assemble finite element Galerkin matrix
    // Dimension of finite element space`
    const lf::assemble::size_type N_dofs(dofh.NoDofs());
    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
    // Actual assembly
    lf::fe::SecOrdBVPLagrFEBoundaryGalMat(fe_space, eta, bd_edge_sel, A);

    // Computation of energy norm of interpolant
    double energy = coeff_vec.dot(A.MatVecMult(1.0, coeff_vec));
    energies.push_back(energy);
  }
  return energies;
}

template <typename SCALAR, typename FFUNC, typename IMP_COEFF,
          typename EDGESELECTOR>
std::vector<SCALAR> BoundaryEnergiesOfInterpolants(
    lf::refinement::MeshHierarchy &multi_mesh, FFUNC f, IMP_COEFF eta,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_edge_p,
    EDGESELECTOR edge_sel) {
  return BoundaryEnergiesOfInterpolants<SCALAR>(multi_mesh.getMeshes(), f, eta,
                                                rfs_tria_p, rfs_quad_p,
                                                rfs_edge_p, edge_sel);
}

/**
 * @brief Evaluates a right hand side functional for the finite element
 *        interpolant of a given function
 *
 * @tparam SCALAR scalar type for the computations
 * @tparam FFUNC functor type providing the scalar-valued function to be
 *               interpolated ("reference function")
 * @tparam SOURCE_FUNC functor type for describing source function
 *
 * @param mesh_ptrs array of pointers to hybrid 2D meshes
 * @param v object encoding the reference function
 * @param f  object providing the source function @g$f\in L^2(\Omega)@f$
 * @param rfs_tria_p pointer to description of local FE space on triangles
 * @param rfs_quad_p pointer to description of local FE space on quadrilaterals
 * @return vector of approximate functional values on all meshes
 *
 * Given a scalar valued function, we compute its finite element projection
 * every mesh of the provided sequence. On every mesh we assemble the
 * finite element load vector and multiply it with the coefficient
 * vector of the finite element interpolant. This yields the value
 * of the right hand side functional.
 *
 * This function can be used for debugging of finite element implementations
 * by comparing the functional values for the interpolant with the exact value
 * for the passed function.
 *
 * ### type requirements
 *
 * - FFUNC must provide an evaluation operator taking a single Eigen::VectorXd
 * argument and returning a scalar
 * - SOURCE_FUNC must behave like `std::function<SCALAR(Eigen::Vector2d)>`.
 *
 */
template <typename SCALAR, typename FFUNC, typename SOURCE_FUNC>
std::vector<SCALAR> RHSFunctionalForInterpolants(
    std::vector<std::shared_ptr<const mesh::Mesh>> mesh_ptrs, FFUNC v,
    SOURCE_FUNC f,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p) {
  // Vector for returning the energies
  std::vector<SCALAR> ell_vals{};

  // Loop over all meshes
  for (auto mesh_p : mesh_ptrs) {
    // Build finite element space and set up local-to-global index map
    UniformScalarFiniteElementSpace fe_space{mesh_p, rfs_tria_p, rfs_quad_p};
    const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

    // I: Perform (nodal) projection of the passed function onto the finite
    // element space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(fe_space, v, base::PredicateTrue{})};

    // II: Assemble finite element right-hand-side vector
    // Dimension of finite element space`
    const lf::assemble::size_type N_dofs(dofh.NoDofs());
    Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_dofs);
    // Actual assembly
    lf::fe::LagrFEVolumeRightHandSideVector(fe_space, f, phi);
    // Evaluation of right-hand-side functional for interpolant
    double ell_val = coeff_vec.dot(phi);
    ell_vals.push_back(ell_val);
  }
  return ell_vals;
}

template <typename SCALAR, typename FFUNC, typename SOURCE_FUNC>
std::vector<SCALAR> RHSFunctionalForInterpolants(
    lf::refinement::MeshHierarchy &multi_mesh, FFUNC v, SOURCE_FUNC f,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p) {
  return RHSFunctionalForInterpolants<SCALAR>(multi_mesh.getMeshes(), v, f,
                                              rfs_tria_p, rfs_quad_p);
}

/**
 * @brief Evaluates boundary contributions to a right-hand-side functional for
 * the finite element interpolant of a given function
 *
 * @tparam SCALAR scalar type for the computations
 * @tparam FFUNC functor type providing the scalar-valued function to be
 *               interpolated ("reference function")
 * @tparam SOURCE_FUNC functor type for describing boundary source function
 * @tparam EDGESELECTOR predicate for the selection of active edges
 *
 * @param mesh_ptrs array of pointers to hybrid 2D meshes
 * @param v object encoding the reference function
 * @param f  object providing the boundary source function @g$f\in
 * L^2(\partial\Omega)@f$
 * @param rfs_tria_p pointer to description of local FE space on triangles
 * @param rfs_quad_p pointer to description of local FE space on quadrilaterals
 * @param rfs_edge_p pointer to description of local FE space on segments
 * @param edge_sel selector for relevant edges on the boundary
 * @return vector of approximate functional values on all meshes
 *
 * Given a scalar valued function, we compute its finite element projection
 * every mesh of the provided sequence. On every mesh we assemble the
 * boundary part of the finite element load vector and multiply it with the
 * coefficient vector of the finite element interpolant. This yields the value
 * of the boundary part of the right-hand-side functional.
 *
 * This function can be used for debugging of finite element implementations
 * in the case of inhomogeneous Neumann problems by comparing the functional
 * values for the interpolant with the exact value for the passed function.
 *
 * ### type requirements
 *
 * - FFUNC must provide an evaluation operator taking a single Eigen::VectorXd
 * argument and returning a scalar
 * - SOURCE_FUNC must behave like `std::function<SCALAR(Eigen::Vector2d)>`.
 *
 */
template <typename SCALAR, typename FFUNC, typename SOURCE_FUNC,
          typename EDGESELECTOR>
std::vector<SCALAR> RHSBoundaryFunctionalForInterpolants(
    std::vector<std::shared_ptr<const mesh::Mesh>> mesh_ptrs, FFUNC v,
    SOURCE_FUNC f,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_edge_p,
    EDGESELECTOR edge_sel) {
  // Vector for returning the energies
  std::vector<SCALAR> ell_vals{};

  // Loop over all meshes
  for (auto mesh_p : mesh_ptrs) {
    // Build finite element space and set up local-to-global index map
    UniformScalarFiniteElementSpace fe_space{mesh_p, rfs_tria_p, rfs_quad_p,
                                             rfs_edge_p};
    const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

    // I: Collect flags for edges on the boundary
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space.Mesh(), 1)};
    auto bd_edge_sel = [&bd_flags,
                        &edge_sel](const lf::mesh::Entity &edge) -> bool {
      return (bd_flags(edge) && edge_sel(edge));
    };

    // II: Perform (nodal) projection of the passed function onto the finite
    // element space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(fe_space, v, base::PredicateTrue{})};

    // II: Assemble finite element right-hand-side vector
    // Dimension of finite element space`
    const lf::assemble::size_type N_dofs(dofh.NoDofs());
    Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_dofs);
    // Actual assembly
    lf::fe::LagrFEBoundaryRightHandSideVector(fe_space, f, bd_edge_sel, phi);
    // Evaluation of right-hand-side functional for interpolant
    double ell_val = coeff_vec.dot(phi);
    ell_vals.push_back(ell_val);
  }
  return ell_vals;
}

template <typename SCALAR, typename FFUNC, typename SOURCE_FUNC,
          typename EDGESELECTOR>
std::vector<SCALAR> RHSBoundaryFunctionalForInterpolants(
    lf::refinement::MeshHierarchy &multi_mesh, FFUNC v, SOURCE_FUNC f,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_edge_p,
    EDGESELECTOR edge_sel) {
  return RHSBoundaryFunctionalForInterpolants<SCALAR>(multi_mesh.getMeshes(), v,
                                                      f, rfs_tria_p, rfs_quad_p,
                                                      rfs_edge_p, edge_sel);
}

}  // namespace lf::fe

#endif
