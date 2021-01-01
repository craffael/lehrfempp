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
#include <lf/uscalfe/uscalfe.h>

namespace lf::uscalfe::test {

/** @brief Incremental assembly of global finite element Galerkin matrix
 *
 * @tparam SCALAR a scalar type
 * @tparam TMPMATRIX matrix type suitable for assembly
 * @tparam DIFF_COEFF a functor providing point evaluation for the diffusion
 * tensor
 * @tparam REACTION_COEFF a functor for point evaluation of the reaction
 * coefficient
 * @param fe_space a Lagrangian finite element space of uniform polynomial
 * degree
 * @param alpha diffusion coefficient
 * @param gamma reaction coefficient
 * @param A a mutable reference to the Galerkin matrix. This argument is used to
 *        return the matrix. The new entries will be added to any previous set
 *        entries.
 *
 * ### Template parameter type requirements
 * - TMPMATRIX is a rudimentary matrix type and must
 *   + provide a constructor taking two matrix dimension arguments
 *   + have a method `AddtoEntry(i,j,value_to_add)` for adding to a matrix entry
 *   A model type is lf::assemble::COOMatrix.
 * - DIFF_COEFF must comply with `std::function<T(Eigen::Vector2d)>`, where `T`
 * is either a SCALAR compatible type of a matrix type like
 * `Eigen::Matrix<SCALAR,...>`.
 * - REACTION_COEFF must behave like `std::function<SCALAR(Eigen::Vector2d)>`.
 *
 * This function can be used to assemble the global finite element Galerkin
 * matrix for a second-order elliptic boundary value problem. Essential
 * conditions are **not** taken into account.
 *
 * @note the matrix A passed as the last argument is _updated_ by the function.
 */
template <typename SCALAR, typename TMPMATRIX, typename DIFF_COEFF,
          typename REACTION_COEFF>
void SecOrdBVPLagrFEFullInteriorGalMat(
    std::shared_ptr<UniformScalarFESpace<SCALAR>> fe_space, DIFF_COEFF alpha,
    REACTION_COEFF gamma, TMPMATRIX &A) {
  using scalar_t = typename TMPMATRIX::Scalar;
  // The underlying finite element mesh
  const lf::mesh::Mesh &mesh{*fe_space->Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  // Object taking care of local computations. No selection of a subset
  // of cells is specified.
  ReactionDiffusionElementMatrixProvider<scalar_t, decltype(alpha),
                                         decltype(gamma)>
      elmat_builder(fe_space, alpha, gamma);
  // Invoke assembly on cells
  AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
}

/** @brief Incremental assembly of parts of a finite element Galerkin matrix
 *         due to boundary contributions
 *
 * @tparam SCALAR a scalar type
 * @tparam TMPMATRIX matrix type suitable for assembly
 * @tparam COEFF a functor providing the impedance function
 * @tparam EDGESELECTOR predicate for selection of edges
 *
 * @param fe_space a Lagrangian finite element space of uniform polynomial
 *         degree
 * @param eta impedance coefficient
 * @param edge_selector object with an evaluation operator returning true
 *        for all edges on the impedance boundary part.
 * @param A a mutable reference to the Galerkin matrix. This argument is used to
 *        return the matrix. The new entries will be added to any previous set
 *        entries.
 *
 * ### Template parameter type requirements
 * - TMPMATRIX is a rudimentary matrix type and must
 *   + provide a constructor taking two matrix dimension arguments
 *   + have a method `AddtoEntry(i,j,value_to_add)` for adding to a matrix entry
 *   A model type is lf::assemble::COOMatrix.
 * - COEFF must behave like `std::function<SCALAR(Eigen::Vector2d)>`.
 * - EDGESELECTOR needs an evaluation operator
 *                `std::function<bool(const Entity &)>`
 *
 * This function can be used to assemble the contribution from an impedance
 * boundary part to the Galerkin matrix for a second order elliptic boundary
 * value problem.
 */
template <typename SCALAR, typename TMPMATRIX, typename COEFF,
          typename EDGESELECTOR>
void SecOrdBVPLagrFEBoundaryGalMat(
    std::shared_ptr<UniformScalarFESpace<SCALAR>> fe_space, COEFF eta,
    EDGESELECTOR edge_sel, TMPMATRIX &A) {
  using scalar_t = typename TMPMATRIX::Scalar;
  // The underlying finite element mesh
  const lf::mesh::Mesh &mesh{*fe_space->Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  // Object taking care of local computations.
  MassEdgeMatrixProvider<scalar_t, decltype(eta), decltype(edge_sel)>
      edgemat_builder(fe_space, eta, edge_sel);
  // Invoke assembly on edges by specifying co-dimension = 1
  AssembleMatrixLocally(1, dofh, dofh, edgemat_builder, A);
}

/**
 * @brief Incremental assembly of the right-hand side vector for a second-order
 *        elliptic boundary value problem based on Lagrangian Finite Elements
 *
 * @tparam SCALAR a scalar type
 * @tparam VECTOR a generic vector type with component access through []
 * @tparam FUNCTOR object with an evaluation operator of signature
 *         std::function<SCALAR(const Eigen::VectorXd &)>, which supplies
 *         the source function
 * @param fe_space underlying finite element space providing index mappings and
 * mesh
 * @param f functor object for source function
 * @param phi mutable reference to a vector with scalar entries
 *
 * This function relies on the the class \ref ScalarLoadElementVectorProvider
 * for local computations and the function \ref
 * lf::assemble::AssembleVectorLocally() for assembly.
 *
 * @note the functions performs an update of the vector
 */
template <typename SCALAR, typename VECTOR, typename FUNCTOR>
void LagrFEVolumeRightHandSideVector(
    std::shared_ptr<UniformScalarFESpace<SCALAR>> fe_space, FUNCTOR f,
    VECTOR &phi) {
  using scalar_t = typename VECTOR::value_type;
  // The underlying finite element mesh
  const lf::mesh::Mesh &mesh{*fe_space->Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  // Object taking care of local computations. No selection of a subset
  // of cells is specified.
  ScalarLoadElementVectorProvider<scalar_t, FUNCTOR> elvec_builder(fe_space, f);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);
}

/**
 * @brief Incremental assembly of the boundary contributions to the right-hand
 * side vector for a second-order elliptic boundary value problem based on
 * Lagrangian Finite Elements
 *
 * @tparam SCALAR a scalar type
 * @tparam VECTOR a generic vector type with component access through []
 * @tparam FUNCTOR object with an evaluation operator of signature
 *         std::function<SCALAR(const Eigen::VectorXd &)>, which supplies
 *         data on a set of active edges.
 * @tparam EDGESELECTOR predicate for selection of edges
 *
 * @param fe_space underlying finite element space providing index mappings and
 * mesh
 * @param data functor object for boundary data
 * @param edge_selector object with an evaluation operator returning true
 *        for all edges on the impedance boundary part.
 * @param phi mutable reference to a vector with scalar entries
 *
 * This function relies on the the class \ref ScalarLoadEdgeMatrixProvider
 * for local computations and the function \ref
 * lf::assemble::AssembleVectorLocally() for assembly.
 *
 * @note the functions performs an update of the vector
 */
template <typename SCALAR, typename VECTOR, typename FUNCTOR,
          typename EDGESELECTOR>
void LagrFEBoundaryRightHandSideVector(
    std::shared_ptr<UniformScalarFESpace<SCALAR>> fe_space, FUNCTOR data,
    EDGESELECTOR edge_sel, VECTOR &phi) {
  using scalar_t = typename VECTOR::value_type;
  // The underlying finite element mesh
  const lf::mesh::Mesh &mesh{*fe_space->Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  // Object taking care of local computations. A predicate selects the edges to
  // be processed
  ScalarLoadEdgeVectorProvider<scalar_t, FUNCTOR, EDGESELECTOR> elvec_builder(
      fe_space, data, edge_sel);
  // Invoke assembly on edges (codim == 1), update vector
  AssembleVectorLocally(1, dofh, elvec_builder, phi);
}

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

  // Loop over all meshes
  for (auto mesh_p : mesh_ptrs) {
    // Build finite element space and set up local-to-global index map
    auto fe_space_p = std::make_shared<UniformScalarFESpace<double>>(
        mesh_p, rfs_tria_p, rfs_quad_p);

    const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
    // Perform (nodal) projection of the passed function onto the finite element
    // space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(*fe_space_p, f)};
    // Compute norms of interpolation error by means of numerical quadrature
    // whose order is controlled by the polynomials degree of the FE space
    auto mf_fe = MeshFunctionFE<double, double>(fe_space_p, coeff_vec);
    auto mf_grad_fe = MeshFunctionGradFE<double, double>(fe_space_p, coeff_vec);
    double L2err =
        std::sqrt(IntegrateMeshFunction(*mesh_p, squaredNorm(f - mf_fe), 2));
    double H1serr = std::sqrt(
        IntegrateMeshFunction(*mesh_p, squaredNorm(grad_f - mf_grad_fe), 2));
    err_norms.emplace_back(L2err, H1serr);
  }
  return err_norms;
}

template <typename FFUNC, typename GRADFUNC>
inline std::vector<std::pair<double, double>> InterpolationErrors(
    refinement::MeshHierarchy &multi_mesh, FFUNC f, GRADFUNC grad_f,
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
  static_assert(mesh::utils::isMeshFunction<DIFF_COEFF>);
  static_assert(mesh::utils::isMeshFunction<REAC_COEFF>);
  static_assert(mesh::utils::isMeshFunction<FFUNC>);
  // Vector for returning the energies
  std::vector<SCALAR> energies{};

  // Loop over all meshes
  for (auto mesh_p : mesh_ptrs) {
    auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
        mesh_p, rfs_tria_p, rfs_quad_p);
    // Build finite element space and set up local-to-global index map
    const assemble::DofHandler &dofh{fe_space->LocGlobMap()};

    // I: Perform (nodal) projection of the passed function onto the finite
    // element space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(*fe_space, f, base::PredicateTrue{})};

    // II: Assemble finite element Galerkin matrix
    // Dimension of finite element space`
    const lf::assemble::size_type N_dofs(dofh.NumDofs());
    // Matrix in triplet format holding Galerkin matrix, zero initially.
    assemble::COOMatrix<double> A(N_dofs, N_dofs);
    // Actual assembly
    SecOrdBVPLagrFEFullInteriorGalMat(fe_space, alpha, gamma, A);

    // Computation of energy norm of interpolant
    double energy = coeff_vec.dot(A.MatVecMult(1.0, coeff_vec));
    energies.push_back(energy);
  }
  return energies;
}

template <typename SCALAR, typename FFUNC, typename DIFF_COEFF,
          typename REAC_COEFF>
std::vector<SCALAR> EnergiesOfInterpolants(
    refinement::MeshHierarchy &multi_mesh, FFUNC f, DIFF_COEFF alpha,
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
    auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
        mesh_p, rfs_tria_p, rfs_quad_p, rfs_edge_p);
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

    // I: Collect flags for edges on the boundary
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
    auto bd_edge_sel = [&bd_flags,
                        &edge_sel](const lf::mesh::Entity &edge) -> bool {
      return (bd_flags(edge) && edge_sel(edge));
    };

    // II: Perform (nodal) projection of the passed function onto the finite
    // element space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(*fe_space, f, base::PredicateTrue{})};

    // III: Assemble finite element Galerkin matrix
    // Dimension of finite element space`
    const lf::assemble::size_type N_dofs(dofh.NumDofs());
    // Matrix in triplet format holding Galerkin matrix, zero initially.
    lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
    // Actual assembly
    SecOrdBVPLagrFEBoundaryGalMat(fe_space, eta, bd_edge_sel, A);

    // Computation of energy norm of interpolant
    double energy = coeff_vec.dot(A.MatVecMult(1.0, coeff_vec));
    energies.push_back(energy);
  }
  return energies;
}

template <typename SCALAR, typename FFUNC, typename IMP_COEFF,
          typename EDGESELECTOR>
std::vector<SCALAR> BoundaryEnergiesOfInterpolants(
    refinement::MeshHierarchy &multi_mesh, FFUNC f, IMP_COEFF eta,
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
    auto fe_space = std::make_shared<UniformScalarFESpace<SCALAR>>(
        mesh_p, rfs_tria_p, rfs_quad_p);
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

    // I: Perform (nodal) projection of the passed function onto the finite
    // element space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(*fe_space, v, base::PredicateTrue{})};

    // II: Assemble finite element right-hand-side vector
    // Dimension of finite element space`
    const lf::assemble::size_type N_dofs(dofh.NumDofs());
    Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_dofs);
    // Actual assembly
    LagrFEVolumeRightHandSideVector(fe_space, f, phi);
    // Evaluation of right-hand-side functional for interpolant
    double ell_val = coeff_vec.dot(phi);
    ell_vals.push_back(ell_val);
  }
  return ell_vals;
}

template <typename SCALAR, typename FFUNC, typename SOURCE_FUNC>
std::vector<SCALAR> RHSFunctionalForInterpolants(
    refinement::MeshHierarchy &multi_mesh, FFUNC v, SOURCE_FUNC f,
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
    auto fe_space = std::make_shared<UniformScalarFESpace<SCALAR>>(
        mesh_p, rfs_tria_p, rfs_quad_p, rfs_edge_p);
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

    // I: Collect flags for edges on the boundary
    auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
    auto bd_edge_sel = [&bd_flags,
                        &edge_sel](const lf::mesh::Entity &edge) -> bool {
      return (bd_flags(edge) && edge_sel(edge));
    };

    // II: Perform (nodal) projection of the passed function onto the finite
    // element space and obtain basis expansion coefficient vector
    auto coeff_vec{NodalProjection(*fe_space, v, base::PredicateTrue{})};

    // II: Assemble finite element right-hand-side vector
    // Dimension of finite element space`
    const lf::assemble::size_type N_dofs(dofh.NumDofs());
    Eigen::VectorXd phi = Eigen::VectorXd::Zero(N_dofs);
    // Actual assembly
    LagrFEBoundaryRightHandSideVector(fe_space, f, bd_edge_sel, phi);
    // Evaluation of right-hand-side functional for interpolant
    double ell_val = coeff_vec.dot(phi);
    ell_vals.push_back(ell_val);
  }
  return ell_vals;
}

template <typename SCALAR, typename FFUNC, typename SOURCE_FUNC,
          typename EDGESELECTOR>
std::vector<SCALAR> RHSBoundaryFunctionalForInterpolants(
    refinement::MeshHierarchy &multi_mesh, FFUNC v, SOURCE_FUNC f,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_tria_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_quad_p,
    std::shared_ptr<const ScalarReferenceFiniteElement<double>> rfs_edge_p,
    EDGESELECTOR edge_sel) {
  return RHSBoundaryFunctionalForInterpolants<SCALAR>(multi_mesh.getMeshes(), v,
                                                      f, rfs_tria_p, rfs_quad_p,
                                                      rfs_edge_p, edge_sel);
}

/**
 * @brief Wraps another ScalarReferenceFiniteElement and multiplies the shape
 * functions with the imaginary unit to create complex valued shape functions.
 * @tparam SCALAR Scalar type of the wrapped FiniteElement.
 */
template <class SCALAR>
class ComplexScalarReferenceFiniteElement
    : public ScalarReferenceFiniteElement<std::complex<SCALAR>> {
 public:
  ComplexScalarReferenceFiniteElement(
      std::unique_ptr<ScalarReferenceFiniteElement<SCALAR>> fe)
      : inner_(std::move(fe)) {}

  [[nodiscard]] base::RefEl RefEl() const override { return inner_->RefEl(); }
  [[nodiscard]] unsigned Degree() const override { return inner_->Degree(); }
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t subidx) const override {
    return inner_->NumRefShapeFunctions(codim, subidx);
  }
  [[nodiscard]] Eigen::Matrix<std::complex<SCALAR>, Eigen::Dynamic,
                              Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    return std::complex<SCALAR>(0, 1) *
           inner_->EvalReferenceShapeFunctions(refcoords);
  }
  [[nodiscard]] Eigen::Matrix<std::complex<SCALAR>, Eigen::Dynamic,
                              Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd &refcoords) const override {
    return std::complex<SCALAR>(0, 1) *
           inner_->GradientsReferenceShapeFunctions(refcoords);
  }
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return inner_->EvaluationNodes();
  }
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return inner_->NumEvaluationNodes();
  }

  [[nodiscard]] Eigen::Matrix<std::complex<SCALAR>, 1, Eigen::Dynamic>
  NodalValuesToDofs(const Eigen::Matrix<std::complex<SCALAR>, 1, Eigen::Dynamic>
                        &nodvals) const override {
    return inner_->NodalValuesToDofs(
        (nodvals / std::complex<SCALAR>(0, 1)).real());
  }

  [[nodiscard]] size_type NumRefShapeFunctions(dim_t dim) const override {
    return inner_->NumRefShapeFunctions(dim);
  }

 private:
  std::unique_ptr<const ScalarReferenceFiniteElement<SCALAR>> inner_;
};

/**
 * @brief Returns a UniformScalarFESpace that is made up of "complexified" (via
 * ComplexScalarReferenceFiniteElement) FeLagrangeO1 finite elements.
 */
inline std::shared_ptr<UniformScalarFESpace<std::complex<double>>>
MakeComplexLagrangeO1FeSpace(std::shared_ptr<const mesh::Mesh> mesh_p) {
  return std::make_shared<UniformScalarFESpace<std::complex<double>>>(
      mesh_p,
      std::make_shared<ComplexScalarReferenceFiniteElement<double>>(
          std::make_unique<FeLagrangeO1Tria<double>>()),
      std::make_shared<ComplexScalarReferenceFiniteElement<double>>(
          std::make_unique<FeLagrangeO1Quad<double>>()),
      std::make_shared<ComplexScalarReferenceFiniteElement<double>>(
          std::make_unique<FeLagrangeO1Segment<double>>()),
      nullptr);
}

}  // namespace lf::uscalfe::test

#endif
