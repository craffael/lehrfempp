#ifndef LF_FETOOLS_H
#define LF_FETOOLS_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @file
 * @brief Functions acting on finite element spaces
 * @author Ralf Hiptmair
 * @date November 2018
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>

#include "scalar_fe_space.h"

namespace lf::fe {

static const unsigned int kout_l2_qr = 1;
static const unsigned int kout_l2_rsfvals = 2;

namespace internal {

// TODO(raffael) convert this method into a lambda function of
// IntegrateMeshFunction() once
// https://developercommunity.visualstudio.com/content/problem/200017/if-constexpr-in-lambda.html
// is resolved.
template <class MF, class QR_SELECTOR>
auto LocalIntegral(const mesh::Entity &e, const QR_SELECTOR &qr_selector,
                   const MF &mf) -> mesh::utils::MeshFunctionReturnType<MF> {
  using MfType = mesh::utils::MeshFunctionReturnType<MF>;
  auto qr = qr_selector(e);
  auto values = mf(e, qr.Points());
  auto weights_ie =
      (qr.Weights().cwiseProduct(e.Geometry()->IntegrationElement(qr.Points())))
          .eval();
  LF_ASSERT_MSG(values.size() == qr.NumPoints(),
                "mf returns vector with wrong size.");
  if constexpr (std::is_arithmetic_v<MfType>) {  // NOLINT
    auto value_m = Eigen::Map<Eigen::Matrix<MfType, 1, Eigen::Dynamic>>(
        &values[0], 1, values.size());
    return (value_m * weights_ie)(0);
  }

  if constexpr (base::EigenMatrix<MfType>) {  // NOLINT
    constexpr int size = MfType::SizeAtCompileTime;
    if constexpr (size != Eigen::Dynamic) {
      auto value_m = Eigen::Map<
          Eigen::Matrix<typename MfType::Scalar, size, Eigen::Dynamic>>(
          &values[0](0, 0), size, values.size());
      MfType result;
      auto result_m =
          Eigen::Map<Eigen::Matrix<typename MfType::Scalar, size, 1>>(
              &result(0, 0));
      result_m = value_m * weights_ie;
      return result;
    }
  }
  // fallback: we cannot make any optimizations:
  MfType temp = weights_ie(0) * values[0];
  for (Eigen::Index i = 1; i < qr.NumPoints(); ++i) {
    temp = temp + weights_ie(i) * values[i];
  }
  return temp;
}
};  // namespace internal

/**
 * @brief Integrate a \ref mesh_function over a mesh (with quadrature rules)
 * @tparam MF The type of the \ref mesh_function "mesh function".
 * @tparam QR_SELECTOR The type of qr_selector (see below)
 * @tparam ENTITY_PREDICATE The type of the entity predicate (see below)
 * @param mesh The mesh to integrate over
 * @param mf The mesh function to integrate
 * @param qr_selector Provides the quadrature rule for every entity of the mesh.
 * @param ep Selects the entities over which `mf` is integrated (default: all
 * entities)
 * @param codim The codimension of the entities over which `mf` is integrated.
 * @return The integrated value
 *
 * ### Requirements for QR_SELECTOR
 * `QR_SELECTOR` should overload `operator()` as follows:
 * ```
 * quad::QuadRule operator()(const mesh::Entity& e) const
 * ```
 * i.e. it should return the quadrature rule for every entity `e` of the mesh
 * that is to be used for computing the integral of `mf` over `e`.
 *
 * ### Requirements for ENTITY_PREDICATE
 * The entity predicate should overload `operator()` as follows:
 * ```
 * bool operator()(const mesh::Entity& e) const
 * ```
 * It should return `true`, if `e` is part of the integration domain and `false`
 * if it is not.
 *
 * ### Example
 * @snippet fe_tools.cc integrateMeshFunction2
 */
template <class MF, class QR_SELECTOR,
          class ENTITY_PREDICATE = base::PredicateTrue,
          class = std::enable_if_t<
              std::is_invocable_v<QR_SELECTOR, const mesh::Entity &>>>
auto IntegrateMeshFunction(const lf::mesh::Mesh &mesh, const MF &mf,
                           const QR_SELECTOR &qr_selector,
                           const ENTITY_PREDICATE &ep = base::PredicateTrue{},
                           int codim = 0)
    -> mesh::utils::MeshFunctionReturnType<MF> {
  static_assert(mesh::utils::isMeshFunction<MF>);
  using MfType = mesh::utils::MeshFunctionReturnType<MF>;

  auto entities = mesh.Entities(codim);
  auto result = internal::LocalIntegral(**entities.begin(), qr_selector, mf);
  for (auto i = entities.begin() + 1; i != entities.end(); ++i) {
    if (!ep(**i)) {
      continue;
    }
    result = result + internal::LocalIntegral(**i, qr_selector, mf);
  }
  return result;
}

/**
 * @brief Integrate a \ref mesh_function "mesh function" over a mesh using
 * quadrature rules of uniform order.
 * @tparam MF type of \ref mesh_function "mesh function" to integrate
 * @tparam ENTITY_PREDICATE type of entity predicate (see below)
 * @param mesh The mesh over which `mf` is integrated.
 * @param mf The \ref mesh_function "mesh function" which is integrated
 * @param quad_degree The quadrature degree of the quadrature rules that are to
 * be used for integration. Internally Gauss-rules created by
 * `quad::make_QuadRule` are used.
 * @param ep The entity predicate selecting the entities over which `mf` is
 * integrated.
 * @param codim The codimension of the entities over which `mf` is integrated.
 * @return `mf` integrated over the entities `e` of `mf` where `ep(e)==true`.
 *
 * ### Requirements for ENTITY_PREDICATE
 * The entity predicate should overload `operator()` as follows:
 * ```
 * bool operator()(const mesh::Entity& e) const
 * ```
 * It should return `true`, if `e` is part of the integration domain and `false`
 * if it is not.
 *
 * ### Example
 * @snippet fe_tools.cc integrateMeshFunction
 */
template <class MF, class ENTITY_PREDICATE = base::PredicateTrue>
auto IntegrateMeshFunction(const lf::mesh::Mesh &mesh, const MF &mf,
                           int quad_degree,
                           const ENTITY_PREDICATE &ep = base::PredicateTrue{},
                           int codim = 0) {
  std::array<quad::QuadRule, 5> qrs;
  for (auto ref_el :
       {base::RefEl::kSegment(), base::RefEl::kTria(), base::RefEl::kQuad()}) {
    qrs[ref_el.Id()] = quad::make_QuadRule(ref_el, quad_degree);
  }
  return IntegrateMeshFunction(
      mesh, mf, [&](const mesh::Entity &e) { return qrs[e.RefEl().Id()]; }, ep,
      codim);
}

// ******************************************************************************

/**
 * @brief Computes nodal projection of a mesh function and returns the
 * finite element basis expansion coefficients of the result
 *
 * @tparam SCALAR a scalar type
 * @tparam MF a \ref mesh_function "MeshFunction" representing the scalar
 * valued function that should be projected
 * @tparam SELECTOR predicate type for selecting cells to be visited
 *
 * @param fe_space a uniform Lagrangian finite element space, providing
 * finite element specifications for the cells of the mesh
 * @param u functor object supplying a scalar-valued function that is to be
 * projected
 * @param pred predicate object for the selection of relevant cells
 * @return column vector of basis expansion coefficients for the resulting
 *         finite element function
 *
 * The implementation relies on the method
 * ScalarReferenceFiniteElement::NodalValuesToDofs(). Refer to its
 * documentation. This method is called for each active cell to **set** the
 * coefficients for the global shape functions associated with that cell.
 *
 * ### Example
 * @snippet fe_tools.cc nodalProjection
 */
template <typename SCALAR, typename MF, typename SELECTOR = base::PredicateTrue>
auto NodalProjection(const lf::fe::ScalarFESpace<SCALAR> &fe_space, MF &&u,
                     SELECTOR &&pred = base::PredicateTrue{}) {
  static_assert(mesh::utils::isMeshFunction<std::remove_reference_t<MF>>);
  // choose scalar type so it can hold the scalar type of u as well as
  // SCALAR
  using scalarMF_t =
      mesh::utils::MeshFunctionReturnType<std::remove_reference_t<MF>>;
  using scalar_t = decltype(SCALAR(0) * scalarMF_t(0));  // NOLINT
  // Return type, type for FE coefficient vector
  using dof_vec_t = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

  // Underlying mesh instance
  const lf::mesh::Mesh &mesh{*fe_space.Mesh()};
  // Fetch local-to-global index mapping for shape functions
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

  // Vector for returning expansion coefficients
  dof_vec_t glob_dofvec(dofh.NumDofs());
  glob_dofvec.setZero();

  // Loop over all cells
  for (const lf::mesh::Entity *cell : mesh.Entities(0)) {
    if (!pred(*cell)) {
      continue;
    }
    // Topological type of the cell
    const lf::base::RefEl ref_el{cell->RefEl()};

    // Information about local shape functions on reference element
    auto ref_shape_fns = fe_space.ShapeFunctionLayout(*cell);
    LF_ASSERT_MSG(ref_shape_fns, "reference shape function for "
                                     << ref_el << " not available.");
    // Obtain reference coordinates for evaluation nodes
    const Eigen::MatrixXd ref_nodes(ref_shape_fns->EvaluationNodes());

    // Collect values of function to be projected in a row vector
    auto uvalvec = u(*cell, ref_nodes);

    // Compute the resulting local degrees of freedom
    auto dofvec(ref_shape_fns->NodalValuesToDofs(
        Eigen::Map<Eigen::Matrix<scalar_t, 1, Eigen::Dynamic>>(
            &uvalvec[0], 1, uvalvec.size())));

    // Set the corresponing global degrees of freedom
    // Note: "Setting", not "adding to"
    // Number of local shape functions
    const size_type num_loc_dofs = dofh.NumLocalDofs(*cell);
    LF_ASSERT_MSG(
        dofvec.size() == num_loc_dofs,
        "Size mismatch: " << dofvec.size() << " <-> " << num_loc_dofs);
    // Fetch global numbers of local shape functions
    const std::span<const lf::assemble::gdof_idx_t> ldof_gidx(
        dofh.GlobalDofIndices(*cell));
    // Insert dof values into the global coefficient vector
    for (int j = 0; j < num_loc_dofs; ++j) {
      glob_dofvec[ldof_gidx[j]] = dofvec[j];
    }
  }
  return glob_dofvec;
}

/**
 * @brief Initialization of flags/values for dofs of a *Scalar*
 * finite element space whose values are imposed by a specified function.
 *
 * @tparam SCALAR scalar type of the basis functions of the finite element space
 * @tparam EDGESELECTOR predicate returning true for edges with fixed dofs
 * @tparam FUNCTION \ref mesh_function "MeshFunction" which defines the
 * imposed values on the edges
 *
 * @param fes The `ScalarFESpace` whose dofs should be fixed.
 * @param esscondflag predicate object whose evaluation operator returns
 *        true for all edges whose associated degrees of freedom should be
 *        set a fixed value.
 * @param g A scalar valued \ref mesh_function "MeshFunction" which describes
 * the values on the edges to which the dofs should be fixed.
 * @return a vector of flag-value pairs, a `true` first component indicating
 *         a fixed dof, with the second component providing the value in
 * this case.
 *
 * This function interpolates a scalar-valued function into a Scalar
 * finite element space restricted to a set of edges. It relies on the
 * method ScalarReferenceFiniteElement::NodalValuesToDofs().
 *
 * The main use of this function is the interpolation of Dirichlet data on
 * the Dirichlet part of the boundary of a domain.
 *
 * ### Template parameter type requirements
 * - SCALAR must be a type like `complex<double>`
 * - EDGESELECTOR must be compatible with
 *                `std::function<bool(const Entity &)>`
 * - FUNCTION is a scalar valued \ref mesh_function "MeshFunction" which is
 * evaluated on edges
 *
 * This function is meant to supply the information needed for the
 * elimination of Dirichlet boundary conditions by means of the function
 * lf::assemble::FixFlaggedSolutionComponents().
 *
 * ### Example
 * @snippet fe_tools.cc InitEssentialConditionFromFunction
 */
template <typename SCALAR, typename EDGESELECTOR, typename FUNCTION>
std::vector<std::pair<bool, SCALAR>> InitEssentialConditionFromFunction(
    const lf::fe::ScalarFESpace<SCALAR> &fes, EDGESELECTOR &&esscondflag,
    FUNCTION &&g) {
  static_assert(mesh::utils::isMeshFunction<std::remove_reference_t<FUNCTION>>);
  // static_assert(isMeshFunction<FUNCTION>, "g must by a MeshFunction object");

  // *** I: Preprocessing ***

  // Underlying mesh
  const lf::mesh::Mesh &mesh{*fes.Mesh()};

  // Vector for returning flags and fixed values for dofs
  const size_type num_coeffs = fes.LocGlobMap().NumDofs();
  std::vector<std::pair<bool, SCALAR>> flag_val_vec(num_coeffs,
                                                    {false, SCALAR{}});

  // *** II: Local computations ****
  // Visit all edges of the mesh (codim-1 entities)
  for (const lf::mesh::Entity *edge : mesh.Entities(1)) {
    // Check whether the current edge carries dofs to be imposed by the
    // function g. The decision relies on the predicate `esscondflag`
    if (esscondflag(*edge)) {
      // retrieve reference finite element:
      auto fe = fes.ShapeFunctionLayout(*edge);
      LF_ASSERT_MSG(fe != nullptr,
                    "No ScalarReferenceFiniteElement for this edge!");
      auto ref_eval_pts = fe->EvaluationNodes();

      // Evaluate mesh function at several points specified by their
      // reference coordinates.
      auto g_vals = g(*edge, ref_eval_pts);

      // Compute degrees of freedom from function values in evaluation
      // points
      auto dof_vals = fe->NodalValuesToDofs(
          Eigen::Map<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>>(&g_vals[0], 1,
                                                               g_vals.size()));
      LF_ASSERT_MSG(dof_vals.size() == fe->NumRefShapeFunctions(),
                    "Mismatch " << dof_vals.size() << " <-> "
                                << fe->NumRefShapeFunctions());
      LF_ASSERT_MSG(fes.LocGlobMap().NumLocalDofs(*edge) == dof_vals.size(),
                    "Mismatch " << fes.LocGlobMap().NumLocalDofs(*edge)
                                << " <-> " << dof_vals.size());
      // Fetch indices of global shape functions associated with current
      // edge
      auto gdof_indices{fes.LocGlobMap().GlobalDofIndices(*edge)};
      int k = 0;
      // Set flags and values; setting, no accumulation here!
      for (const lf::assemble::gdof_idx_t gdof_idx : gdof_indices) {
        flag_val_vec[gdof_idx] = {true, dof_vals[k++]};
      }
    }
  }
  return flag_val_vec;
}

}  // namespace lf::fe

#endif
