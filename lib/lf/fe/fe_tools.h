#ifndef LF_FETOOLS_H
#define LF_FETOOLS_H
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Functions acting on finite element spaces
 * @author Ralf Hiptmair
 * @date November 2018
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/quad/quad.h>
#include "fe_space.h"
#include "loc_comp_norms.h"

namespace lf::fe {

CONTROLDECLAREINFO(ctrl_l2, "ctrl_l2", "Output control for SumCellFEContrib()");
static const unsigned int kout_l2_qr = 1;
static const unsigned int kout_l2_rsfvals = 2;

/**
 * @brief Summation of cell values computed based on a finite element
 *        function specified through its basis expansion coefficients
 *
 * @tparam LOC_COMP type taking care of local computations
 * @tparam COEFFVECTOR a vanilla vector type
 * @tparam SELECTOR a predicate for selecting cells
 *
 * ### type requirements
 *
 * - The type LOC_COMP must feature an `isActive()` method for the
 *   selection of cells to be visited, must provide a type `dofvector_t`
 *   for passing coefficients of local shape functions and a method
 * ~~~
 * double operator ()(const lf::mesh::Entity &cell, const DOFVECTOR &dofs);
 * ~~~
 *   that performs the local computations.
 * - The type COEFFVECTOR must provide component access through `[]`,
 *   a `size()` method telling the vector length, and `value_type` typedef
 *   telling the  component type.
 * - The type functor must provide an evaluation operator `()` taking
 *   a point coordinate vector `Eigen::Vector2d` as an argument.
 *
 * @param dofh Local-to-global index mapping belonging to a finite element space
 * @param loc_comp reference to helper object for local computations
 *        This object must be aware of the shape functions!
 * @param uh coefficient vector of finite element function
 *
 */
template <typename LOC_COMP, typename COEFFVECTOR, typename SELECTOR>
double SumCellFEContrib(const lf::assemble::DofHandler &dofh,
                        LOC_COMP &loc_comp, const COEFFVECTOR &uh,
                        SELECTOR &&pred) {
  // Type for passing local coefficient vectors
  using dofvector_t = std::vector<typename COEFFVECTOR::value_type>;
  // Retrieve the mesh
  const lf::mesh::Mesh &mesh{*dofh.Mesh()};

  // Check whether sufficiently large vector uh
  LF_ASSERT_MSG(uh.size() >= dofh.NoDofs(), "uh vector too short!");

  double sum = 0.0;
  /// loop over cells ( = codim-0 entities)
  for (const lf::mesh::Entity &cell : mesh.Entities(0)) {
    if (pred(cell)) {
      // Number of local shape functions
      const size_type num_loc_dofs = dofh.NoLocalDofs(cell);
      // Allocate a temporary vector of appropriate size
      dofvector_t loc_coeffs(num_loc_dofs);
      // Fetch global numbers of local shape functions
      lf::base::RandomAccessRange<const gdof_idx_t> ldof_gidx(
          dofh.GlobalDofIndices(cell));
      // Copy degrees of freedom into temporory vector
      for (int j = 0; j < num_loc_dofs; ++j) {
        loc_coeffs[j] = uh[ldof_gidx[j]];
      }
      // Sum local contribution
      sum += loc_comp(cell, loc_coeffs);
    }  // end if (isACtive())
  }
  return sum;
}

/** @brief predicate for non-discriminatory positive selection of entities */
struct DefaultEntitySelector {
  bool operator()(const lf::mesh::Entity & /*cell*/) { return true; }
};

/**
 * @brief Summation of local contributions over _all_ cells
 *
 * @sa SumCellFEContrib(const lf::assemble::DofHandler,LOC_COMP, const
 * COEFFVECTOR &uh, SELECTOR)
 */

template <typename LOC_COMP, typename COEFFVECTOR>
double SumCellFEContrib(const lf::assemble::DofHandler &dofh,
                        LOC_COMP &loc_comp, const COEFFVECTOR &uh) {
  return SumCellFEContrib(dofh, loc_comp, uh, DefaultEntitySelector());
}

/**
 * @brief Computation of an inner product norm of the difference of a finite
 * element function and a general functions.
 *
 * @tparam LOC_NORM_COMP helper type like lf::fe::LocalL2NormDifference
 * @tparam COEFFVECTOR a vanilla vector type
 * @tparam SELECTOR a predicate for selecting cells
 *
 * ### type requirements
 *
 * - The type LOC_NORM_COMP must feature an `isActive()` method for the
 *   selection of cells to be visited, must provide a type `dofvector_t`
 *   for passing coefficients of local shape functions and a method
 * ~~~
 * double operator () (const lf::mesh::Entity &cell, const DOFVECTOR &dofs);
 * ~~~
 *   that performs the local computations and returns the **square** of
 *   the local norm of the difference function
 * - The type COEFFVECTOR must provide component access through `[]`,
 *   a `size()` method telling the vector length, and `value_type` typedef
 *   telling the  component type.
 * - The type functor must provide an evaluation operator `()` taking
 *   a point coordinate vector `Eigen::Vector2d` as an argument.
 *
 * @param dofh Local-to-global index mapping belonging to a finite element space
 * @param loc_comp reference to helper object for local computations
 *        This object must be aware of the shape functions!
 * @param uh coefficient vector of finite element function
 *
 */
template <typename LOC_NORM_COMP, typename COEFFVECTOR, typename SELECTOR>
double NormOfDifference(const lf::assemble::DofHandler &dofh,
                        LOC_NORM_COMP &loc_comp, const COEFFVECTOR &uh,
                        SELECTOR &&pred) {
  const double norm_sq = SumCellFEContrib(dofh, loc_comp, uh, pred);
  return std::sqrt(norm_sq);
}
template <typename LOC_NORM_COMP, typename COEFFVECTOR>
double NormOfDifference(const lf::assemble::DofHandler &dofh,
                        LOC_NORM_COMP &loc_comp, const COEFFVECTOR &uh) {
  return NormOfDifference(dofh, loc_comp, uh, DefaultEntitySelector());
}

// ******************************************************************************

// Output control for nodal projection
CONTROLDECLAREINFO(ctrl_prj, "ctrl_prj",
                   "Output control for  NodalProjection()");
static const unsigned int kout_prj_cell = 1;
static const unsigned int kout_prj_vals = 2;

/*
 * @brief Computes nodal projection of a function and returns the finite
 * element basis expansion coefficients of the result
 *
 */
template <typename FUNCTOR, typename SELECTOR>
Eigen::Matrix<typename std::invoke_result<FUNCTOR, Eigen::VectorXd>::type,
              Eigen::Dynamic, 1>
NodalProjection(const UniformScalarFiniteElementSpace &fe_space, FUNCTOR u,
                SELECTOR &&pred) {
  // Scalar type determined by the function
  using scalar_t = typename std::invoke_result<FUNCTOR, Eigen::VectorXd>::type;
  // Return type, type for FE coefficient vector
  using dof_vec_t = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

  // Underlying mesh instance
  const lf::mesh::Mesh &mesh{*fe_space.Mesh()};
  // Fetch local-to-global index mapping for shape functions
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

  // Vector for returning expansion coefficients
  dof_vec_t glob_dofvec(dofh.NoDofs());
  glob_dofvec.setZero();

  // Loop over all cells
  for (const lf::mesh::Entity &cell : mesh.Entities(0)) {
    // Topological type of the cell
    const lf::base::RefEl ref_el{cell.RefEl()};
    // Query the shape of the cell
    const lf::geometry::Geometry *geo_ptr = cell.Geometry();
    LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
    LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                  "Only 2D implementation available!");
    SWITCHEDSTATEMENT(ctrl_prj, kout_prj_cell,
                      std::cout << ref_el << ", shape = \n"
                                << geo_ptr->Global(ref_el.NodeCoords())
                                << std::endl);

    // Information about local shape functions on reference element
    const ScalarReferenceFiniteElement<double> &ref_shape_fns{
        *fe_space.ShapeFunctionLayout(ref_el)};
    // Number of evaluation nodes
    const size_type num_eval_nodes = ref_shape_fns.NumEvaluationNodes();
    // Obtain reference coordinates for evaluation nodes
    const Eigen::MatrixXd ref_nodes(ref_shape_fns.EvaluationNodes());
    // Obtain world coordinates of evaluation points
    const Eigen::MatrixXd mapped_evlpts(geo_ptr->Global(ref_nodes));
    // Collect values of function to be projected in a row vector
    Eigen::Matrix<scalar_t, 1, Eigen::Dynamic> uvalvec(num_eval_nodes);
    for (int j = 0; j < num_eval_nodes; ++j) {
      uvalvec[j] = u(mapped_evlpts.col(j));
    }
    // Compute the resulting local degrees of freedom
    auto dofvec(ref_shape_fns.NodalValuesToDofs(uvalvec));

    // Set the corresponing global degrees of freedom
    // Note: "Setting", not "adding to"
    // Number of local shape functions
    const size_type num_loc_dofs = dofh.NoLocalDofs(cell);
    LF_ASSERT_MSG(
        dofvec.size() == num_loc_dofs,
        "Size mismatch: " << dofvec.size() << " <-> " << num_loc_dofs);
    // Fetch global numbers of local shape functions
    lf::base::RandomAccessRange<const gdof_idx_t> ldof_gidx(
        dofh.GlobalDofIndices(cell));
    // Insert dof values into the global coefficient vector
    for (int j = 0; j < num_loc_dofs; ++j) {
      glob_dofvec[ldof_gidx[j]] = dofvec[j];
    }
  }
  return glob_dofvec;
}

/** @sa NodalProjection(const UniformScalarFiniteElementSpace &fe_space, FUNCTOR
 * u, SELECTOR &&pred)
 *
 * Nodal projection applied to all cells
 */
template <typename FUNCTOR>
Eigen::Matrix<typename std::invoke_result<FUNCTOR, Eigen::VectorXd>::type,
              Eigen::Dynamic, 1>
NodalProjection(const UniformScalarFiniteElementSpace &fe_space, FUNCTOR u) {
  return NodalProjection(fe_space, u, DefaultEntitySelector());
}

/** @brief Incremental assembly of global finite element Galerkin matrix
 *
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
 * matrix for a second-order elliptic boundary value problem.
 */
template <typename TMPMATRIX, typename DIFF_COEFF, typename REACTION_COEFF>
void SecOrdBVPLagrFEFullInteriorGalMat(
    const UniformScalarFiniteElementSpace &fe_space, DIFF_COEFF alpha,
    REACTION_COEFF gamma, TMPMATRIX &A) {
  // The underlying finite element mesh
  const lf::mesh::Mesh &mesh{*fe_space.Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
  // Object taking care of local computations. No selection of a subset
  // of cells is specified.
  LagrangeFEEllBVPElementMatrix<double, decltype(alpha), decltype(gamma)>
      elmat_builder(fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()),
                    fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()),
                    alpha, gamma);
  // Invoke assembly on cells
  AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
}

}  // namespace lf::fe

#endif
