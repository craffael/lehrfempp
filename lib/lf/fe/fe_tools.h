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
#include "fe_space_lagrange_uniform.h"
#include "loc_comp_norms.h"

namespace lf::fe {

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

/**
 * @brief Summation of local contributions over _all_ cells
 *
 * @sa SumCellFEContrib()
 */

template <typename LOC_COMP, typename COEFFVECTOR>
double SumCellFEContrib(const lf::assemble::DofHandler &dofh,
                        LOC_COMP &loc_comp, const COEFFVECTOR &uh) {
  return SumCellFEContrib(dofh, loc_comp, uh, base::PredicateTrue{});
}

/**
 * @brief Computation of an inner product norm of the difference of a finite
 * element function and a general functions.
 *
 * @tparam LOC_NORM_COMP helper type like lf::fe::MeshFunctionL2NormDifference
 * @tparam COEFFVECTOR a vanilla vector type, `std::vector<SCALAR>`
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
  return NormOfDifference(dofh, loc_comp, uh, base::PredicateTrue{});
}

// ******************************************************************************

// Output control for nodal projection
// TODO(ralfh) putting this in a header file leads to multiple symbols with the
// name ctrl_l2
// CONTROLDECLAREINFO(ctrl_prj, "ctrl_prj",
//                    "Output control for  NodalProjection()");
static const unsigned int kout_prj_cell = 1;
static const unsigned int kout_prj_vals = 2;

/*
 * @brief Computes nodal projection of a function and returns the finite
 * element basis expansion coefficients of the result
 *
 * @tparam SCALAR a scalar type
 * @tparam FUNCTOR a \ref mesh_function "MeshFunction" representing the scalar
 * valued function that should be projected
 * @tparam SELECTOR predicate type for selecting cells to be visited
 *
 * @param fe_space a uniform Lagrangian finite element space, providing finite
 *        element specifications for the cells of the mesh
 * @param u functor object supplying a scalar-valued function
 * @param pred predicate object for the selection of relevant cells
 * @return column vector of basis expansion coefficients for the resulting
 *         finite element function
 *
 * The implemetation relies on the method
 * ScalarReferenceFiniteElement::NodalValuesToDofs(). Refer to its
 * documentation. This method is called for each active cell to **set** the
 * coefficients for the global shape functions associated with that cell.
 *
 */
template <typename SCALAR, typename FUNCTOR,
          typename SELECTOR = base::PredicateTrue>
auto NodalProjection(const FeSpaceUniformScalar<SCALAR> &fe_space, FUNCTOR &&u,
                     SELECTOR &&pred = base::PredicateTrue{}) {
  static_assert(isMeshFunction<std::remove_reference_t<FUNCTOR>>);
  // choose scalar type so it can hold the scalar type of u as well as SCALAR
  using scalarMF_t = MeshFunctionReturnType<std::remove_reference_t<FUNCTOR>>;
  using scalar_t = decltype(SCALAR(0) * scalarMF_t(0));
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
    // TODO(ralfh) uncommend when ctrl_prj is well-defined.
    // SWITCHEDSTATEMENT(ctrl_prj, kout_prj_cell,
    //                   std::cout << ref_el << ", shape = \n"
    //                             << geo_ptr->Global(ref_el.NodeCoords())
    //                             << std::endl);

    // Information about local shape functions on reference element
    const ScalarReferenceFiniteElement<double> &ref_shape_fns{
        *fe_space.ShapeFunctionLayout(ref_el)};
    // Number of evaluation nodes
    const size_type num_eval_nodes = ref_shape_fns.NumEvaluationNodes();
    // Obtain reference coordinates for evaluation nodes
    const Eigen::MatrixXd ref_nodes(ref_shape_fns.EvaluationNodes());

    // Collect values of function to be projected in a row vector
    auto uvalvec = u(cell, ref_nodes);

    // Compute the resulting local degrees of freedom
    auto dofvec(ref_shape_fns.NodalValuesToDofs(
        Eigen::Map<Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>>(
            &uvalvec[0], uvalvec.size(), 1)));

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
    const FeSpaceUniformScalar<SCALAR> &fe_space, DIFF_COEFF alpha,
    REACTION_COEFF gamma, TMPMATRIX &A) {
  using scalar_t = typename TMPMATRIX::Scalar;
  // The underlying finite element mesh
  const lf::mesh::Mesh &mesh{*fe_space.Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
  // Object taking care of local computations. No selection of a subset
  // of cells is specified.
  LagrangeFEEllBVPElementMatrix<scalar_t, decltype(alpha), decltype(gamma)>
      elmat_builder(fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()),
                    fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()),
                    alpha, gamma);
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
void SecOrdBVPLagrFEBoundaryGalMat(const FeSpaceUniformScalar<SCALAR> &fe_space,
                                   COEFF eta, EDGESELECTOR edge_sel,
                                   TMPMATRIX &A) {
  using scalar_t = typename TMPMATRIX::Scalar;
  // The underlying finite element mesh
  const lf::mesh::Mesh &mesh{*fe_space.Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

  // Object taking care of local computations.
  LagrangeFEEdgeMassMatrix<scalar_t, decltype(eta), decltype(edge_sel)>
      edgemat_builder(fe_space.ShapeFunctionLayout(lf::base::RefEl::kSegment()),
                      eta, edge_sel);
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
 * This function relies on the the class \ref ScalarFELocalLoadVector for local
 * computations and the function \ref lf::assemble::AssembleVectorLocally()
 * for assembly.
 *
 * @note the functions performs an update of the vector
 */
template <typename SCALAR, typename VECTOR, typename FUNCTOR>
void LagrFEVolumeRightHandSideVector(
    const FeSpaceUniformScalar<SCALAR> &fe_space, FUNCTOR f, VECTOR &phi) {
  using scalar_t = typename VECTOR::value_type;
  // The underlying finite element mesh
  const lf::mesh::Mesh &mesh{*fe_space.Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
  // Object taking care of local computations. No selection of a subset
  // of cells is specified.
  ScalarFELocalLoadVector<scalar_t, FUNCTOR> elvec_builder(
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()),
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()), f);
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
 * This function relies on the the class \ref ScalarFEEdgeLocalLoadVector
 * for local computations and the function \ref
 * lf::assemble::AssembleVectorLocally() for assembly.
 *
 * @note the functions performs an update of the vector
 */
template <typename SCALAR, typename VECTOR, typename FUNCTOR,
          typename EDGESELECTOR>
void LagrFEBoundaryRightHandSideVector(
    const FeSpaceUniformScalar<SCALAR> &fe_space, FUNCTOR data,
    EDGESELECTOR edge_sel, VECTOR &phi) {
  using scalar_t = typename VECTOR::value_type;
  // The underlying finite element mesh
  const lf::mesh::Mesh &mesh{*fe_space.Mesh()};
  // The local-to-global index map for the finite element space
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
  // Object taking care of local computations. No selection of a subset
  // of cells is specified.
  ScalarFEEdgeLocalLoadVector<scalar_t, FUNCTOR, EDGESELECTOR> elvec_builder(
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kSegment()), data,
      edge_sel);
  // Invoke assembly on edges (codim == 1), update vector
  AssembleVectorLocally(1, dofh, elvec_builder, phi);
}

/**
 * @brief Initialization of flags/values for dofs of a uniform Lagrangian
 *        finite element space whose values are imposed by a specified function.
 *
 * @tparam SCALAR scalar type for BVP = return type of the function g
 * @tparam EDGESELECTOR predicate returning true for edges with fixed dofs
 * @tparam FUNCTION functor type for object providing scalar-valued function
 *
 * @param dofh local-to-global index mapping
 * @param fe_spec_edge description of arrangement for local shape functions
 *        for a `kSegment`-type entity
 * @param esscondflag predicate object whose evaluation operator returns
 *        true for all edges whose associated degrees of freedom should be
 *        set a fixed value.
 * @return a vector of flag-value pairs, a `true` first component indicating
 *         a fixed dof, with the second component providing the value in this
 *         case.
 *
 * This function interpolates a scalar-valued function into a Lagrangian
 * finite element space restricted to a set of edges. It relies on the
 * method ScalarReferenceFiniteElement::NodalValuesToDofs().
 *
 * The main use of this function is the interpolation of Dirichet data on the
 * Dirichlet part of vthe boundary of a domain.
 *
 * ### Template parameter type requirements
 * - SCALAR must be a type like `complex<double>`
 * - EDGESELECTOR must be compatible with
 *                `std::function<bool(const Entity &)>`
 * - FUNCTION must be a type like 'std::function<SCALAR(VectorXd)>`
 *
 * This function is meant to supply the information needed for the elimination
 * of Dirichlet boundary conditions by means of the function
 * lf::assemble::fix_flagged_solution_components().
 */
template <typename SCALAR, typename EDGESELECTOR, typename FUNCTION>
std::vector<std::pair<bool, SCALAR>> InitEssentialConditionFromFunction(
    const lf::assemble::DofHandler &dofh,
    const ScalarReferenceFiniteElement<SCALAR> &fe_spec_edge,
    EDGESELECTOR &&esscondflag, FUNCTION &&g) {
  LF_ASSERT_MSG(fe_spec_edge.RefEl() == lf::base::RefEl::kSegment(),
                "finite element specification must be for an edge!");

  // *** I: Preprocessing ***
  // Fetch numbers of evaluation nodes and local shape functions
  const size_type num_eval_pts = fe_spec_edge.NumEvaluationNodes();
  const size_type num_rsf = fe_spec_edge.NumRefShapeFunctions();

  // Preprocessing: obtain evaluation nodes on reference segment [0,1]
  const Eigen::MatrixXd ref_eval_pts{fe_spec_edge.EvaluationNodes()};
  Eigen::MatrixXd eval_pts(2, num_eval_pts);
  Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dof_vals(num_eval_pts);

  // Underlying mesh
  const lf::mesh::Mesh &mesh{*dofh.Mesh()};

  // Vector for returning flags and fixed values for dofs
  const size_type num_coeffs = dofh.NoDofs();
  std::vector<std::pair<bool, SCALAR>> flag_val_vec(num_coeffs,
                                                    {false, SCALAR{}});

  // *** II: Local computations ****
  // Visit all edges of the mesh (codim-1 entities)
  for (const lf::mesh::Entity &edge : mesh.Entities(1)) {
    // Check whether the current edge carries dofs to be imposed by the
    // function g.
    if (esscondflag(edge) == true) {
      // Fetch the shape of the edge
      const lf::geometry::Geometry *edge_geo_p{edge.Geometry()};
      auto g_vals = g(edge, ref_eval_pts);

      // Compute degrees of freedom from function values in evaluation points
      dof_vals = fe_spec_edge.NodalValuesToDofs(
          Eigen::Map<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>>(&g_vals[0], 1,
                                                               g_vals.size()));
      LF_ASSERT_MSG(dof_vals.size() == num_rsf,
                    "Mismatch " << dof_vals.size() << " <-> " << num_rsf);
      LF_ASSERT_MSG(
          dofh.NoLocalDofs(edge) == num_rsf,
          "Mismatch " << dofh.NoLocalDofs(edge) << " <-> " << num_rsf);
      // Fetch indices of global shape functions associated with current edge
      auto gdof_indices{dofh.GlobalDofIndices(edge)};
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
