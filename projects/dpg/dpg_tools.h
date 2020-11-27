#ifndef PROJECTS_DPG_DPG_TOOLS_H
#define PROJECTS_DPG_DPG_TOOLS_H

/**
  @file
 * @headerfile projects/dpg/dpg_tools.h
 * @brief contains helper functions and classes used in dpg computations
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <algorithm>
#include <cmath>
#include <vector>

#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>

#include "product_dofhandler.h"
#include "product_element_matrix_provider.h"
#include "product_element_vector_provider.h"
#include "product_fe_space.h"

namespace projects::dpg {

/**
 * @brief Constructs a quadrature rule  on the boundary of a reference element
 * @param ref_el The reference element on whose boundary  a quadrule
 * should be constructed (triangle or quadrilateral).
 * @param qr a quadrule on the reference line segment.
 * @return a vector x s.t. x[i] is a valid quadrule on edge i of the referene
 * element ref_el
 *
 * This function constructs a quadrature rule on the boundary of a reference
 * element based on another quadrature rule on the reference line segment which
 * is transformed to all edges of the reference element.
 *
 */
std::vector<lf::quad::QuadRule> BoundaryQuadRule(lf::base::RefEl ref_el,
                                                 const lf::quad::QuadRule& qr);

/**
 * @brief Compute the outer normals of a geometry object
 * @param geometry The  geometry of the cell for which the outer normals are
 * computed. Should be either a triangle or quadrilateral
 * @return The outer normals, packed into the columns of a 2xn matrix,
 *         where n is the number of edges of the  geometry .
 */
Eigen::MatrixXd OuterNormals(const lf::geometry::Geometry& geometry);

/**
 * @brief Class which allows evaluation of the \f$ \text{sgn}_K \f$ function.
 *
 * DPG formulations can
 * involve the notion of prescribed edge normals \f$ n_e \f$. These can be
represented by a  unit
 * vector field on the mesh skeleton that is normal to the edges and points
outwards on the boundary.
 * The prescribed edge normals are fully specified by the definition of  a
function \f$ \text{sgn}_K: K \rightarrow \{ -1,1 \} \f$
 * that distinguishes them from the outer normals \f$n_K\f$ on each cell \f$ K
\f$ of the mesh.
 *
 * \f[
 \text{sgn}_K = \begin{cases}
    1  & \text{if }  n_K    =  n_e \\
    -1 & \text{if }  n_K    = -  n_e
\end{cases} \f]
 * This class provides an implementation of the \f$ \text{sgn}_K \f$ function
that further guarantues that
 * the function is constant on any edge \f$ e \f$  of \f$ \partial K \f$. In
particular the function is specified
 * fully on any cell \f$ K \f$ by evaluating it in all midpoints of the edges of
\f$ K \f$.
 */
class PrescribedSignProvider {
 public:
  /** default constructors and destructors:
   * @note creates an invalid object that cannot be used.*/
  PrescribedSignProvider() = default;
  PrescribedSignProvider(const PrescribedSignProvider&) = delete;
  PrescribedSignProvider(PrescribedSignProvider&&) noexcept = default;
  PrescribedSignProvider& operator=(const PrescribedSignProvider&) = delete;
  PrescribedSignProvider& operator=(PrescribedSignProvider&&) noexcept =
      default;
  virtual ~PrescribedSignProvider() = default;

  /**
   * @brief main constructor, initializes the \f$\text{sgn}_K\f$ function for
   * all cells \f$ K \f$  of the passed mesh
   * @param mesh_ptr: underlying mesh
   */
  explicit PrescribedSignProvider(
      std::shared_ptr<const lf::mesh::Mesh> mesh_ptr)
      : mesh_ptr_(std::move(mesh_ptr)),
        maxElement_ptr_(std::move(
            lf::mesh::utils::make_CodimMeshDataSet(mesh_ptr_, 1, -1))) {
    init();
  }

  /**
   * @brief evaluates the the \f$ \text{sgn}_K \f$ function at the midpoint of
   * the edge \f$ e \f$ of cell \f$ K \f$
   * @param element the cell \f$ K \f$ on which the  function is evauated
   * @param edge the edge \f$ e \subset \partial K \f$  for which the function
   * is evaluated.
   * @return value of \f$ \text{sgn}_K \f$ at the midpoint of \f$ e \f$.
   *
   * @note By construction the value of the function is constant along an edge
   * \f$ e \subset \partial K \f$  and this is the value  function on the
   * complete edge \f$ e \f$
   *
   */
  [[nodiscard]] int PrescribedSign(const lf::mesh::Entity& element,
                                   const lf::mesh::Entity& edge) const;

 private:
  /** Initialization of the maxElement data set*/
  void init();

  /** the underlying mesh */
  std::shared_ptr<const lf::mesh::Mesh> mesh_ptr_;

  /** a mesh dataset for edges. maxElement_ptr_(edge) equals the maximum entity
   *index of a cell whose subentity edge is. We define the prescribed normal on
   *that edge to be the outer normal of that cell with the bigger index.
   * @note An edge can only be a subentity of maximum two cells, so this
   *uniquely determines the sign.
   **/
  std::shared_ptr<lf::mesh::utils::CodimMeshDataSet<int>> maxElement_ptr_;
};

/**
 * @brief flag all edges located on the inflow boundary
 * @param mesh_p the underlying mesh
 * @param beta the prescribed convection field \f$ \beta \f$ .
 * @return a datastructure of boolean values which specifies for each edge, if
 * it lies on the inflow boundary or not.
 *
 * An edge is located on the inflow boundary, if it is located on the boundary
 * \f$ \partial \Omega \f$ and if it holds that \f$ \beta \cdot n_{\Omega} \f$
 * < 0 on that edge.
 *
 * @note For non constant convection fields this condition may or may not hold
 * on the whole edge, so we only check this condition in the _midpoint_ of the
 * edge.
 */

template <typename CONVECTION_COEFF>
lf::mesh::utils::CodimMeshDataSet<bool> flagEntitiesOnInflowBoundary(
    const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
    CONVECTION_COEFF beta) {
  // flag edges on the boundary
  auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1);

  // initialize inflow boundary flags to fals
  lf::mesh::utils::CodimMeshDataSet<bool> inflow_boundary{mesh_p, 1, false};

  // iterate over all cells:
  for (const auto* const cell : mesh_p->Entities(0)) {
    // query geometry
    auto geo_ptr = cell->Geometry();
    auto sub_entities = cell->SubEntities(1);

    // calcualte normals
    Eigen::MatrixXd normals = OuterNormals(*geo_ptr);

    // iterate over all edges of the cell
    for (int i = 0; i < cell->RefEl().NumSubEntities(1); ++i) {
      const auto* const edge = sub_entities[i];

      // check the inflow condition at the center of the edge
      bool inflow = normals.col(i).transpose() *
                        beta(*(sub_entities[i]),
                             (Eigen::MatrixXd(1, 1) << 0.5).finished())[0] <
                    0.0;
      inflow_boundary(*edge) = bd_flags(*edge) && inflow;
    }
  }
  return inflow_boundary;
}

/**
 * @brief flag all edges located on the outflow boundary
 * @param mesh_p the underlying mesh
 * @param beta the prescribed convection field  \f$ \beta \f$
 * @return a datastructure of boolean values, which specifies for each edge, if
 * it lies on the outflow boundary or not
 *
 * An edge is located on the outflow boundary, if it is located on the boundary
 * \f$ \partial \Omega \f$ and if it holds that \f$ \beta \cdot n_{\Omega}  \geq
 * 0  \f$ on that edge.
 *
 * @note For non constant convection fields this condition may or may not hold
 * on the whole edge, so we only check this condition in the _midpoint_ of the
 * edge.
 */
template <typename CONVECTION_COEFF>
lf::mesh::utils::CodimMeshDataSet<bool> flagEntitiesOnOutflowBoundary(
    const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
    CONVECTION_COEFF beta) {
  // flag all entities on the boundary:
  auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1);

  // flag all entities on the inflow boundary:
  auto inflow_boundary = flagEntitiesOnInflowBoundary(mesh_p, beta);

  // initialize outflow flags:
  lf::mesh::utils::CodimMeshDataSet<bool> outflow_boundary{mesh_p, 1, false};

  // outflow boundary dofs are those, that lie on the boundary, but not on the
  // inflow boundary This guarantues a partition of the boundary edges into
  // inflow and outflow edges.
  for (const auto* const edge : mesh_p->Entities(1)) {
    outflow_boundary(*edge) = bd_flags(*edge) && (!inflow_boundary(*edge));
  }
  return outflow_boundary;
}

/**
 * @brief Initialization of flags/values for dofs of fluxes and traces
 * of a on a product finite element space  whose values are imposed by
 * two specified functions on two specified parts of the boundary
 *
 * @tparam SCALAR scalar type like 'double'
 * @tparam EDGESELECTOR_DIRICHLET and EDGESELECTOR_NEUMANN predicates returning
 * true for edges with fixed dofs
 * @tparam FUNCTION_G and FUNCTION_H \ref mesh_function "MeshFunction" which
 * define the imposed values on the edges
 *
 * @brief InitEssentialConditionsFromFunctions
 * @param dofh the local_to-global index mapping of the product space
 * @param fe_spec_edge_trace description of the  local shape functions for a
 * 'kSegment' type entity for the trace component
 * @param fe_spec_edge_flux description of the local shape functions for a
 * 'kSegment' type entity for the flux component
 * @param dirichletcondflag predicate object whose evaluation operator returns
 * true for all edges on which the degrees of freedom for the trace component
 * should be set to a fixed value
 * @param neumanncondflag predicate object whose evaluation operator returns
 * true for all edges on which the degrees of freedom for the flux component
 * should be set to a fixed value.
 * @param g functor providing the evaluation of the essential conditions for the
 * trace component (dirichlet data)
 * @param h functor providing the evalutation of the essential conditions for
 * the flux component (neumann data)
 * @param trace_component the index of the trace component in the product space
 * @param flux_component the index of the flux component in the product space
 * @return  a vector of flag-value pairs, a 'true' first component indicates a
 * fixed dof, the second component provides the value in this case.
 *
 * The dpg formulation usually involves variables, which represent fluxes or
 * traces. This allows us to impose both Dirichlet and Neumann boundary
 * conditions as essential boundary conditions on the corresponding components
 * and is a generalization of the lehrfempp function
 * lf::uscalfe::InitEssentialConditionFromFunction to product spaces.
 *
 *
 *  This function is meant to supply information for the
 *  elimination of essential dofs by means of the function
 * lf::assembl::fix_flagged_solution_components().
 *
 */
template <typename SCALAR, typename EDGESELECTOR_DIRICHLET,
          typename EDGESELECTOR_NEUMANN, typename FUNCTION_G,
          typename FUNCTION_H>
std::vector<std::pair<bool, SCALAR>> InitEssentialConditionsFromFunctions(
    const ProductUniformFEDofHandler& dofh,
    const lf::fe::ScalarReferenceFiniteElement<SCALAR>& fe_spec_edge_trace,
    const lf::fe::ScalarReferenceFiniteElement<SCALAR>& fe_spec_edge_flux,
    EDGESELECTOR_DIRICHLET&& dirichletcondflag,
    EDGESELECTOR_NEUMANN&& neumanncondflag, FUNCTION_G&& g, FUNCTION_H&& h,
    size_type trace_component, size_type flux_component) {
  // check, that the two components actually differ (necessary?)
  LF_ASSERT_MSG(trace_component != flux_component,
                "Boundary conditions imposed on same component");

  // retrive the dofhandlers for the corresponding  trace and flux components
  const auto& trace_dofh = dofh.ComponentDofHandler(trace_component);
  const auto& flux_dofh = dofh.ComponentDofHandler(flux_component);

  // flag dofs on the component spaces
  std::vector<std::pair<bool, SCALAR>> componentTraceConditions =
      lf::fe::InitEssentialConditionFromFunction(trace_dofh, fe_spec_edge_trace,
                                                 dirichletcondflag, g);
  std::vector<std::pair<bool, SCALAR>> componentFluxConditions =
      lf::fe::InitEssentialConditionFromFunction(flux_dofh, fe_spec_edge_flux,
                                                 neumanncondflag, h);

  // calculate offsets to transform flags on component spaces
  // to the full product space
  size_type trace_offset = dofh.Offset(trace_component);
  size_type flux_offset = dofh.Offset(flux_component);

  std::vector<std::pair<bool, SCALAR>> EssentialConditions(dofh.NumDofs(),
                                                           {false, SCALAR{}});

  // flag dofs on product space corresponding to fixed trace dofs
  for (size_type dof = 0; dof < trace_dofh.NumDofs(); dof++) {
    EssentialConditions[trace_offset + dof] = componentTraceConditions[dof];
  }

  // flag dofs on product space corresponding to fixed flux dofs
  for (size_type dof = 0; dof < flux_dofh.NumDofs(); dof++) {
    EssentialConditions[flux_offset + dof] = componentFluxConditions[dof];
  }

  return EssentialConditions;
}

/**
 * @brief Evaluation of the local DPG error estimators on a mesh
 * @param trial_dofh local-to-global index mapping for the trial space of
 * the DPG method
 * @param stiffness_provider evaluates the extended element stiffness matrix \f$
 * B \f$
 * @param gramian_provider evaluates the local Gramian matrix \f$ G \f$
 * @param load_provider evaluates the extended element load vector \f$ l \f$
 * @param sol_vec finite element vector of a function in the trial space
 * @return A datastructure  of scalar values, which for each cell contains the
 * loacl DPG error estimate for the function represented by sol_vec.
 *
 * @tparam SCALAR a scalar type like 'double'
 */
template <typename SCALAR>
lf::mesh::utils::CodimMeshDataSet<SCALAR> ElementErrorEstimators(
    const ProductUniformFEDofHandler& trial_dofh,
    // const ProductUniformFEDofHandler& test_dofh,
    std::shared_ptr<ProductElementMatrixProvider<SCALAR>> stiffness_provider,
    std::shared_ptr<ProductElementMatrixProvider<SCALAR>> gramian_provider,
    std::shared_ptr<ProductElementVectorProvider<SCALAR>> load_provider,
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>& sol_vec) {
  // retrive the underlying mesh:
  auto mesh_p = trial_dofh.Mesh();

  // retrive the number of cells, initialize the return vector
  lf::mesh::utils::CodimMeshDataSet<SCALAR> element_errors{mesh_p, 0, 0.0};

  // compute the element wise errors
  for (const auto* const cell : mesh_p->Entities(0)) {
    // evaluate the providers on the current cell
    auto G = gramian_provider->Eval(*cell);
    auto B = stiffness_provider->Eval(*cell);
    auto l = load_provider->Eval(*cell);

    // retrive local solution
    Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> local_solution(
        trial_dofh.NumLocalDofs(*cell));
    auto global_dofs = trial_dofh.GlobalDofIndices(*cell);
    for (size_type i = 0; i < trial_dofh.NumLocalDofs(*cell); ++i) {
      local_solution(i) = sol_vec(global_dofs[i]);
    }

    // evaluate posterior error on current element:
    Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> r = l - B * local_solution;

    SCALAR e = std::sqrt(r.transpose() * G.ldlt().solve(r));
    element_errors(*cell) = e;
  }
  return element_errors;
}

/**
 * @brief Evaluation of the DPG error estimator
 * @param mesh_p the underlying mesh
 * @param local_errors data structure containing the local DPG error estimator
 * for each cell of the mesh
 * @return The global DPG error estimator.
 *
 * @tparam Scalar a scalar type like 'double'
 */
template <typename SCALAR>
SCALAR EvalPosteriorError(
    const std::shared_ptr<const lf::mesh::Mesh>& mesh_p,
    const lf::mesh::utils::CodimMeshDataSet<SCALAR>& local_errors) {
  SCALAR error = 0.0;
  for (const auto* const cell : mesh_p->Entities(0)) {
    error += local_errors(*cell) * local_errors(*cell);
  }
  return std::sqrt(error);
}

}  // namespace projects::dpg
#endif  // PROJECTS_DPG_DPG_TOOLS_H
