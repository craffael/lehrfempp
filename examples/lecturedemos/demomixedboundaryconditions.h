/** @file
 * @brief Solving a second-order elliptic BVP with special mixed
 * boundary conditions
 * @author Ralf Hiptmair
 * @date July 2020
 * @copyright MIT License
 */

#ifndef DMXBC_H_
#define DMXBC_H_

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <array>
#include <filesystem>
#include <iostream>
#include <map>
#include <vector>

namespace dmxbc {

/** @brief Computation of gradients of barycentric coordinate functions,
 * exterior edge-length-weighted normals , and of the area of a triangle.
 *
 * @param vertices 2x3 matrix whose columns contain the vertex coordinates of
 * the triangle
 * @return tuple containing the gradients of barycentric coordinate functions,
 * exterior edge-length-weighted normals , and of the area of a triangle. The
 * gradients and the normals are stored in the columns of a 2x3 matrix.
 */
std::tuple<Eigen::Matrix<double, 2, 3>, Eigen::Matrix<double, 2, 3>, double>
getTriangleGradLambdaNormals(Eigen::Matrix<double, 2, 3> vertices);

/** @brief Compute exterior edge-length-weighted normals for triangular or
 * quadrilateral cells
 *
 * @param corners 2xn-matrix whose columns contain the vertex coordinates of the
 * n vertices
 */
Eigen::MatrixXd exteriorCellNormals(const Eigen::MatrixXd& corners);

/** @brief Compute exterior edge-length-weighted normals for edges on the
 * boundary
 *
 * @param mesh_p pointer to finite element mesh containing only triangles with
 * straight edges.
 * @return Edge-indeed array containing edge-length-weighted normals for every
 * edge on the boundary , the zero vector for internal edges
 */
lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> exteriorEdgeWeightedNormals(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p);

/** @brief Reading mesh from Gmsh file and flagging edges on contact boundary
 *         parts
 *
 * @param filename path to Gmsh mesh file (.msh-file)
 * @return pointer to read mesh and an array of integers index by the edges of
 * the mesh
 *
 * The .msh file must contain boundary edges with the two physical names
 * "Contact0" and "Contact1". In the returned array edges tagged with "Context0"
 * recive the id 0, edges belonging to the group "Contact1" get id 1. All other
 * edges are assigned -1.
 */
std::pair<std::shared_ptr<const lf::mesh::Mesh>,
          lf::mesh::utils::CodimMeshDataSet<int>>
readMeshWithTags(std::string filename);

/** @brief Spreading positive edge ids to endpoint nodes
 *
 * @param mesh_p pointer to underlying FE mesh
 * @param edgeids ids for the edges of the mesh
 * @return Array index by the nodes of the mesh containing their indices
 *
 * If an edge carries a positive id, then this id is set for both its endpoints.
 */
lf::mesh::utils::CodimMeshDataSet<int> tagNodes(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p,
    lf::mesh::utils::CodimMeshDataSet<int> edgeids);

/** @brief Solving mixed BVP based on lowest-order Lagrangian finite elements
 *
 * @tparam SIGMAFUNCTOR type for functor providing diffusion coefficient.
 * @param fe_space pointer to finite element space object
 * @param nodeflags array of consecutive ids 0,...,m for nodes of the mesh
 * @param voltvals array of values to which the finite element solution should
 * be set in the nodes with the corresponding id, which serves as an index into
 * the array.
 * @sigma diffusion coefficient
 *
 * Solution of a second-order elliptic scalar boundary value problem with fixed
 * solution values in particular nodes. Standard LehrFEM++ implementation.
 */
template <typename SIGMAFUNCTOR>
Eigen::VectorXd solveMixedBVP(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    lf::mesh::utils::CodimMeshDataSet<int>& nodeflags,
    std::vector<double>&& voltvals, SIGMAFUNCTOR&& sigma) {
  // Mesh functions for coefficients
  lf::mesh::utils::MeshFunctionGlobal mf_sigma{sigma};
  // Reference to current mesh
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space = number of nodes of the mesh
  const lf::base::size_type N_dofs(dofh.NumDofs());
  LF_ASSERT_MSG(N_dofs == mesh.NumEntities(2),
                " N_dofs must agree with number of nodes");
  std::cout << "Solving BVP with " << N_dofs << " FE dofs" << std::endl;
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Assembly of Galerkin matrix for full finite element space
  lf::fe::DiffusionElementMatrixProvider<double, decltype(mf_sigma)>
      elmat_builder(fe_space, mf_sigma);
  // Cell-oriented assembly
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  // Right-hand side vector is zero
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  // Number of Dirichlet boundary parts
  const unsigned int NContacts = voltvals.size();
  // Selector functor for Dirichlet data
  auto selector = [&](lf::assemble::gdof_idx_t idx) -> std::pair<bool, double> {
    const lf::mesh::Entity& node{dofh.Entity(idx)};
    const int ids = nodeflags(node);
    if ((ids >= 0) && (ids < NContacts)) {
      return {true, voltvals[ids]};
    }
    return {false, 42.0};
  };
  // Eliminate Dirichlet dofs from linear system
  lf::assemble::FixFlaggedSolutionComponents<double>(selector, A, phi);
  // Assembly completed: Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Solve linear system using Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);  // LU decomposition
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  Eigen::VectorXd sol_vec = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");
  return sol_vec;
}  // end solveMixedBVP

/** @brief Evaluation of boundary formula for contact flux
 *
 * @tparam SIGMAFUNCTION functor type for diffusion coefficient
 * @param fe_space pointer to lowest-order FE space object
 * @param sigma diffusion coefficient
 * @param edgeids edge-indexed array of id numbers marking contacts
 * @param contact_id
 *
 * Direct computation of flux through a contact based on a finite-element
 * solution of the mixed BVP.
 */
template <typename SIGMAFUNCTION>
double contactFlux(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& sol_vec, SIGMAFUNCTION&& sigma,
    const lf::mesh::utils::CodimMeshDataSet<int>& edgeids, int contact_id = 0) {
  // Obtain object managing dof indexing
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Check whether coefficient vector matches dof handler
  LF_ASSERT_MSG(sol_vec.size() == dofh.NumDofs(),
                "Size mismatch for coefficient vector");
  // Summation variable for returning result
  double s = 0.0;
  // Counter for edges on selected contact
  unsigned int ed_cnt = 0;
  const lf::mesh::Mesh& mesh{*dofh.Mesh()};
  // We cannot loop over edges, because information from dofs not located on the
  // boundary is also required. Therefore we have to loop over all cells and
  // check whether they abut the relevant boundary part.
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    // Implemented for triangles only
    // Make sure the cell is of triangular shape
    const lf::base::RefEl ref_el_type{cell->RefEl()};
    LF_ASSERT_MSG(ref_el_type == lf::base::RefEl::kTria(),
                  "contactFlux: implemented for triangles only");
    // Obtain array of edge pointers (sub-entities of co-dimension 1)
    std::span<const lf::mesh::Entity* const> sub_ent_range{
        cell->SubEntities(1)};
    // Must be three edges
    LF_ASSERT_MSG(sub_ent_range.size() == 3, "Triangle must have three edges!");
    // Check whether a relevant contact edge belongs to the cell
    std::array<bool, 3> on_contact{false};
    unsigned int cnt = 0;
    // loop over the edges and check whether they belong to the boundary
    for (lf::base::sub_idx_t j = 0; j < ref_el_type.NumSubEntities(1); ++j) {
      const lf::mesh::Entity& edge{*sub_ent_range[j]};
      if (edgeids(edge) == contact_id) {
        on_contact[j] = true;
        cnt++;
        ed_cnt++;
      }
    }
    if (cnt > 0) {
      // A contact edge belongs to the current cell
      // Compute the gradients, edge-weighted exterior normals and area
      const lf::geometry::Geometry& geo{*(cell->Geometry())};
      auto [grad_bary_coords, normals, area] =
          getTriangleGradLambdaNormals(lf::geometry::Corners(geo));

      // Compute (constant) local gradient of the finite element solution
      // DofHandler must provide three degrees of freedom per cell  piecewise
      // linear Lagrangian finite elements on triangles
      LF_ASSERT_MSG(dofh.NumLocalDofs(*cell) == 3,
                    "contactFlux: 3 dofs per triangle mandatory!");
      const auto glob_dof_idx{dofh.GlobalDofIndices(*cell)};
      const Eigen::Vector2d local_gradient{
          grad_bary_coords * (Eigen::Vector3d() << sol_vec[glob_dof_idx[0]],
                              sol_vec[glob_dof_idx[1]],
                              sol_vec[glob_dof_idx[2]])
                                 .finished()};
      // Midpoints of edges in reference coordinates
      const Eigen::MatrixXd mp_ref{
          (Eigen::MatrixXd(2, 3) << 0.5, 0.5, 0.0, 0.0, 0.5, 0.5).finished()};
      // Physical coordinates of midpoints of edges
      const Eigen::MatrixXd mp_phys{geo.Global(mp_ref)};
      for (lf::base::sub_idx_t j = 0; j < ref_el_type.NumSubEntities(1); ++j) {
        if (on_contact[j]) {
          s += normals.col(j).dot(sigma(mp_phys.col(j)) * local_gradient);
        }
      }
    }
  }  // end loop over cells
  std::cout << "Summed flux for " << ed_cnt << " edges." << std::endl;
  return s;
}  // end contact flux

/** @see \ref contactFlux
 *
 * Alternative implementation amking use of @ref MeshFunction
 */
template <typename SIGMAFUNCTION>
double contactFluxMF(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& sol_vec, SIGMAFUNCTION&& sigma,
    const lf::mesh::utils::CodimMeshDataSet<int>& edgeids, int contact_id = 0) {
  // The underlying finite element mesh
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Compute exterior edge-weighted normals
  lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> normals{
      exteriorEdgeWeightedNormals(fe_space->Mesh())};
  // Build a MeshFunction representing the gradient of the finite element
  // solution
  const lf::fe::MeshFunctionGradFE mf_grad(fe_space, sol_vec);
  // Reference coordinates of edge midpoints of a triangle
  const Eigen::MatrixXd mp_refc{
      (Eigen::Matrix<double, 2, 3>() << 0.5, 0.5, 0.0, 0.0, 0.5, 0.5)
          .finished()};
  // Variable for summing boundary flux
  double s = 0.0;
  // Counter for edges on selected contact
  unsigned int ed_cnt = 0;
  // Loop over all cells
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    const lf::base::RefEl ref_el_type{cell->RefEl()};
    LF_ASSERT_MSG(ref_el_type == lf::base::RefEl::kTria(),
                  "contactFlux: implemented for triangles only");
    // Obtain array of edge pointers (sub-entities of co-dimension 1)
    std::span<const lf::mesh::Entity* const> sub_ent_range{
        cell->SubEntities(1)};
    // Must be three edges
    LF_ASSERT_MSG(sub_ent_range.size() == 3, "Triangle must have three edges!");
    // Obtain gradient values at midpoints of edges
    const std::vector<Eigen::VectorXd> grad_at_mp{mf_grad(*cell, mp_refc)};
    // Visit edges, check flags, and add contribution to flux integral
    for (lf::base::sub_idx_t j = 0; j < ref_el_type.NumSubEntities(1); ++j) {
      const lf::mesh::Entity& edge{*sub_ent_range[j]};
      LF_ASSERT_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
                    "Not an edge!");
      if (edgeids(edge) == contact_id) {
        // Exterior normal vector
        LF_ASSERT_MSG(normals.DefinedOn(edge),
                      "Normal vector not available for " << edge);
        const Eigen::Vector2d n{normals(edge)};
        LF_ASSERT_MSG(n.norm() > 0, " zero normal vector ?");
        const auto sigma_val{sigma((cell->Geometry())->Global(mp_refc.col(j)))};
        s += n.dot(sigma_val * grad_at_mp[j]);
        ed_cnt++;
      }
    }  // end loop over edges
  }    // end loop over cells
  std::cout << "Summed flux for " << ed_cnt << " edges." << std::endl;
  return s;
}  // end contactFluxMF

/** @brief Volume based formula for the evaluation of contact fluxes
 *
 * @param fe_space pointer to an object describing a lowest-order Lagrangian
 * finite element space
 * @param sol_vec basis expansion coefficient vector of FE solution
 * @param sigma diffusion coefficient: involved in the definition of the flux
 * @param gradpsi gradient of weighting function for volume formula.
 */
template <typename SIGMAFUNCTION, typename PSIGRAD>
double stabFlux(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& sol_vec, SIGMAFUNCTION&& sigma, PSIGRAD&& gradpsi) {
  // Underlying FE mesh
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Local-to-Global map for local/global shape function indices
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Obtain quadrature rule
  const lf::quad::QuadRule quadrule{
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2)};
  // Summation variable
  double s = 0.0;
  // Loop over all cells
  for (const lf::mesh::Entity* cell : mesh.Entities(0)) {
    // Check matching of reference element (unit triangle)
    LF_VERIFY_MSG(cell->RefEl() == quadrule.RefEl(),
                  "Mismatch of reference element for " << *cell);
    // Obtain geometry information for entity
    const lf::geometry::Geometry& geo{*cell->Geometry()};
    // Compute the gradients, edge-weighted exterior normals and area
    auto [grad_bary_coords, normals, area] =
        getTriangleGradLambdaNormals(lf::geometry::Corners(geo));
    // Compute (constant) local gradient of the finite element solution
    // DofHandler must provide three degrees of freedom per cell  piecewise
    // linear Lagrangian finite elements on triangles
    LF_ASSERT_MSG(dofh.NumLocalDofs(*cell) == 3,
                  "contactFlux: 3 dofs per triangle mandatory!");
    const auto glob_dof_idx{dofh.GlobalDofIndices(*cell)};
    const Eigen::Vector2d local_gradient{
        grad_bary_coords * (Eigen::Vector3d() << sol_vec[glob_dof_idx[0]],
                            sol_vec[glob_dof_idx[1]], sol_vec[glob_dof_idx[2]])
                               .finished()};

    // Number of quadrature points
    const int P = quadrule.NumPoints();
    // Quadrature points
    const Eigen::MatrixXd zeta_ref{quadrule.Points()};
    // Map quadrature points to physical/world coordinates
    const Eigen::MatrixXd zeta{geo.Global(zeta_ref)};
    // Quadrature weights
    const Eigen::VectorXd w_ref{quadrule.Weights()};
    // Gramian determinants
    const Eigen::VectorXd gram_dets{geo.IntegrationElement(zeta_ref)};
    // Iterate over the quadrature points
    for (int l = 0; l < P; ++l) {
      const auto quadnode{zeta.col(l)};
      s += w_ref[l] * gram_dets[l] *
           (gradpsi(quadnode)).dot(sigma(quadnode) * local_gradient);
    }
  }
  return s;
}

/** @see @ref stabFlux
 *
 * Alternative implementation based on @ref MeshFunction
 */
template <typename SIGMAFUNCTION, typename PSIGRAD>
double stabFluxMF(
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space,
    const Eigen::VectorXd& sol_vec, SIGMAFUNCTION&& sigma, PSIGRAD&& gradpsi) {
  std::shared_ptr<const lf::mesh::Mesh> mesh_p{fe_space->Mesh()};
  // Coefficient function and weight function
  const lf::mesh::utils::MeshFunctionGlobal mf_sigma(sigma);
  const lf::mesh::utils::MeshFunctionGlobal mf_gradpsi(gradpsi);
  // Build a MeshFunction representing the gradient of the finite element
  // solution
  const lf::fe::MeshFunctionGradFE mf_grad(fe_space, sol_vec);
  // Mesh function representing the integrand
  const auto mf_itg{lf::mesh::utils::transpose(mf_sigma * mf_grad) *
                    mf_gradpsi};
  const double s = lf::fe::IntegrateMeshFunction(
      *mesh_p, mf_itg, [](const lf::mesh::Entity& e) {
        return lf::quad::make_QuadRule(e.RefEl(), 2);
      })(0, 0);
  return s;
}  // end stabFlux

/** @brief Main driver function
 *
 * @param basename Filename of .msh-file without suffix
 */
std::pair<double, double> computePotential(std::string basename);

}  // namespace dmxbc

#endif
