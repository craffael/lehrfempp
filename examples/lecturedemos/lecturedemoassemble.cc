/**
 * @file
 * @brief demonstration of assembly of Galerkin linear system in LehrFEM++
 * assemble module; meant to provide sample codes for lecture document
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include "lecturedemoassemble.h"

#include <Eigen/Eigen>

namespace lecturedemo {

// clang-format off
/* SAM_LISTING_BEGIN_1 */
  Eigen::Matrix<double, 3, 3> LinFELaplaceElemMatProvider::Eval( // NOLINT
    const lf::mesh::Entity &tria) {
    // Throw error in case no triangular cell
    LF_VERIFY_MSG(tria.RefEl() == lf::base::RefEl::kTria(),
		  "Unsupported cell type " << tria.RefEl());
    // Obtain vertex coordinates of the triangle in a 2x3 matrix
    const auto vertices = lf::geometry::Corners(*(tria.Geometry()));
    LF_ASSERT_MSG((vertices.cols() == 3) && (vertices.rows() == 2),
		  "Invalid vertex coordinate " << vertices.rows() << "x"
		  << vertices.cols() << " matrix");

    // Set up an auxiliary 3x3-matrix with a leading column 1 and
    // the vertex coordinates in its right 3x2 block
    Eigen::Matrix<double, 3, 3> X;  // temporary matrix
    X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
    X.block<3, 2>(0, 1) = vertices.transpose();
    // The determinant of the auxiliary matrix also supplies the determinant
    const double area = 0.5 * std::abs(X.determinant());
    // Compute the gradients of the barycentric coordinate functions
    // and store them in the columns of a 2x3 matrix grad\_bary\_coords
    Eigen::Matrix<double, 2, 3>
      grad_bary_coords{X.inverse().block<2, 3>(1, 0)};

    // Since the gradients are constant, local integration is easy
    return ((area * grad_bary_coords.transpose()) * grad_bary_coords);
  }

/* SAM_LISTING_END_1 */
// clang-format on

// clang-format off
/* SAM_LISTING_BEGIN_2 */
  Eigen::Matrix2d LinFEMassEdgeMatProvider::Eval(const lf::mesh::Entity &edge) { // NOLINT
    LF_VERIFY_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
		  "Unsupported edge type " << edge.RefEl());
    // Obtain endpoint coordinates of the triangle in a 2x3 matrix
    const auto endpoints = lf::geometry::Corners(*(edge.Geometry()));
    // Compute length of edge
    const double edge_length = (endpoints.col(1) - endpoints.col(0)).norm();
    // Diagonal and off-diagonal entries of edge mass matrix
    const double m1 = edge_length * 1.0 / 3.0;
    const double m2 = edge_length * 1.0 / 6.0;
    return ((Eigen::Matrix2d(2, 2) << m1, m2, m2, m1).finished());
  }

/* SAM_LISTING_END_2 */
// clang-format on

// Main driver function
void lecturedemoassemble() {
  std::cout << "LehrFEM++ demo: assembly of Galerkin linear system" << '\n';
  // Obtain a purely triangular mesh from the collection of LehrFEM++'s
  // built-in meshes
  const std::shared_ptr<lf::mesh::Mesh> mesh_p{
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3)};
  // Define the right hand side function as a simple lambda function of
  // point world coordinates (passed as a column vector)
  auto f = [](const Eigen::Vector2d &x) { return x[0] - x[1]; };

  // Initialization of local-to-global index mapping for linear finite elements
  const lf::assemble::UniformFEDofHandler dof_handler(
      mesh_p, {{lf::base::RefEl::kPoint(), 1}});
  // clang-format off
  /* SAM_LISTING_BEGIN_3 */
  // Query dimension of the finite element space, equal to the number of nodes
  const size_type N_dofs(dof_handler.NumDofs());
  // Matrix in \samemp{triplet format} holding temporary Galerkin matrix 
  lf::assemble::COOMatrix<double> mat(N_dofs, N_dofs);

  // Initialize objects for local computations
  LinFELaplaceElemMatProvider loc_mat_laplace{};
  LinFEElemVecProvider<decltype(f)> loc_vec_sample(f);

  // Building the Galerkin matrix (trial space = test space)
  // for the pure Neumann Laplacian, assembly over cells
  mat = lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
        0, dof_handler, loc_mat_laplace);
  // Filling the right-hand-side vector, assembly over cells
  auto rhsvec = lf::assemble::AssembleVectorLocally<Eigen::VectorXd>(
        0, dof_handler, loc_vec_sample);
  /* SAM_LISTING_END_3 */
  // clang-format on

  // Add boundary mass matrix
  /* SAM_LISTING_BEGIN_5 */
  // Flag all edge (co-dimension-1 entities) on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  LinFEMassEdgeMatProvider loc_edge_mass(bd_flags);
  // Add boundary contributions to Galerkin matrix stored in 'mat'
  // Assembly covers edges (co-dimensions-1 entities)!
  lf::assemble::AssembleMatrixLocally(1, dof_handler, dof_handler,
                                      loc_edge_mass, mat);
  /* SAM_LISTING_END_5 */
  /* SAM_LISTING_BEGIN_4 */
  // Convert the matrix from triplet format to CRS format
  const Eigen::SparseMatrix<double> A(mat.makeSparse());  // NOLINT
  /* SAM_LISTING_END_4 */

  // Print Galerkin matrix, posible because it is small
  std::cout << "Galerkin matrix: " << '\n' << Eigen::MatrixXd(A) << '\n';
  // Print right-hand-side vector (as row vector)
  std::cout << "R.h.s. vector: " << '\n'
            << '[' << rhsvec.transpose() << ']' << '\n';

  // Solve linear system using a direct sparse elimination solver
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  Eigen::VectorXd sol_vec = solver.solve(rhsvec);
  if (solver.info() != Eigen::Success) {
    std::cout << "solver failed!" << '\n';
  } else {
    // Compute energy norm of solution
    const double energy = sol_vec.transpose() * A * sol_vec;
    std::cout << "Energy norm of solution = " << energy << '\n';
  }
}

// Main driver function
void lecturedemoDirichlet() {
  std::cout << "LehrFEM++ demo: treatment of Dirichlet boundary conditions"
            << '\n';
  // This is the exact solution, a harmonic function, from which we also derive
  // the Dirichlet data
  auto u_sol = [](Eigen::Vector2d x) -> double { return (x[0] - x[1]); };

  // Obtain a purely triangular mesh from the collection of LehrFEM++'s
  // built-in meshes
  const std::shared_ptr<lf::mesh::Mesh> mesh_p{
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3)};

  // In this demo routine we just solve a boundary value problem for the
  // Laplace equation with non-zero Dirichlet boundary conditions, that is the
  // right-hand-side source function vanishes and the exact solution will
  // be a harmonic function.

  // clang-format off
  /* SAM_LISTING_BEGIN_7 */
  // Initialization of local-to-global index mapping for linear finite elements
  const lf::assemble::UniformFEDofHandler dof_handler(
      mesh_p, {{lf::base::RefEl::kPoint(), 1}});
  // Query dimension of the finite element space, equal to the number of nodes
  const size_type N_dofs(dof_handler.NumDofs());
  // Matrix in \samemp{triplet format} holding temporary Galerkin matrix
  lf::assemble::COOMatrix<double> mat(N_dofs, N_dofs);
  // Initialize objects for local computation of element matrices for $\cob{-\Delta}$
  LinFELaplaceElemMatProvider loc_mat_laplace{};

  // Building the Galerkin matrix (trial space = test space)
  // for the pure Neumann Laplacian, assembly over cells
  mat = lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dof_handler, loc_mat_laplace);
  // Zero right-hand side vector
  Eigen::VectorXd rhsvec(N_dofs);
  rhsvec.setZero();

  // Treatment of Dirichlet boundary conditions $\cob{g=\rst{u}{\partial\Omega}}$
  // Flag all nodes on the boundary (and only those)
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Set up predicate: Run through all global shape functions and check whether
  // they are associated with an entity on the boundary, store Dirichlet data.
  std::vector<std::pair<bool, double>> ess_dof_select{};
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N_dofs; ++dofnum) {
    const lf::mesh::Entity &dof_node{dof_handler.Entity(dofnum)};
    const Eigen::Vector2d node_pos{
        lf::geometry::Corners(*dof_node.Geometry()).col(0)};
    const double g_val = u_sol(node_pos);
    if (bd_flags(dof_node)) {
      // Dof associated with a entity on the boundary: "essential dof"
      // The value of the dof should be set to the value of the function
      // u at the location of the node.
      ess_dof_select.emplace_back(true, g_val);
    } else {
      // Interior node, also store value of solution for comparison purposes
      ess_dof_select.emplace_back(false, g_val);
    }
  }
  // modify linear system of equations 
  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&ess_dof_select](glb_idx_t dof_idx) -> std::pair<bool, double> {
        return ess_dof_select[dof_idx];
      },
      mat, rhsvec);

  // Convert the matrix from triplet format to CRS format
  const Eigen::SparseMatrix<double> A(mat.makeSparse());
  /* SAM_LISTING_END_7 */
  // clang-format on

  // Solve linear system using a direct sparse elimination solver
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  Eigen::VectorXd sol_vec = solver.solve(rhsvec);
  if (solver.info() != Eigen::Success) {
    std::cout << "solver failed!" << '\n';
  } else {
    // Compute the L2 norm of nodal error cell by cell by the 2D trapezoidal
    // rule
    double nodal_err = 0.0;
    for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
      const std::span<const lf::assemble::gdof_idx_t> cell_dof_idx(
          dof_handler.GlobalDofIndices(*cell));
      LF_ASSERT_MSG(dof_handler.NumLocalDofs(*cell) == cell->RefEl().NumNodes(),
                    "Inconsistent node number");
      const lf::base::size_type num_nodes = cell->RefEl().NumNodes();
      double sum = 0.0;
      for (int k = 0; k < num_nodes; ++k) {
        sum += std::pow(
            (sol_vec[cell_dof_idx[k]] - ess_dof_select[cell_dof_idx[k]].second),
            2);
      }
      nodal_err += lf::geometry::Volume(*cell->Geometry()) * (sum / num_nodes);
    }
    std::cout << "L2 error = " << std::sqrt(nodal_err) << '\n';
  }
}

}  // namespace lecturedemo
