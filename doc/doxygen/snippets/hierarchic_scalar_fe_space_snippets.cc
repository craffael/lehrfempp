/**
 * @file
 * @brief Snippets that illustrate how HierarchicScalarFESpace can be used.
 * @author Raffael Casagrande
 * @date   2021-01-25 01:37:04
 * @copyright MIT License
 */

#include <lf/fe/fe.h>
#include <lf/io/io.h>

#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::fe {

void L2Projection() {
  //! [Laplace]
  // Get a mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Create HierarchicalFESpace with uniform degree 5.
  const unsigned degree = 5;
  const auto fe_space =
      std::make_shared<lf::fe::HierarchicScalarFESpace<double>>(mesh_p, degree);

  // define diffusion coefficient
  const lf::mesh::utils::MeshFunctionConstant mf_alpha(1);
  // define rhs load
  auto mf_load = lf::mesh::utils::MeshFunctionConstant(1.);

  // Assemble the system matrix and right hand side
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(fe_space->LocGlobMap().NumDofs());
  const assemble::DofHandler &dofh = fe_space->LocGlobMap();
  assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());

  DiffusionElementMatrixProvider element_matrix_provider(fe_space, mf_alpha);
  AssembleMatrixLocally(0, dofh, dofh, element_matrix_provider, A_COO);
  ScalarLoadElementVectorProvider element_vector_provider(fe_space, mf_load);
  AssembleVectorLocally(0, dofh, element_vector_provider, rhs);

  // Enforce zero dirichlet boundary conditions

  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p);
  const auto selector = [&](unsigned int idx) -> std::pair<bool, double> {
    const auto &entity = dofh.Entity(idx);
    return {entity.Codim() > 0 && boundary(entity), 0};
  };
  FixFlaggedSolutionComponents(selector, A_COO, rhs);

  // Solve the LSE using the cholesky decomposition
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
  const Eigen::VectorXd solution = solver.solve(rhs);

  // visualize the solution:
  io::VtkWriter vtk(mesh_p, "solution.vtk", 0, 5);
  vtk.WritePointData("solution", MeshFunctionFE(fe_space, solution));
  //! [Laplace]
}

}  // namespace lf::fe
