/**
 * @file
 * @brief Unit tests for local element matrix builders
 * @author Tobias Rohner
 * @date January 2021
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <iostream>

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>

namespace lf::fe::test {

TEST(lf_fe, diffusion_mat_test) {
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Use a relatively high polynomial degree
  const unsigned degree = 15;
  const auto fe_space =
      std::make_shared<lf::fe::FeSpaceHierarchic<double>>(mesh_p, degree);

  // The analytic solution
  const auto u = [](const Eigen::VectorXd &x) -> double {
    return std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
  };
  const lf::mesh::utils::MeshFunctionGlobal mf_u(u);

  // Define the load function of the manufactured solution
  const auto load = [](const Eigen::Vector2d &x) -> double {
    return 2 * M_PI * M_PI * std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
  };
  const lf::mesh::utils::MeshFunctionGlobal mf_load(load);

  // Assemble the system matrix and right hand side
  const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
  lf::assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());
  std::cout << "> Assembling System Matrix" << std::endl;
  const lf::mesh::utils::MeshFunctionConstant<double> mf_alpha(1);
  lf::fe::DiffusionElementMatrixProvider element_matrix_provider(fe_space,
                                                                 mf_alpha);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, element_matrix_provider,
                                      A_COO);
  std::cout << "> Assembling right Hand Side" << std::endl;
  lf::fe::ScalarLoadElementVectorProvider element_vector_provider(fe_space,
                                                                  mf_load);
  lf::assemble::AssembleVectorLocally(0, dofh, element_vector_provider, rhs);

  // Enforce zero dirichlet boundary conditions
  std::cout << "> Enforcing Boundary Conditions" << std::endl;
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p);
  const auto selector = [&](unsigned int idx) -> std::pair<bool, double> {
    const auto &entity = dofh.Entity(idx);
    return {entity.Codim() > 0 && boundary(entity), 0};
  };
  lf::assemble::FixFlaggedSolutionComponents(selector, A_COO, rhs);

  // Solve the LSE using the cholesky decomposition
  std::cout << "> Solving LSE" << std::endl;
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
  const Eigen::VectorXd solution = solver.solve(rhs);
  const lf::fe::MeshFunctionFE<double, double> mf_numeric(fe_space, solution);
  const lf::fe::MeshFunctionGradFE<double, double> mf_numeric_grad(fe_space,
                                                                   solution);

  // Compute the L2 error
  std::cout << "> Computing Error Norms" << std::endl;
  const auto qr_segment =
      lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), 2 * degree - 1);
  const auto qr_tria =
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2 * degree - 1);
  const auto qr_quad =
      lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 2 * degree - 1);
  const auto quadrule_provider = [&](const lf::mesh::Entity &entity) {
    const lf::base::RefEl refel = entity.RefEl();
    switch (refel) {
      case lf::base::RefEl::kTria():
        return qr_tria;
      case lf::base::RefEl::kSegment():
        return qr_segment;
      case lf::base::RefEl::kQuad():
        return qr_quad;
      default:
        return lf::quad::make_QuadRule(refel, 2 * degree - 1);
    }
  };
  const double L2_err = std::sqrt(lf::fe::IntegrateMeshFunction(
      *mesh_p, lf::mesh::utils::squaredNorm(mf_u - mf_numeric),
      quadrule_provider));

  // Assert that the L2 error is small enough
  ASSERT_NEAR(L2_err, 0, 1e-5);
}

TEST(lf_fe, mass_mat_test) {
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Use a relatively high polynomial degree
  const unsigned degree = 15;
  const auto fe_space =
      std::make_shared<lf::fe::FeSpaceHierarchic<double>>(mesh_p, degree);

  // Define some right hand side
  const auto load = [](const Eigen::Vector2d &x) -> double {
    return 2 * M_PI * M_PI * std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
  };
  const lf::mesh::utils::MeshFunctionGlobal mf_load(load);

  // Assemble the system matrix and right hand side
  const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
  lf::assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());
  std::cout << "> Assembling System Matrix" << std::endl;
  const lf::mesh::utils::MeshFunctionConstant<double> mf_gamma(1);
  lf::fe::MassElementMatrixProvider element_matrix_provider(fe_space, mf_gamma);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, element_matrix_provider,
                                      A_COO);
  std::cout << "> Assembling right Hand Side" << std::endl;
  lf::fe::ScalarLoadElementVectorProvider element_vector_provider(fe_space,
                                                                  mf_load);
  lf::assemble::AssembleVectorLocally(0, dofh, element_vector_provider, rhs);

  // Enforce zero dirichlet boundary conditions
  std::cout << "> Enforcing Boundary Conditions" << std::endl;
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p);
  const auto selector = [&](unsigned int idx) -> std::pair<bool, double> {
    const auto &entity = dofh.Entity(idx);
    return {entity.Codim() > 0 && boundary(entity), 0};
  };
  lf::assemble::FixFlaggedSolutionComponents(selector, A_COO, rhs);

  // Solve the LSE using the cholesky decomposition
  std::cout << "> Solving LSE" << std::endl;
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
  const Eigen::VectorXd solution = solver.solve(rhs);
  const lf::fe::MeshFunctionFE<double, double> mf_numeric(fe_space, solution);
  const lf::fe::MeshFunctionGradFE<double, double> mf_numeric_grad(fe_space,
                                                                   solution);

  // Compute the L2 error
  std::cout << "> Computing Error Norms" << std::endl;
  const auto qr_segment =
      lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), 2 * degree - 1);
  const auto qr_tria =
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2 * degree - 1);
  const auto qr_quad =
      lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 2 * degree - 1);
  const auto quadrule_provider = [&](const lf::mesh::Entity &entity) {
    const lf::base::RefEl refel = entity.RefEl();
    switch (refel) {
      case lf::base::RefEl::kTria():
        return qr_tria;
      case lf::base::RefEl::kSegment():
        return qr_segment;
      case lf::base::RefEl::kQuad():
        return qr_quad;
      default:
        return lf::quad::make_QuadRule(refel, 2 * degree - 1);
    }
  };
  const double L2_err = std::sqrt(lf::fe::IntegrateMeshFunction(
      *mesh_p, lf::mesh::utils::squaredNorm(mf_load - mf_numeric),
      quadrule_provider));

  // Assert that the L2 error is small enough
  ASSERT_NEAR(L2_err, 0, 1e-5);
}

}  // namespace lf::fe::test
