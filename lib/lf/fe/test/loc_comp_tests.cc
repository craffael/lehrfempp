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
#include <lf/fe/test_utils/test_utils.h>
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
      std::make_shared<lf::fe::HierarchicScalarFESpace<double>>(mesh_p, degree);

  // The analytic solution
  const auto u = [](const Eigen::VectorXd &x) -> double {
    return std::sin(base::kPi * x[0]) * std::sin(base::kPi * x[1]);
  };
  const lf::mesh::utils::MeshFunctionGlobal mf_u(u);

  // Define the load function of the manufactured solution
  const auto load = [](const Eigen::Vector2d &x) -> double {
    return 2 * base::kPi * base::kPi * std::sin(base::kPi * x[0]) *
           std::sin(base::kPi * x[1]);
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
      std::make_shared<lf::fe::HierarchicScalarFESpace<double>>(mesh_p, degree);

  // Define some right hand side
  const auto load = [](const Eigen::Vector2d &x) -> double {
    return 2 * base::kPi * base::kPi * std::sin(base::kPi * x[0]) *
           std::sin(base::kPi * x[1]);
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

template <class SCALAR_TEST, class SCALAR_TRIAL, class EMP, class EVALUATOR>
void CheckElementMatrixProvider(const ScalarFESpace<SCALAR_TEST> &fes_test,
                                const ScalarFESpace<SCALAR_TRIAL> &fes_trial,
                                const EMP &emp, EVALUATOR &&evaluator,
                                base::dim_t codim = 0) {
  Eigen::Matrix<SCALAR_TEST, Eigen::Dynamic, 1> coeff_test(
      fes_test.LocGlobMap().NumDofs());
  Eigen::Matrix<SCALAR_TRIAL, Eigen::Dynamic, 1> coeff_trial(
      fes_trial.LocGlobMap().NumDofs());
  coeff_test.setZero();
  coeff_trial.setZero();

  using scalar_matrix_t = decltype(evaluator(coeff_test, coeff_trial));

  // assemble the matrix:
  assemble::COOMatrix<scalar_matrix_t> coo_matrix(coeff_test.rows(),
                                                  coeff_trial.rows());
  assemble::AssembleMatrixLocally(codim, fes_trial.LocGlobMap(),
                                  fes_test.LocGlobMap(), emp, coo_matrix);
  Eigen::Matrix<scalar_matrix_t, Eigen::Dynamic, Eigen::Dynamic> system_matrix =
      coo_matrix.makeDense();

  for (base::size_type i = 0; i < coeff_test.rows(); ++i) {
    coeff_test(i) = 1.;
    for (base::size_type j = 0; j < coeff_trial.rows(); ++j) {
      coeff_trial(j) = 1.;
      scalar_matrix_t result = evaluator(coeff_test, coeff_trial);
      auto entry = system_matrix(i, j);
      ASSERT_LT(std::abs(result - system_matrix(i, j)), 1e-7);
      coeff_trial(j) = 0.;
    }
    coeff_test(i) = 0;
  }
}

template <class SCALAR, class EVP, class EVALUATOR>
void CheckEntityVectorProvider(const ScalarFESpace<SCALAR> &fes_test,
                               const EVP &evp, EVALUATOR &&evaluator,
                               base::dim_t codim) {
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> coeff_test(
      fes_test.LocGlobMap().NumDofs());
  coeff_test.setZero();

  using scalar_vector_t = decltype(evaluator(coeff_test));

  // assemble the vector:
  Eigen::Matrix<scalar_vector_t, Eigen::Dynamic, 1> vector(coeff_test.rows());
  vector.setZero();
  assemble::AssembleVectorLocally(codim, fes_test.LocGlobMap(), evp, vector);

  for (base::size_type i = 0; i < coeff_test.rows(); ++i) {
    coeff_test(i) = 1;
    auto result = evaluator(coeff_test);
    auto entry = vector(i);
    ASSERT_LT(std::abs(result - entry), 1e-7);
    coeff_test(i) = 0;
  }
}

TEST(lf_fe, DiffusionElementMatrixProviderComplexCoeff) {
  // get a mesh of affine elements so that the element matrices can be
  // calculated exactly.
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(5);

  //// Generate a variable order hierarchic fe space:
  auto fe_space = std::make_shared<HierarchicScalarFESpace<double>>(
      mesh, [&](const mesh::Entity &e) { return (mesh->Index(e)) % 2 + 1; });

  auto max_degree = 2;

  //// complex diffusion coefficient:
  auto alpha = [&](const mesh::Entity &e, const Eigen::MatrixXd &local) {
    return std::vector<double>(local.cols(), mesh->Index(e) + 10.);
  };

  DiffusionElementMatrixProvider emp(fe_space, alpha);

  // calculate every element of the matrix with MeshFunctions:
  auto evaluator = [&](const Eigen::VectorXd &test_coeff,
                       const Eigen::VectorXd &trial_coeff) {
    auto mf_grad_test = MeshFunctionGradFE(fe_space, test_coeff);
    auto mf_grad_trial = MeshFunctionGradFE(fe_space, trial_coeff);
    return lf::fe::IntegrateMeshFunction(
        *mesh, transpose(mf_grad_test) * (alpha * mf_grad_trial),
        2 * (max_degree))(0);
  };

  CheckElementMatrixProvider(*fe_space, *fe_space, emp, evaluator);
}

TEST(lf_fe, DiffusionElementMatrixProviderTensorCoeff) {
  // get a mesh of affine elements so that the element matrices can be
  // calculated exactly.
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(5);
  //// Generate a variable order hierarchic fe space:
  auto fe_space = std::make_shared<HierarchicScalarFESpace<double>>(
      mesh, [&](const mesh::Entity &e) { return (mesh->Index(e)) % 2 + 1; });

  auto max_degree = 2;

  auto alpha2 = mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d &x) {
    return (Eigen::Matrix2d() << 1, x.x(), x.y(), 2).finished();
  });
  DiffusionElementMatrixProvider emp2(fe_space, alpha2);
  auto evaluator2 = [&](const Eigen::VectorXd &test_coeff,
                        const Eigen::VectorXd &trial_coeff) {
    auto mf_grad_test = MeshFunctionGradFE(fe_space, test_coeff);
    auto mf_grad_trial = MeshFunctionGradFE(fe_space, trial_coeff);
    return lf::fe::IntegrateMeshFunction(
        *mesh, transpose(mf_grad_test) * (alpha2 * mf_grad_trial),
        2 * (max_degree))(0);
  };
  CheckElementMatrixProvider(*fe_space, *fe_space, emp2, evaluator2);
}

TEST(lf_fe, DiffusionElementMatrixProviderComplexFESpace) {
  // get a mesh of affine elements so that the element matrices can be
  // calculated exactly.
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(5);

  //// complex diffusion coefficient:
  auto alpha = [&](const mesh::Entity &e, const Eigen::MatrixXd &local) {
    return std::vector<double>(local.cols(), mesh->Index(e) + 10.);
  };

  // with a complex fe_space:
  auto fes_complex = test_utils::MakeComplexLagrangeO1FeSpace(mesh);
  DiffusionElementMatrixProvider emp3(fes_complex, alpha);
  auto evaluator3 = [&](const Eigen::VectorXcd &test_coeff,
                        const Eigen::VectorXcd &trial_coeff) {
    auto mf_grad_test = MeshFunctionGradFE(fes_complex, test_coeff);
    auto mf_grad_trial = MeshFunctionGradFE(fes_complex, trial_coeff);
    return lf::fe::IntegrateMeshFunction(
        *mesh, adjoint(mf_grad_test) * (alpha * mf_grad_trial), 2)(0);
  };
  CheckElementMatrixProvider(*fes_complex, *fes_complex, emp3, evaluator3);
}

TEST(lf_fe, MassElementMatrixProvider) {
  auto mesh = mesh::test_utils::GenerateHybrid2DTestMesh(5);

  auto g = [&](const mesh::Entity &e, const Eigen::MatrixXd &local) {
    return std::vector<double>(local.cols(), mesh->Index(e) + 10.);
  };

  auto fes = std::make_shared<HierarchicScalarFESpace<double>>(mesh, 1);
  MassElementMatrixProvider emp(fes, g);
  auto evaluator = [&](const Eigen::VectorXd &test_coeff,
                       Eigen::VectorXd &trial_coeff) {
    auto mf_test = MeshFunctionFE(fes, test_coeff);
    auto mf_trial = MeshFunctionFE(fes, trial_coeff);
    return IntegrateMeshFunction(*mesh, mf_test * mf_trial * g, 2);
  };
  CheckElementMatrixProvider(*fes, *fes, emp, evaluator);
}

TEST(lf_fe, MassElementMatrixProviderComplex) {
  auto mesh = mesh::test_utils::GenerateHybrid2DTestMesh(5);

  auto g = mesh::utils::MeshFunctionConstant(std::complex<double>{1, 2});

  auto fes = fe::test_utils::MakeComplexLagrangeO1FeSpace(mesh);
  MassElementMatrixProvider emp(fes, g);
  auto evaluator = [&](const Eigen::VectorXcd &test_coeff,
                       Eigen::VectorXcd &trial_coeff) {
    auto mf_test = MeshFunctionFE(fes, test_coeff);
    auto mf_trial = MeshFunctionFE(fes, trial_coeff);
    return IntegrateMeshFunction(*mesh, conjugate(mf_test) * mf_trial * g, 2);
  };
  CheckElementMatrixProvider(*fes, *fes, emp, evaluator);
}

TEST(lf_fe, MassEdgeMatrixProvider) {
  auto mesh = mesh::test_utils::GenerateHybrid2DTestMesh(5);

  auto g = [&](const mesh::Entity &e, const Eigen::MatrixXd &local) {
    return std::vector<double>(local.cols(), mesh->Index(e) + 10.);
  };
  auto fes = std::make_shared<HierarchicScalarFESpace<double>>(mesh, 1);
  MassEdgeMatrixProvider emp(fes, g);
  auto evaluator = [&](const Eigen::VectorXd &test_coeff,
                       const Eigen::VectorXd &trial_coeff) {
    auto mf_test = MeshFunctionFE(fes, test_coeff);
    auto mf_trial = MeshFunctionFE(fes, trial_coeff);
    return IntegrateMeshFunction(*mesh, mf_test * g * mf_trial, 2,
                                 base::PredicateTrue{}, 1);
  };
  CheckElementMatrixProvider(*fes, *fes, emp, evaluator, 1);
}

TEST(lf_fe, MassEdgeMatrixProviderComplex) {
  auto mesh = mesh::test_utils::GenerateHybrid2DTestMesh(5);

  auto g = mesh::utils::MeshFunctionConstant(std::complex<double>{1, 2});
  auto fes = test_utils::MakeComplexLagrangeO1FeSpace(mesh);
  auto selector = [&](const mesh::Entity &e) {
    LF_ASSERT_MSG(e.Codim() == 1, "This is not an edge.");
    return mesh->Index(e) < 7;
  };
  MassEdgeMatrixProvider emp(fes, g, selector);
  auto evaluator = [&](const Eigen::VectorXcd &test_coeff,
                       const Eigen::VectorXcd &trial_coeff) {
    auto mf_test = MeshFunctionFE(fes, test_coeff);
    auto mf_trial = MeshFunctionFE(fes, trial_coeff);
    return IntegrateMeshFunction(*mesh, conjugate(mf_test) * g * mf_trial, 2,
                                 selector, 1);
  };
  CheckElementMatrixProvider(*fes, *fes, emp, evaluator, 1);
}

TEST(lf_fe, ScalarLoadElementVectorProvider_complex) {
  auto mesh = mesh::test_utils::GenerateHybrid2DTestMesh(5);

  auto g = mesh::utils::MeshFunctionConstant(std::complex<double>{1, 2});
  auto fes = test_utils::MakeComplexLagrangeO1FeSpace(mesh);
  ScalarLoadElementVectorProvider evp(fes, g);
  auto evaluator = [&](const Eigen::VectorXcd &test_coeff) {
    auto mf_test = MeshFunctionFE(fes, test_coeff);
    return IntegrateMeshFunction(*mesh, conjugate(mf_test) * g, 2);
  };
  CheckEntityVectorProvider(*fes, evp, evaluator, 0);
}

TEST(lf_fe, ScalarLoadEdgeVectorProvider_complex) {
  auto mesh = mesh::test_utils::GenerateHybrid2DTestMesh(5);

  auto g = mesh::utils::MeshFunctionConstant(std::complex<double>{1, 2});
  auto fes = test_utils::MakeComplexLagrangeO1FeSpace(mesh);
  auto selector = [&](const mesh::Entity &e) {
    LF_ASSERT_MSG(e.Codim() == 1, "This is not an edge.");
    return mesh->Index(e) < 7;
  };
  ScalarLoadEdgeVectorProvider evp(fes, g, selector);
  auto evaluator = [&](const Eigen::VectorXcd &test_coeff) {
    auto mf_test = MeshFunctionFE(fes, test_coeff);
    return IntegrateMeshFunction(*mesh, conjugate(mf_test) * g, 2, selector, 1);
  };
  CheckEntityVectorProvider(*fes, evp, evaluator, 1);
}

}  // namespace lf::fe::test
