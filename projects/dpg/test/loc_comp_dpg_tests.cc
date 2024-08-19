/**
 * @file
 * @brief LOCAL tests for the (concrete) element matrix providers.
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>

#include <iostream>

#include "../loc_comp_dpg.h"
#include "../product_fe_space.h"
#include "../product_fe_space_factory.h"

namespace projects::dpg::test {

//*************************************************
// Diffusion element matrix provider
//*************************************************

TEST(DiffusionElementMatrixProvider, mass_mat_test_1) {
  std::cout << ">>> DiffusionElementMatrixProvider: Test of computation of "
               "element matrices (alpha= 1) \n";
  // Build the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Create a trial space of quadratic  H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
  auto u = factory_trial.AddH1Component(2);
  auto fe_space_trial = factory_trial.Build();

  // Create a test space of cubic H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_test(mesh_p);
  auto v = factory_test.AddH1Component(3);
  auto fe_space_test = factory_test.Build();

  auto alpha = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> double { return 1.0; });
  DiffusionElementMatrixProvider provider(fe_space_trial, fe_space_test, u, v,
                                          alpha);

  // compute element matrices on all cells:
  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);
    // multiply matrix from left and right with one vectors
    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);
    const double val = one_vec_r.dot(M * one_vec_c);

    std::cout << *cell << ": (" << M.rows() << ", " << M.cols() << ")" << ": "
              << val << std::endl;
    // The multiplication with one vectors means, that the gradients are zero
    // (constant function), such that the value should be zero to.
    EXPECT_NEAR(val, 0, 1.0E-10) << "missmatch for cell " << *cell;
  }
}

TEST(DiffusionElementMatrixProvider, mass_mat_test_2) {
  std::cout << ">>> DiffusionElementMatrixProvider: Test of computation of "
               "element matrices (alpha=[1,1;1,1]) \n";
  // Build the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Create a trial space of quadratic  H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
  auto u = factory_trial.AddH1Component(2);
  auto fe_space_trial = factory_trial.Build();

  // Create a test space of cubic H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_test(mesh_p);
  auto v = factory_test.AddH1Component(3);
  auto fe_space_test = factory_test.Build();

  auto alpha = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> Eigen::Matrix2d {
        return (Eigen::Matrix2d(2, 2) << 1.0, 1.0, 1.0, 1.0).finished();
      });
  DiffusionElementMatrixProvider provider(fe_space_trial, fe_space_test, u, v,
                                          alpha);

  // compute element matrices on all cells
  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);
    // multiply matrix from left and right with one vectors
    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);
    const double val = one_vec_r.dot(M * one_vec_c);
    std::cout << *cell << ": (" << M.rows() << ", " << M.cols() << "): " << val
              << std::endl;
    // The multiplication with the one vectors means, that the gradients are
    // zero (constant function), so that the value should be zero
    EXPECT_NEAR(val, 0, 1.0E-10) << "missmatch for cell " << *cell;
  }
}

//*************************************************
// Reaction element matrix provider
//*************************************************
TEST(ReactionElementMatrixProvider, mass_mat_test) {
  std::cout << ">>> ReactionElementMatrixProvider: Test of computation of "
               "element matrices (gamma= 1)\n";
  // build the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Create a trial space of quadratic  H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
  auto u = factory_trial.AddH1Component(2);
  auto fe_space_trial = factory_trial.Build();

  // Create a test space of cubic H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_test(mesh_p);
  auto v = factory_test.AddH1Component(3);
  auto fe_space_test = factory_test.Build();

  auto gamma = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> double { return 1.0; });
  ReactionElementMatrixProvider provider(fe_space_trial, fe_space_test, u, v,
                                         gamma);

  // Compute element matrices on all cells
  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);
    // multiply with one vector from left and right
    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);
    const double val = one_vec_r.dot(M * one_vec_c);

    std::cout << *cell << ": (" << M.rows() << ", " << M.cols() << ")" << ": "
              << val << std::endl;
    // We the value corresponds to the integral of the constant one function
    // over the cell, so we expect it to be equal to the size of the cell.
    EXPECT_NEAR(val, lf::geometry::Volume(*(cell->Geometry())), 1.0E-10)
        << "missmatch for cell " << *cell;
  }
}

//*************************************************
// Convection element matrix provider
//*************************************************
TEST(ConvectionElementMatrixProvider, mass_mat_test_1) {
  std::cout << ">>> ConvectionElementMatrixProvider: Test of computation of "
               "element matrices\n";
  std::cout << "beta_1 = [1.0;1.0], beta_2 = [0.0,0.0], v=[1,...,1] "
               "u=Eigen::Random \n";
  // build test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Create a trial space of quadratic  H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
  auto u = factory_trial.AddH1Component(2);
  auto fe_space_trial = factory_trial.Build();

  // Create a test space of cubic H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_test(mesh_p);
  auto v = factory_test.AddH1Component(3);
  auto fe_space_test = factory_test.Build();

  auto beta_1 = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> Eigen::Vector2d {
        return (Eigen::Vector2d(2) << 1.0, 1.0).finished();
      });
  auto beta_2 = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> Eigen::Vector2d {
        return (Eigen::Vector2d(2) << 0.0, 0.0).finished();
      });
  ConvectionElementMatrixProvider provider(fe_space_trial, fe_space_test, u, v,
                                           beta_1, beta_2);

  // compute element matrices on all cells
  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);
    // multiply with vectors from left and right
    Eigen::VectorXd u = Eigen::VectorXd::Random(M.cols());
    Eigen::VectorXd v = Eigen::VectorXd::Constant(M.rows(), 1.0);
    const double val = v.dot(M * u);
    std::cout << *cell << ": (" << M.rows() << ", " << M.cols() << ")" << ": "
              << val << std::endl;
    // v is constant and thus its gradient is zero, thus val should be zero
    EXPECT_NEAR(val, 0.0, 1.0E-10) << "missmatch on cell " << *cell;
  }
}

TEST(ConvectionElementMatrixProvider, mass_mat_test_2) {
  std::cout << ">>> ConvectionElementMatrixProvider: Test of computation of "
               "element matrices\n";
  std::cout << "beta_1 = [0.0;0.0], beta_2 = [1.0,1.0], v=Eigen::Random "
               "u=[1,...,1] \n";
  // build test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Create a trial space of quadratic  H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
  auto u = factory_trial.AddH1Component(2);
  auto fe_space_trial = factory_trial.Build();

  // Create a test space of cubic H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_test(mesh_p);
  auto v = factory_test.AddH1Component(3);
  auto fe_space_test = factory_test.Build();

  auto beta_1 = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> Eigen::Vector2d {
        return (Eigen::Vector2d(2) << 0.0, 0.0).finished();
      });
  auto beta_2 = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> Eigen::Vector2d {
        return (Eigen::Vector2d(2) << 1.0, 1.0).finished();
      });
  ConvectionElementMatrixProvider provider(fe_space_trial, fe_space_test, u, v,
                                           beta_1, beta_2);

  // compute element matrices on all cells
  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);
    // multiply with vectors from left and right
    Eigen::VectorXd u = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd v = Eigen::VectorXd::Random(M.rows());
    const double val = v.dot(M * u);
    std::cout << *cell << ": (" << M.rows() << ", " << M.cols() << ")" << ": "
              << val << std::endl;
    // u is constant and thus its gradient is zero, thus val should be zero.
    EXPECT_NEAR(val, 0.0, 1.0E-10) << "missmatch on cell " << *cell;
  }
}

//*************************************************
// Flux element matrix provider
//*************************************************
TEST(FluxElementMatrixProvider, mass_mat_test) {
  std::cout << ">>> FluxElementMatrixProvider: Test of computation of element "
               "matrices (alpha=1) \n";
  // Build test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Create a trial space involving a flux
  ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
  auto q_n = factory_trial.AddFluxComponent(0);
  auto fe_space_trial = factory_trial.Build();

  // Create a test space of cubic H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_test(mesh_p);
  auto v = factory_test.AddH1Component(3);
  auto fe_space_test = factory_test.Build();

  auto alpha = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> double { return 1.0; });
  FluxElementMatrixProvider provider(fe_space_trial, fe_space_test, q_n, v,
                                     alpha);

  // compute element matrices on all cells
  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);
    std::cout << *cell << ": (" << M.rows() << ", " << M.cols() << ")"
              << std::endl;
    // multiply with vector representin constant test function
    // since we chose q_n to be piecewise constant,
    // the resulting vector should contain +/- the length
    // of the edges.
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);
    Eigen::VectorXd segment_lengths = one_vec_r.transpose() * M;

    for (int k = 0; k < cell->RefEl().NumSubEntities(1); k++) {
      double length =
          lf::geometry::Volume(*(cell->Geometry()->SubGeometry(1, k)));
      EXPECT_NEAR(length, std::abs(segment_lengths(k)), 1.0E-10)
          << "missmatch on cell " << *cell << ", segment " << k;
    }
  }
}

//*************************************************
// Trace element matrix provider
//*************************************************
// The following four tests are based on Gauss divergence theorem:
// Setting u and v to constant one, we expect that
// \int_{\partial K} u v \beta \cdot n ds = \int_K div(\beta) dx
// For functions with zero divergence or constant one divergence, we can
// then verify the provider by comparing the value of 1^T B 1
// with  zero or the volume of the element.
TEST(TraceElementMatrixProvider, mass_mat_test_1) {
  std::cout << ">>> TraceElementMatrixProvider: Test of computation of element "
               "matrices (beta= [1.0, 1.0]) \n";
  // Build test mesh:
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Create a trial space of quadratic  H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
  auto u = factory_trial.AddH1Component(2);
  auto fe_space_trial = factory_trial.Build();

  // Create a test space of cubic H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_test(mesh_p);
  auto v = factory_test.AddH1Component(3);
  auto fe_space_test = factory_test.Build();

  auto beta = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> Eigen::Vector2d {
        return (Eigen::Vector2d(2) << 1.0, 1.0).finished();
      });
  TraceElementMatrixProvider provider(fe_space_trial, fe_space_test, u, v,
                                      beta);

  // compute element matrices on all cells:
  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);
    // multiply with constant one vectors from left and right
    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);
    const double val = one_vec_r.dot(M * one_vec_c);

    std::cout << *cell << ": (" << M.rows() << ", " << M.cols() << ") :" << val
              << std::endl;
    // check using gauss theorem.
    EXPECT_NEAR(val, 0.0, 1E-10) << "missmatch for cell" << *cell;
  }
}

TEST(TraceElementMatrixProvider, mass_mat_test_2) {
  std::cout << ">>> TraceElementMatrixProvider: Test of computation of element "
               "matrices (beta= [x/2, y/2]) \n";
  // Build test mesh:
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Create a trial space of quadratic  H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
  auto u = factory_trial.AddH1Component(2);
  auto fe_space_trial = factory_trial.Build();

  // Create a test space of cubic H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_test(mesh_p);
  auto v = factory_test.AddH1Component(3);
  auto fe_space_test = factory_test.Build();

  auto beta = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> Eigen::Vector2d { return x / 2.0; });
  TraceElementMatrixProvider provider(fe_space_trial, fe_space_test, u, v,
                                      beta);

  // compute element matrices on all cells:
  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);
    // multiply with constant one vectors from left and riths.
    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);
    const double val = one_vec_r.dot(M * one_vec_c);

    std::cout << *cell << ": (" << M.rows() << ", " << M.cols() << ") :" << val
              << std::endl;
    // check using gauss theorem.
    EXPECT_NEAR(val, lf::geometry::Volume(*(cell->Geometry())), 1.0E-10)
        << "missmatch for cell " << *cell;
  }
}

TEST(TraceElementMatrixProvider, mass_mat_test_3) {
  std::cout << ">>> TraceElementMatrixProvider: Test of computation of element "
               "matrices (beta= [x^2*y, -x*y^2]) \n";
  // Build test mesh:
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Create a trial space of quadratic  H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
  auto u = factory_trial.AddH1Component(3);
  auto fe_space_trial = factory_trial.Build();

  // Create a test space of cubic H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_test(mesh_p);
  auto v = factory_test.AddH1Component(3);
  auto fe_space_test = factory_test.Build();

  auto beta = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> Eigen::Vector2d {
        return (Eigen::Vector2d(2) << x[0] * x[0] * x[1], -x[0] * x[1] * x[1])
            .finished();
      });
  TraceElementMatrixProvider provider(fe_space_trial, fe_space_test, u, v,
                                      beta);

  // compute element matrices on all cells:
  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);
    // multiply with constant one vectors from left and right
    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);
    const double val = one_vec_r.dot(M * one_vec_c);

    std::cout << *cell << ": (" << M.rows() << ", " << M.cols() << ") :" << val
              << std::endl;
    // check using gauss theorem
    EXPECT_NEAR(val, 0, 1.0E-8) << "missmatch for cell " << *cell;
  }
}

TEST(TraceElementMatrixProvider, mass_mat_test_4) {
  std::cout << ">>> TraceElementMatrixProvider: Test of computation of element "
               "matrices (alpha= [e^x*cos(y), -e^x*sin(y)]) \n";

  // Build test mesh:
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Create a trial space of quadratic  H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_trial(mesh_p);
  auto u = factory_trial.AddH1Component(3);
  auto fe_space_trial = factory_trial.Build();

  // Create a test space of cubic H1-conforming polynomials
  ProductUniformFESpaceFactory<double> factory_test(mesh_p);
  auto v = factory_test.AddH1Component(3);
  auto fe_space_test = factory_test.Build();

  auto beta = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> Eigen::Vector2d {
        return (Eigen::Vector2d(2) << std::exp(x[0]) * std::cos(x[1]),
                -std::exp(x[0]) * std::sin(x[1]))
            .finished();
      });
  TraceElementMatrixProvider provider(fe_space_trial, fe_space_test, u, v,
                                      beta);

  // compute element matrices on all cells:
  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);
    // multiply with constant one vectors from left and right.
    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);
    const double val = one_vec_r.dot(M * one_vec_c);

    std::cout << *cell << ": (" << M.rows() << ", " << M.cols() << ") :" << val
              << std::endl;
    // check using gauss theorem
    EXPECT_NEAR(val, 0, 1.0E-5) << "missmatch for cell " << *cell;
  }
}

}  // namespace projects::dpg::test
