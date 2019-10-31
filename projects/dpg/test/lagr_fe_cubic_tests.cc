/**
 * @file
 * @brief local and global tests for the implementation of cubic
 * Lagrangian finite element local shape functions.
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */
#include <gtest/gtest.h>
#include <iostream>

#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#include "../lagr_fe_cubic.h"
#include "lagr_test_utils.h"

#define REFLEV 4

namespace projects::dpg::test {

TEST(FeLagrO3, scal_fe_coeff_node) {
  std::cout << ">>> O3 FE: Test of consistency of nodal interpolation \n";

  FeLagrangeO3Tria<double> tfe{};
  EXPECT_TRUE(scalarFEEvalNodeTest(tfe));

  FeLagrangeO3Quad<double> qfe{};
  EXPECT_TRUE(scalarFEEvalNodeTest(qfe));

  FeLagrangeO3Segment<double> sfe{};
  EXPECT_TRUE(scalarFEEvalNodeTest(sfe));
}

TEST(FeLagrO3, scal_fe_val_node) {
  std::cout << ">>> O3 FE: Test of consistency of nodal interpolation \n";

  FeLagrangeO3Tria<double> tfe{};
  EXPECT_TRUE(scalarFEInterpTest(tfe));

  FeLagrangeO3Quad<double> qfe{};
  EXPECT_TRUE(scalarFEInterpTest(qfe));

  FeLagrangeO3Segment<double> sfe{};
  EXPECT_TRUE(scalarFEInterpTest(sfe));
}

TEST(FeLagrO3, shape_function_information) {
  std::cout << ">>> O3 Fe: Test of number of reference shape functions \n";

  FeLagrangeO3Tria<double> tfe{};
  EXPECT_EQ(tfe.NumRefShapeFunctions(), 10);
  EXPECT_EQ(tfe.NumRefShapeFunctions(0, 0), 1);
  EXPECT_EQ(tfe.NumRefShapeFunctions(1, 0), 2);
  EXPECT_EQ(tfe.NumRefShapeFunctions(2, 0), 1);

  FeLagrangeO3Quad<double> qfe{};
  EXPECT_EQ(qfe.NumRefShapeFunctions(), 16);
  EXPECT_EQ(qfe.NumRefShapeFunctions(0, 0), 4);
  EXPECT_EQ(qfe.NumRefShapeFunctions(1, 0), 2);
  EXPECT_EQ(qfe.NumRefShapeFunctions(2, 0), 1);

  FeLagrangeO3Segment<double> sfe{};
  EXPECT_EQ(sfe.NumRefShapeFunctions(), 4);
  EXPECT_EQ(sfe.NumRefShapeFunctions(0, 0), 2);
  EXPECT_EQ(sfe.NumRefShapeFunctions(1, 0), 1);
}

TEST(FeLagrO3, shape_function_cardinality) {
  std::cout << ">>> O3 Fe: Test of cardinal basis property of reference shape "
               "functions \n";

  FeLagrangeO3Tria<double> tfe{};
  EXPECT_TRUE(tfe.EvalReferenceShapeFunctions(tfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(10, 10)));

  FeLagrangeO3Quad<double> qfe{};
  EXPECT_TRUE(qfe.EvalReferenceShapeFunctions(qfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(16, 16)));

  FeLagrangeO3Segment<double> sfe{};
  EXPECT_TRUE(sfe.EvalReferenceShapeFunctions(sfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(4, 4)));
}

TEST(FeLagrO3, shape_function_gradient_sum) {
  std::cout << ">>> O3 Fe: Test of gadients summing up to 0 \n";

  FeLagrangeO3Tria<double> tfe{};
  EXPECT_TRUE(scalarFEEvalGRadTest(tfe));

  FeLagrangeO3Quad<double> qfe{};
  EXPECT_TRUE(scalarFEEvalGRadTest(qfe));

  FeLagrangeO3Segment<double> sfe{};
  EXPECT_TRUE(scalarFEEvalGRadTest(sfe));
}

// lehrfem tests for the ReactionDiffusionElementMatrixProvider
// are used to test if the computation of the local shape functions and
// their gradients are correct.
TEST(FeLagrO3, mass_mat_test) {
  std::cout << " >>> O3 FE: Test of computation of element matrices \n";

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto fe_space = std::make_shared<lf::uscalfe::UniformScalarFESpace<double>>(
      mesh_p, rfs_tria, rfs_quad);

  // coefficients:
  auto alpha = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> double { return 1.0; });
  auto gamma = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> double { return 1.0; });

  lf::uscalfe::quad_rule_collection_t quad_rules{
      {lf::base::RefEl::kTria(),
       lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 6)},
      {lf::base::RefEl::kQuad(),
       lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 6)}};

  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                                      decltype(gamma)>
      provider(fe_space, alpha, gamma, quad_rules);

  for (const lf::mesh::Entity* const cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);

    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);

    const double vol = one_vec_r.dot(M * one_vec_c);
    EXPECT_NEAR(vol, lf::geometry::Volume(*(cell->Geometry())), 1.0E-10)
        << "missmatch for cell " << cell;
  }
}

//***************************************************************************
// full galerkin tests: Run some of the full galerkin test of lf::uscalfe::test
//*****************************************************************************

TEST(FeLagrO3, a_dir_dbg_1) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Bilinear form test 1, constant reference function"
            << std::endl;
  std::cout << "alpha = 1, gamma = 0, => energy = 0 " << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  // Diffusion coefficient
  auto alpha = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  // Reaction coefficient
  auto gamma = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();

  auto en{FEEnergyTest(reflevels, v, alpha, gamma, rfs_tria, rfs_quad)};
  Eigen::VectorXd energies =
      Eigen::Map<Eigen::VectorXd>(en.data(), en.size(), 1);
  EXPECT_NEAR(energies.norm(), 0.0, 1.0E-10);
}

TEST(FeLagrO3, a_dir_dbg_2) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Bilinear form test 2, linear reference function"
            << std::endl;
  std::cout << "alpha = 1+x, gamma = 0, => energy = 7.5 " << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Diffusion coefficient
  auto alpha = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (1.0 + x[0]); });
  // Reaction coefficient
  auto gamma = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto en{FEEnergyTest(reflevels, v, alpha, gamma, rfs_tria, rfs_quad)};
  EXPECT_NEAR(en.back(), 7.5, 1.0E-4);
}

TEST(FeLagrO3, a_dir_dbg_3) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Bilinear form test 3, linear reference function"
            << std::endl;
  std::cout << "alpha = 1+x^2+y^2, gamma = 0, => energy = 25/3 = 8.333.. "
            << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Diffusion coefficient
  auto alpha = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return (1.0 + x[0] * x[0] + x[1] * x[1]);
  });
  // Reaction coefficient
  auto gamma = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto en{FEEnergyTest(reflevels, v, alpha, gamma, rfs_tria, rfs_quad)};
  EXPECT_NEAR(en.back(), 25.0 / 3.0, 1.0E-4);
}

TEST(FeLagrO3, a_dir_dbg_4) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Bilinear form test 4, linear reference function"
            << std::endl;
  std::cout << "alpha = 0, gamma = 1, => energy = 8/3 = 2.6666.." << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Diffusion coefficient
  auto alpha = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });
  // Reaction coefficient
  auto gamma = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto en{FEEnergyTest(reflevels, v, alpha, gamma, rfs_tria, rfs_quad)};
  EXPECT_NEAR(en.back(), 8.0 / 3.0, 1.0E-4);
}

TEST(FeLagrO3, a_dir_dbg_5) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Bilinear form test 5, linear reference function"
            << std::endl;
  std::cout << "alpha = 0, gamma = 1/(x^2+y^2), => energy = 1.42447.."
            << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Diffusion coefficient
  auto alpha = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });
  // Reaction coefficient
  auto gamma = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return (1.0 / (1 + x[0] * x[0] + x[1] * x[1]));
  });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto en{FEEnergyTest(reflevels, v, alpha, gamma, rfs_tria, rfs_quad)};
  EXPECT_NEAR(en.back(), 1.42447, 1.0E-4);
}

TEST(FeLagrO3, a_dir_dbg_6) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Bilinear form test 6, f = e^(x*y)" << std::endl;
  std::cout << "alpha = 1, gamma = 0, => energy = 1.59726.." << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  // Diffusion coefficient
  auto alpha = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  // Reaction coefficient
  auto gamma = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto en{FEEnergyTest(reflevels, v, alpha, gamma, rfs_tria, rfs_quad)};
  EXPECT_NEAR(en.back(), 1.59726, 2.0E-3);
}

TEST(FeLagrO3, a_dir_dbg_7) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Bilinear form test 7, f = e^(x*y)" << std::endl;
  std::cout << "alpha = 1+x^2+y^2, gamma = 0, => energy = 3.39057.."
            << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  // Diffusion coefficient
  auto alpha = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return (1.0 + x[0] * x[0] + x[1] * x[1]);
  });
  // Reaction coefficient
  auto gamma = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto en{FEEnergyTest(reflevels, v, alpha, gamma, rfs_tria, rfs_quad)};
  EXPECT_NEAR(en.back(), 3.39057, 2.0E-3);
}

TEST(FeLagrO3, a_dir_dbg_8) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Bilinear form test 8, f = e^(x*y)" << std::endl;
  std::cout << "alpha = 1+x^2+y^2, gamma = 1/(x^2+y^2), => energy = 4.44757.."
            << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  // Diffusion coefficient
  auto alpha = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return (1.0 + x[0] * x[0] + x[1] * x[1]);
  });
  // Reaction coefficient
  auto gamma = lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return 1.0 / (1.0 + x[0] * x[0] + x[1] * x[1]);
  });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto en{FEEnergyTest(reflevels, v, alpha, gamma, rfs_tria, rfs_quad)};
  EXPECT_NEAR(en.back(), 4.44757, 2.0E-3);
}

TEST(FeLagrO3, a_dir_dbg_9) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Bilinear form test 9, linear f" << std::endl;
  std::cout << "alpha = [1+x_2,0;0,1+x_1], gamma = 0, => energy = 7.5.."
            << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Diffusion coefficient
  auto alpha =
      lf::uscalfe::MeshFunctionGlobal([](Eigen::Vector2d x) -> Eigen::Matrix2d {
        return (Eigen::Matrix2d(2, 2) << 1.0 + x[1], 0.0, 0.0, 1.0 + x[0])
            .finished();
      });
  // Reaction coefficient
  auto gamma = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto en{FEEnergyTest(reflevels, v, alpha, gamma, rfs_tria, rfs_quad)};
  EXPECT_NEAR(en.back(), 7.5, 2.0E-3);
}

//*******************************************************************
// full galerkin tests: run the lf::uscalfe::test boundary contributions
// assembly tests
//

TEST(FeLagrO3, bd_a_1) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Boundary energy test 1, linear f" << std::endl;

  std::cout << "eta = 1 => energy = 32/3 " << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Impedance coefficient
  auto eta = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto rfs_segment = std::make_shared<FeLagrangeO3Segment<double>>();
  auto en{FEInterfaceEnergyTest(reflevels, v, eta, rfs_tria, rfs_quad,
                                rfs_segment)};
  EXPECT_NEAR(en.back(), 32.0 / 3.0, 1.0E-4);
}

TEST(FeLagrO3, bd_a_2) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << ">>> O3 FE: Boundary energy test 2, linear f" << std::endl;

  std::cout << "eta = x^2+y^2 => energy = 46/3" << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Impedance coefficient
  auto eta = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (x[0] * x[0] + x[1] * x[1]); });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto rfs_segment = std::make_shared<FeLagrangeO3Segment<double>>();
  auto en{FEInterfaceEnergyTest(reflevels, v, eta, rfs_tria, rfs_quad,
                                rfs_segment)};
  EXPECT_NEAR(en.back(), 46.0 / 3.0, 1.0E-4);
}

TEST(FeLagrO3, bd_a_3) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "O3 FE: Boundary energy test 3, f = exp(x*y)" << std::endl;

  std::cout << "eta = x^2+y^2 => energy = 9.58358.." << std::endl;
  // Synthetic function
  auto v = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  // Impedance coefficient
  auto eta = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (x[0] * x[0] + x[1] * x[1]); });

  auto rfs_tria = std::make_shared<FeLagrangeO3Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO3Quad<double>>();
  auto rfs_segment = std::make_shared<FeLagrangeO3Segment<double>>();
  auto en{FEInterfaceEnergyTest(reflevels, v, eta, rfs_tria, rfs_quad,
                                rfs_segment)};
  EXPECT_NEAR(en.back(), 9.58358, 1.0E-3);
}

}  // namespace projects::dpg::test
