/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Unit tests for parametric finite element facilities
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <iostream>
#include "lf/fe/fe_testutils.h"
#include "lf/fe/lin_fe.h"
#include "lf/fe/sec_ord_ell_bvp.h"

#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::fe::test {

template <typename FFUNC, typename DIFF_COEFF, typename REAC_COEFF>
std::vector<double> LinFEEnergyTest(int reflevels, FFUNC v, DIFF_COEFF alpha,
                                    REAC_COEFF gamma) {
  // Generate hierarchy of meshes of the unit square
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(
          lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3),
          reflevels);
  // Compute energies of the pw. linear interpolants of v on the different
  // refinement levels
  std::vector<double> energies{EnergiesOfInterpolants<double>(
      *multi_mesh_p, v, alpha, gamma,
      std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>())};
  size_type L = energies.size();
  for (int l = 0; l < L; ++l) {
    std::cout << "Level" << l << ": energy = " << energies[l] << std::endl;
  }
  return energies;
}

TEST(lf_gfe, a_dir_dbg_1) {
  // Four levels of refinement
  const int reflevels = 4;
  std::cout << "#### TEST: Bilinear form test 1, constant reference function"
            << std::endl;
  std::cout << "alpha = 1, gamma =, => energy = 0 " << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return 1.0; };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double { return 1.0; };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double { return 0.0; };

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  Eigen::VectorXd energies =
      Eigen::Map<Eigen::VectorXd>(en.data(), en.size(), 1);
}

TEST(lf_gfe, a_dir_dbg_2) {
  // Four levels of refinement
  const int reflevels = 4;
  std::cout << "#### TEST: Bilinear form test 2, linear reference function"
            << std::endl;
  std::cout << "alpha = 1+x, gamma =, => energy = 7.5 " << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double { return (1.0 + x[0]); };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double { return 0.0; };

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
}

TEST(lf_gfe, a_dir_dbg_3) {
  // Four levels of refinement
  const int reflevels = 4;
  std::cout << "#### TEST: Bilinear form test 3, linear reference function"
            << std::endl;
  std::cout << "alpha = 1+x^2+y^2, gamma =, => energy = 25/3 " << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double {
    return (1.0 + x[0] * x[0] + x[1] * x[1]);
  };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double { return 0.0; };

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
}

TEST(lf_gfe, a_dir_dbg_4) {
  // Four levels of refinement
  const int reflevels = 4;
  std::cout << "#### TEST: Bilinear form test 4, linear reference function"
            << std::endl;
  std::cout << "alpha = 0, gamma = 1, => energy = 8/3 = 2.6666.." << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double { return 0.0; };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double { return 1.0; };

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
}

TEST(lf_gfe, a_dir_dbg_5) {
  // Four levels of refinement
  const int reflevels = 4;
  std::cout << "#### TEST: Bilinear form test 5, linear reference function"
            << std::endl;
  std::cout << "alpha = 0, gamma = 1/(x^2+y^2), => energy = 8/3 = 2.6666.."
            << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double { return 0.0; };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double {
    return (1.0 / (1 + x[0] * x[0] + x[1] * x[1]));
  };

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
}

TEST(lf_gfe, a_dir_dbg_6) {
  // Four levels of refinement
  const int reflevels = 4;
  std::cout << "#### TEST: Bilinear form test 6, f = e^(x*y)" << std::endl;
  std::cout << "alpha = 1, gamma = 0, => energy = 1.59726.." << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double { return 1.0; };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double { return 0.0; };

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
}

TEST(lf_gfe, a_dir_dbg_7) {
  // Four levels of refinement
  const int reflevels = 4;
  std::cout << "#### TEST: Bilinear form test 7, f = e^(x*y)" << std::endl;
  std::cout << "alpha = 1+x^2+y^2, gamma = 0, => energy = 3.39057.."
            << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double {
    return (1.0 + x[0] * x[0] + x[1] * x[1]);
  };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double { return 0.0; };

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
}

TEST(lf_gfe, a_dir_dbg_8) {
  // Four levels of refinement
  const int reflevels = 4;
  std::cout << "#### TEST: Bilinear form test 8, f = e^(x*y)" << std::endl;
  std::cout << "alpha = 1+x^2+y^2, gamma = 1/(x^2+y^2), => energy = 3.39057.."
            << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double {
    return (1.0 + x[0] * x[0] + x[1] * x[1]);
  };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double {
    return 1.0 / (1.0 + x[0] * x[0] + x[1] * x[1]);
  };

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
}

TEST(lf_gfe, Neu_BVP_ass) {
  std::cout << "#### TEST: Assembly for Neumann BVP ###" << std::endl;
  // Right-hand side source function
  auto f = [](Eigen::Vector2d x) -> double {
    return (std::sin(2 * M_PI * x[0]) * std::sin(2 * M_PI * x[1]));
  };
  // Neumann data
  auto h = [](Eigen::Vector2d /*x*/) -> double { return 1.0; };

  // Class describing the Neumann boundary value problem
  std::shared_ptr<PureNeumannProblemLaplacian<decltype(f), decltype(h)>> bvp_p =
      std::make_shared<PureNeumannProblemLaplacian<decltype(f), decltype(h)>>(
          f, h);

  // Obtain finite Element space
  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Set up global FE space; lowest order Lagrangian finite elements
  LinearLagrangianFESpace fe_space(mesh_p);

  // Generate linear system
  auto [A, phi] = SecOrdEllBVPLagrFELinSys<double>(fe_space, bvp_p);

  Eigen::MatrixXd A_dense = A;
  std::cout << "Galerkin matrix = " << A_dense << std::endl;
  std::cout << "rhs. vector = " << phi.transpose() << std::endl;
  std::cout << "Row sums = " << A_dense.rowwise().sum().transpose()
            << std::endl;
  std::cout << "Col sums = " << A_dense.colwise().sum() << std::endl;
}

}  // namespace lf::fe::test
