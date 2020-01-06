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
#include <lf/uscalfe/uscalfe.h>
#include <iostream>
#include "fe_testutils.h"

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>

#define REFLEV 6

namespace lf::uscalfe::test {

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
      std::make_shared<FeLagrangeO1Tria<double>>(),
      std::make_shared<FeLagrangeO1Quad<double>>())};
  size_type L = energies.size();
  for (int l = 0; l < L; ++l) {
    std::cout << "Level" << l << ": energy = " << energies[l] << std::endl;
  }
  return energies;
}

TEST(lf_gfe, a_dir_dbg_1) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 1, constant reference function"
            << std::endl;
  std::cout << "alpha = 1, gamma =, => energy = 0 " << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  // Diffusion coefficient
  auto alpha = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  // Reaction coefficient
  auto gamma = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  Eigen::VectorXd energies =
      Eigen::Map<Eigen::VectorXd>(en.data(), en.size(), 1);
  EXPECT_NEAR(energies.norm(), 0.0, 1.0E-10);
}

TEST(lf_gfe, a_dir_dbg_2) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 2, linear reference function"
            << std::endl;
  std::cout << "alpha = 1+x, gamma =, => energy = 7.5 " << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Diffusion coefficient
  auto alpha = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (1.0 + x[0]); });
  // Reaction coefficient
  auto gamma = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  EXPECT_NEAR(en.back(), 7.5, 1.0E-4);
}

TEST(lf_gfe, a_dir_dbg_3) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 3, linear reference function"
            << std::endl;
  std::cout << "alpha = 1+x^2+y^2, gamma =, => energy = 25/3 " << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Diffusion coefficient
  auto alpha = mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return (1.0 + x[0] * x[0] + x[1] * x[1]);
  });
  // Reaction coefficient
  auto gamma = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  EXPECT_NEAR(en.back(), 25.0 / 3.0, 1.0E-4);
}

TEST(lf_gfe, a_dir_dbg_4) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 4, linear reference function"
            << std::endl;
  std::cout << "alpha = 0, gamma = 1, => energy = 8/3 = 2.6666.." << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Diffusion coefficient
  auto alpha = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });
  // Reaction coefficient
  auto gamma = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  EXPECT_NEAR(en.back(), 8.0 / 3.0, 1.0E-4);
}

TEST(lf_gfe, a_dir_dbg_5) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 5, linear reference function"
            << std::endl;
  std::cout << "alpha = 0, gamma = 1/(x^2+y^2), => energy = 1.42447.."
            << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Diffusion coefficient
  auto alpha = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });
  // Reaction coefficient
  auto gamma = mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return (1.0 / (1 + x[0] * x[0] + x[1] * x[1]));
  });

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  EXPECT_NEAR(en.back(), 1.42447, 1.0E-4);
}

TEST(lf_gfe, a_dir_dbg_6) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 6, f = e^(x*y)" << std::endl;
  std::cout << "alpha = 1, gamma = 0, => energy = 1.59726.." << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  // Diffusion coefficient
  auto alpha = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  // Reaction coefficient
  auto gamma = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  EXPECT_NEAR(en.back(), 1.59726, 2.0E-3);
}

TEST(lf_gfe, a_dir_dbg_7) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 7, f = e^(x*y)" << std::endl;
  std::cout << "alpha = 1+x^2+y^2, gamma = 0, => energy = 3.39057.."
            << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  // Diffusion coefficient
  auto alpha = mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return (1.0 + x[0] * x[0] + x[1] * x[1]);
  });
  // Reaction coefficient
  auto gamma = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  EXPECT_NEAR(en.back(), 3.39057, 2.0E-3);
}

TEST(lf_gfe, a_dir_dbg_8) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 8, f = e^(x*y)" << std::endl;
  std::cout << "alpha = 1+x^2+y^2, gamma = 1/(x^2+y^2), => energy = 4.44757.."
            << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  // Diffusion coefficient
  auto alpha = mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return (1.0 + x[0] * x[0] + x[1] * x[1]);
  });
  // Reaction coefficient
  auto gamma = mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return 1.0 / (1.0 + x[0] * x[0] + x[1] * x[1]);
  });

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  EXPECT_NEAR(en.back(), 4.44757, 2.0E-3);
}

TEST(lf_gfe, a_dir_dbg_9) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 9, linear f" << std::endl;
  std::cout << "alpha = [1+x_2,0;0,1+x_1], gamma = 0, => energy = 7.5.."
            << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Diffusion coefficient
  auto alpha =
      mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) -> Eigen::Matrix2d {
        return (Eigen::Matrix2d(2, 2) << 1.0 + x[1], 0.0, 0.0, 1.0 + x[0])
            .finished();
      });
  // Reaction coefficient
  auto gamma = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  EXPECT_NEAR(en.back(), 7.5, 2.0E-3);
}

// **********************************************************************
// Tests of assembly of boundary contributions
//

template <typename FFUNC, typename IMP_COEFF>
std::vector<double> LinFEInterfaceEnergyTest(int reflevels, FFUNC v,
                                             IMP_COEFF eta) {
  // Generate hierarchy of meshes of the unit square
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(
          lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3),
          reflevels);
  // Select edges on top and right side of unit square
  auto pos_pred = [](Eigen::Vector2d x) -> bool {
    return ((x[0] + x[1]) > 1.0);
  };
  lf::refinement::EntityCenterPositionSelector edge_sel{pos_pred};

  // Compute energies of the pw. linear interpolants of v on the different
  // refinement levels
  std::vector<double> energies{BoundaryEnergiesOfInterpolants<double>(
      *multi_mesh_p, v, eta, std::make_shared<FeLagrangeO1Tria<double>>(),
      std::make_shared<FeLagrangeO1Quad<double>>(),
      std::make_shared<FeLagrangeO1Segment<double>>(), edge_sel)};
  // Output of energies
  size_type L = energies.size();
  for (int l = 0; l < L; ++l) {
    std::cout << "Level" << l << ": energy = " << energies[l] << std::endl;
  }
  return energies;
}

TEST(lf_full, bd_a_1) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "<<<< TEST: Boundary energy test 1, linear f >>>>" << std::endl;

  std::cout << "eta = 1 => energy = 32/3 " << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Impedance coefficient
  auto eta = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto en{LinFEInterfaceEnergyTest(reflevels, v, eta)};
  EXPECT_NEAR(en.back(), 32.0 / 3.0, 1.0E-4);
}

TEST(lf_full, bd_a_2) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "<<<< TEST: Boundary energy test 2, linear f >>>>" << std::endl;

  std::cout << "eta = x^2+y^2 => energy = 46/3" << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Impedance coefficient
  auto eta = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (x[0] * x[0] + x[1] * x[1]); });

  auto en{LinFEInterfaceEnergyTest(reflevels, v, eta)};
  EXPECT_NEAR(en.back(), 46.0 / 3.0, 1.0E-4);
}

TEST(lf_full, bd_a_3) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "<<<< TEST: Boundary energy test 3, f = exp(x*y) >>>>"
            << std::endl;

  std::cout << "eta = x^2+y^2 => energy = 9.58358.." << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  // Impedance coefficient
  auto eta = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (x[0] * x[0] + x[1] * x[1]); });

  auto en{LinFEInterfaceEnergyTest(reflevels, v, eta)};
  EXPECT_NEAR(en.back(), 9.58358, 1.0E-3);
}

// **********************************************************************
// Tests for load vector
// **********************************************************************

template <typename FFUNC, typename SOURCE_FUNC>
std::vector<double> LinFERHSTest(int reflevels, FFUNC v, SOURCE_FUNC f) {
  // Generate hierarchy of meshes of the unit square
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(
          lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3),
          reflevels);
  // Compute energies of the pw. linear interpolants of v on the different
  // refinement levels
  std::vector<double> fvals{RHSFunctionalForInterpolants<double>(
      *multi_mesh_p, v, f, std::make_shared<FeLagrangeO1Tria<double>>(),
      std::make_shared<FeLagrangeO1Quad<double>>())};

  size_type L = fvals.size();
  for (int l = 0; l < L; ++l) {
    std::cout << "Level" << l << ": fval = " << fvals[l] << std::endl;
  }
  return fvals;
}

TEST(lf_full, ell_1) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "<<<< TEST: RHS test 1, linear v >>>>" << std::endl;

  std::cout << "f = 1 => fval = 1.5 " << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Right hand side source function
  auto f = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto en{LinFERHSTest(reflevels, v, f)};
  EXPECT_NEAR(en.back(), 1.5, 1.0E-4);
}

TEST(lf_full, ell_2) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "<<<< TEST: RHS test 2, v = exp(x*y) >>>>" << std::endl;

  std::cout << "f = x/(1+y) => fval = " << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  // Right hand side source function
  auto f = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (x[0] / (1.0 + x[1])); });

  auto en{LinFERHSTest(reflevels, v, f)};
  EXPECT_NEAR(en.back(), 0.478559, 1.0E-4);
}

// **********************************************************************
// Tests for boundary contributions to RHS functional
// **********************************************************************

template <typename FFUNC, typename SOURCE_FUNC>
std::vector<double> LinFEBoundaryRHSTest(int reflevels, FFUNC v,
                                         SOURCE_FUNC f) {
  // Generate hierarchy of meshes of the unit square
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(
          lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3),
          reflevels);
  // Select edges on top and right side of unit square
  auto pos_pred = [](Eigen::Vector2d x) -> bool {
    return ((x[0] + x[1]) > 1.0);
  };
  lf::refinement::EntityCenterPositionSelector edge_sel{pos_pred};

  // Compute energies of the pw. linear interpolants of v on the different
  // refinement levels
  std::vector<double> fvals{RHSBoundaryFunctionalForInterpolants<double>(
      *multi_mesh_p, v, f, std::make_shared<FeLagrangeO1Tria<double>>(),
      std::make_shared<FeLagrangeO1Quad<double>>(),
      std::make_shared<FeLagrangeO1Segment<double>>(), edge_sel)};
  // Output of energies
  size_type L = fvals.size();
  for (int l = 0; l < L; ++l) {
    std::cout << "Level" << l << ": energy = " << fvals[l] << std::endl;
  }
  return fvals;
}

TEST(lf_full, bd_rhs_1) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "<<<< TEST: Boundary RHS test 1, linear v >>>>" << std::endl;

  std::cout << "h = 1 => fval = 4.5 " << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); });
  // Right hand side source function
  auto f = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto en{LinFEBoundaryRHSTest(reflevels, v, f)};
  EXPECT_NEAR(en.back(), 4.5, 1.0E-4);
}

TEST(lf_full, bd_rhs_2) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "<<<< TEST: Boundary RHS test 1, v  = exp(x*y)>>>>" << std::endl;

  std::cout << "h = x/(1+y) => fval = 1.62539" << std::endl;
  // Synthetic function
  auto v = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  // Right hand side source function
  auto f = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (x[0] / (1.0 + x[1])); });

  auto en{LinFEBoundaryRHSTest(reflevels, v, f)};
  EXPECT_NEAR(en.back(), 1.62539, 1.0E-4);
}

}  // namespace lf::uscalfe::test
