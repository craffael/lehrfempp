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

#define REFLEV 5

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
  const int reflevels = REFLEV;
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
  EXPECT_NEAR(energies.norm(), 0.0, 1.0E-10);
}

TEST(lf_gfe, a_dir_dbg_2) {
  // Four levels of refinement
  const int reflevels = REFLEV;
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
  EXPECT_NEAR(en.back(), 7.5, 1.0E-4);
}

TEST(lf_gfe, a_dir_dbg_3) {
  // Four levels of refinement
  const int reflevels = REFLEV;
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
  EXPECT_NEAR(en.back(), 25.0 / 3.0, 1.0E-4);
}

TEST(lf_gfe, a_dir_dbg_4) {
  // Four levels of refinement
  const int reflevels = REFLEV;
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
  auto v = [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double { return 0.0; };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double {
    return (1.0 / (1 + x[0] * x[0] + x[1] * x[1]));
  };

  auto en{LinFEEnergyTest(reflevels, v, alpha, gamma)};
  EXPECT_NEAR(en.back(), 1.42447, 1.0E-4);
}

TEST(lf_gfe, a_dir_dbg_6) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 6, f = e^(x*y)" << std::endl;
  std::cout << "alpha = 1, gamma = 0, => energy = 1.59726.." << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double { return 1.0; };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double { return 0.0; };

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
  auto v = [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> double {
    return (1.0 + x[0] * x[0] + x[1] * x[1]);
  };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double { return 0.0; };

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
  EXPECT_NEAR(en.back(), 4.44757, 2.0E-3);
}

TEST(lf_gfe, a_dir_dbg_9) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "#### TEST: Bilinear form test 9, linear f" << std::endl;
  std::cout << "alpha = [1+x_2,0;0,1+x_1], gamma = 0, => energy = 7.5.."
            << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); };
  // Diffusion coefficient
  auto alpha = [](Eigen::Vector2d x) -> Eigen::Matrix2d {
    return (Eigen::Matrix2d(2, 2) << 1.0 + x[1], 0.0, 0.0, 1.0 + x[0])
        .finished();
  };
  // Reaction coefficient
  auto gamma = [](Eigen::Vector2d x) -> double { return 0.0; };

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
      *multi_mesh_p, v, eta, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>(),
      std::make_shared<SegmentLinearLagrangeFE<double>>(), edge_sel)};
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
  auto v = [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); };
  // Impedance coefficient
  auto eta = [](Eigen::Vector2d x) -> double { return 1.0; };

  auto en{LinFEInterfaceEnergyTest(reflevels, v, eta)};
  EXPECT_NEAR(en.back(), 32.0 / 3.0, 1.0E-4);
}

TEST(lf_full, bd_a_2) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "<<<< TEST: Boundary energy test 2, linear f >>>>" << std::endl;

  std::cout << "eta = x^2+y^2 => energy = 46/3" << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); };
  // Impedance coefficient
  auto eta = [](Eigen::Vector2d x) -> double {
    return (x[0] * x[0] + x[1] * x[1]);
  };

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
  auto v = [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); };
  // Impedance coefficient
  auto eta = [](Eigen::Vector2d x) -> double {
    return (x[0] * x[0] + x[1] * x[1]);
  };

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
      *multi_mesh_p, v, f, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>())};

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
  auto v = [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); };
  // Right hand side source function
  auto f = [](Eigen::Vector2d x) -> double { return 1.0; };

  auto en{LinFERHSTest(reflevels, v, f)};
  EXPECT_NEAR(en.back(), 1.5, 1.0E-4);
}

TEST(lf_full, ell_2) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "<<<< TEST: RHS test 2, v = exp(x*y) >>>>" << std::endl;

  std::cout << "f = x/(1+y) => fval = " << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); };
  // Right hand side source function
  auto f = [](Eigen::Vector2d x) -> double { return (x[0] / (1.0 + x[1])); };

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
      *multi_mesh_p, v, f, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>(),
      std::make_shared<SegmentLinearLagrangeFE<double>>(), edge_sel)};
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
  auto v = [](Eigen::Vector2d x) -> double { return (2.0 * x[0] + x[1]); };
  // Right hand side source function
  auto f = [](Eigen::Vector2d x) -> double { return 1.0; };

  auto en{LinFEBoundaryRHSTest(reflevels, v, f)};
  EXPECT_NEAR(en.back(), 4.5, 1.0E-4);
}

TEST(lf_full, bd_rhs_2) {
  // Four levels of refinement
  const int reflevels = REFLEV;
  std::cout << "<<<< TEST: Boundary RHS test 1, v  = exp(x*y)>>>>" << std::endl;

  std::cout << "h = x/(1+y) => fval = 1.62539" << std::endl;
  // Synthetic function
  auto v = [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); };
  // Right hand side source function
  auto f = [](Eigen::Vector2d x) -> double { return (x[0] / (1.0 + x[1])); };

  auto en{LinFEBoundaryRHSTest(reflevels, v, f)};
  EXPECT_NEAR(en.back(), 1.62539, 1.0E-4);
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
  LinearLagrangianFESpace<double> fe_space(mesh_p);

  // Generate linear system
  auto [A, phi] = SecOrdEllBVPLagrFELinSys<double>(fe_space, bvp_p);

  Eigen::MatrixXd A_dense = A;
  std::cout << "Galerkin matrix = " << A_dense << std::endl;
  std::cout << "rhs. vector = " << phi.transpose() << std::endl;
  EXPECT_NEAR(A_dense.rowwise().sum().norm(), 0.0, 1.0E-10)
      << "Row sum not zero";
  EXPECT_NEAR(A_dense.colwise().sum().norm(), 0.0, 1.0E-10)
      << "Col sum not zero";
}

// **********************************************************************
// Test: Solution of BVP
// **********************************************************************

template <typename SOLFUNC, typename SOLGRAD>
std::vector<std::pair<double, double>> TestConvergenceEllBVPFESol(
    int reflevels, std::shared_ptr<const SecondOrderEllipticBVP<double>> bvp_p,
    SOLFUNC solution, SOLGRAD sol_grad) {
  std::vector<std::pair<double, double>> errnorms{};

  // Generate hierarchy of meshes of the unit square
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(
          lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3),
          reflevels);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};
  // Number of levels
  size_type L = multi_mesh.NumLevels();

  // Do computations on all levels
  for (int l = 0; l < L; ++l) {
    std::shared_ptr<const mesh::Mesh> mesh_p{multi_mesh.getMesh(l)};
    // Set up global FE space; lowest order Lagrangian finite elements
    LinearLagrangianFESpace<double> fe_space(mesh_p);
    const lf::assemble::DofHandler& dofh{fe_space.LocGlobMap()};

    // Compute Galerkin matrix A and right-hand-side vector
    auto [A, phi] = SecOrdEllBVPLagrFELinSys<double>(fe_space, bvp_p);
    // Solve linear system
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd sol_vec = solver.solve(phi);
    // Compute error norms
    // Helper class for L2 error computation
    LocalL2NormDifference lc_L2(
        fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()),
        fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()), solution);
    // Helper class for H1 semi norm
    LocL2GradientFEDifference lc_H1(
        fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()),
        fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()), sol_grad);

    double L2err = NormOfDifference(dofh, lc_L2, sol_vec);
    double H1serr = NormOfDifference(dofh, lc_H1, sol_vec);
    errnorms.emplace_back(L2err, H1serr);
  }

  return errnorms;
}

// Definition of boundary value problem

class MixedEllipticBVP : public SecondOrderEllipticBVP<double> {
 public:
  MixedEllipticBVP()
      : edge_sel_(
            [](Eigen::Vector2d x) -> bool { return ((x[0] + x[1]) < 1.0); }) {}
  MixedEllipticBVP(const MixedEllipticBVP&) = default;
  MixedEllipticBVP(MixedEllipticBVP&&) noexcept = default;
  MixedEllipticBVP& operator=(const MixedEllipticBVP&) = default;
  MixedEllipticBVP& operator=(MixedEllipticBVP&&) noexcept = default;

  Eigen::Matrix<double, 2, 2> alpha(Eigen::Vector2d x) const override {
    return (1 + x.squaredNorm()) * Eigen::Matrix<double, 2, 2>::Identity(2, 2);
  }
  double gamma(Eigen::Vector2d x) const override { return 0.0; }
  double eta(Eigen::Vector2d x) const override { return 0.0; }
  bool EssentialConditionsOnEdge(const lf::mesh::Entity& edge) const override {
    return edge_sel_(edge);
  }
  bool IsImpedanceEdge(const lf::mesh::Entity& /*edge*/) const override {
    return false;
  }
  double f(Eigen::Vector2d /*x*/) const override { return 2.0; }
  double g(Eigen::Vector2d x) const override {
    return 0.5 * std::log(x.squaredNorm() + 1);
  }
  double h(Eigen::Vector2d x) const override {
    if (x[0] > x[1]) {
      return x[0];
    } else {
      return x[1];
    }
  }

 private:
  /// Selector for Dirichlet boundary conditions
  lf::refinement::EntityCenterPositionSelector<
      std::function<bool(Eigen::Vector2d)>>
      edge_sel_;
};

TEST(lfe_bvp, bvp_1) {
  // The solution and its gradient
  auto u = [](Eigen::Vector2d x) -> double {
    return (0.5 * std::log(x.squaredNorm() + 1.0));
  };
  auto grad_u = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return (1.0 / (x.squaredNorm() + 1.0)) * x;
  };
  auto bvp_p = std::make_shared<MixedEllipticBVP>();
  auto errnorms = TestConvergenceEllBVPFESol(4, bvp_p, u, grad_u);
}

}  // namespace lf::fe::test
