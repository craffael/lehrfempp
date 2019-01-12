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
#include <lf/uscalfe/lagrfe.h>
#include <iostream>
#include "fe_testutils.h"

#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"
#include "sec_ord_ell_bvp.h"

#define REFLEV 6

namespace lf::uscalfe::test {

TEST(lf_bvpfe, Neu_BVP_ass) {
  std::cout << "#### TEST: Assembly for Neumann BVP ###" << std::endl;
  // Right-hand side source function
  auto f = [](Eigen::Vector2d x) -> double {
    return (std::sin(2 * base::kPi * x[0]) * std::sin(2 * base::kPi * x[1]));
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
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

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
  std::cout << multi_mesh;
  // Number of levels
  size_type L = multi_mesh.NumLevels();

  // Do computations on all levels
  for (int l = 0; l < L; ++l) {
    std::shared_ptr<const mesh::Mesh> mesh_p{multi_mesh.getMesh(l)};
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);
    const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

    // Compute Galerkin matrix A and right-hand-side vector
    auto [A, phi] = SecOrdEllBVPLagrFELinSys<double>(fe_space, bvp_p);
    // Solve linear system
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd sol_vec = solver.solve(phi);
    // Compute error norms
    // Helper class for L2 error computation
    MeshFunctionL2NormDifference lc_L2(fe_space, solution, 2);
    // Helper class for H1 semi norm
    MeshFunctionL2GradientDifference lc_H1(fe_space, sol_grad, 2);

    double L2err = NormOfDifference(dofh, lc_L2, sol_vec);
    double H1serr = NormOfDifference(dofh, lc_H1, sol_vec);
    std::cout << "Level " << l << " : errors, L2 = " << L2err
              << ", H1 = " << H1serr << std::endl;
    errnorms.emplace_back(L2err, H1serr);
  }

  return errnorms;
}

// Definition of boundary value problem

class MixedEllipticBVP : public SecondOrderEllipticBVP<double> {
 public:
  MixedEllipticBVP()
      : edge_sel_dir_(
            [](Eigen::Vector2d x) -> bool { return (x[0] < 1.0E-7); }),
        edge_sel_imp_(
            [](Eigen::Vector2d x) -> bool { return (x[1] < 1.0E-7); }) {}
  MixedEllipticBVP(const MixedEllipticBVP&) = default;
  MixedEllipticBVP(MixedEllipticBVP&&) noexcept = default;
  MixedEllipticBVP& operator=(const MixedEllipticBVP&) = default;
  MixedEllipticBVP& operator=(MixedEllipticBVP&&) noexcept = default;

  Eigen::Matrix<double, 2, 2> alpha(Eigen::Vector2d x) const override {
    return (1 + x.squaredNorm()) * Eigen::Matrix<double, 2, 2>::Identity(2, 2);
  }
  double gamma(Eigen::Vector2d x) const override { return 0.0; }
  double eta(Eigen::Vector2d x) const override { return 1.0; }
  bool EssentialConditionsOnEdge(const lf::mesh::Entity& edge) const override {
    return edge_sel_dir_(edge);
  }
  bool IsImpedanceEdge(const lf::mesh::Entity& edge) const override {
    return edge_sel_imp_(edge);
  }
  double f(Eigen::Vector2d /*x*/) const override { return -2.0; }
  double g(Eigen::Vector2d x) const override {
    return 0.5 * std::log(x.squaredNorm() + 1);
  }
  double h(Eigen::Vector2d x) const override {
    if (x[1] < 1E-7) {
      return (0.5 * std::log(x.squaredNorm() + 1.0));
    } else if (x[0] > x[1]) {
      return x[0];
    } else {
      return x[1];
    }
  }

 private:
  /// Selector for Dirichlet boundary conditions
  lf::refinement::EntityCenterPositionSelector<
      std::function<bool(Eigen::Vector2d)>>
      edge_sel_dir_;
  /// Selector for Dirichlet boundary conditions
  lf::refinement::EntityCenterPositionSelector<
      std::function<bool(Eigen::Vector2d)>>
      edge_sel_imp_;
};  // namespace lf::uscalfe::test

TEST(lfe_bvpfe, bvp_DirNeu) {
  // Set debugging flags
  // LFELinSys_ctrl = 0;
  // The solution and its gradient
  auto u = MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return (0.5 * std::log(x.squaredNorm() + 1.0));
  });
  auto grad_u = MeshFunctionGlobal([](Eigen::Vector2d x) -> Eigen::Vector2d {
    return (1.0 / (x.squaredNorm() + 1.0)) * x;
  });
  auto bvp_p = std::make_shared<MixedEllipticBVP>();
  auto errnorms = TestConvergenceEllBVPFESol(REFLEV, bvp_p, u, grad_u);

  std::array<double, 8> exp_l2err{0.025494,    0.00736786,  0.0019407,
                                  0.000493379, 0.000123978, 3.10413e-05,
                                  7.76373e-06, 1.94118e-06};
  std::array<double, 8> exp_h1serr{0.104398,   0.0565337,  0.0290143,
                                   0.014624,   0.00732936, 0.00366718,
                                   0.00183394, 0.00091702};

  int k = 0;
  for (auto& err : errnorms) {
    if (k < 8) {
      EXPECT_NEAR(err.first, exp_l2err[k], 1.0E-5)
          << "Dubious L2 norm on level " << k;
      EXPECT_NEAR(err.second, exp_h1serr[k], 1.0E-5)
          << "Dubious H1 norm on level " << k;
      k++;
    }
  }
}

}  // namespace lf::uscalfe::test
