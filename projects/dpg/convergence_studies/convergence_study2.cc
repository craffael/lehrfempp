#include <cmath>
#include <fstream>
#include <iomanip>

#include "primal_dpg.h"
#include "ultraweak_dpg.h"

lf::base::size_type reflev = 5;

int main() {
  std::cout << "Convergence study 2: Ultraweak DPG method, smooth solution \n";
  std::cout << "Exact solution: u(x,y) = sin(pi*x)sin(pi*y), homogenous "
               "dirichlet b.c. \n";
  std::cout << "BVP 1: epsilon_1=1.0, beta_1=[0.0,0.0] \n";
  std::cout << "BVP 2: epsilon_2=1.0, beta_2=[2.0,1.0] \n ";
  std::cout << "p = {0,1}";

  // specify the two boundary value problems:
  // exact solution:
  const double PI = lf::base::kPi;
  auto u_exact = [PI](const Eigen::Vector2d& x) -> double {
    return std::sin(PI * x[0]) * std::sin(PI * x[1]);
  };
  auto u_grad = [PI](const Eigen::Vector2d& x) -> Eigen::Vector2d {
    return (Eigen::VectorXd(2)
                << PI * std::cos(PI * x[0]) * std::sin(PI * x[1]),
            PI * std::sin(PI * x[0]) * std::cos(PI * x[1]))
        .finished();
  };

  // BVP1:
  auto epsilon_1 = [](const Eigen::Vector2d& /*x*/) -> double { return 1.0; };
  auto beta_1 = [](const Eigen::Vector2d& /*x*/) -> Eigen::Vector2d {
    return (Eigen::VectorXd(2) << 0.0, 0.0).finished();
  };
  auto f_1 = [PI](const Eigen::Vector2d& x) -> double {
    return 2.0 * PI * PI * std::sin(PI * x[0]) * std::sin(PI * x[1]);
  };
  auto g_1 = [](const Eigen::Vector2d& /*x*/) -> double { return 0.0; };

  // BVP 2:
  auto epsilon_2 = [](const Eigen::Vector2d& /*x*/) -> double { return 1.0; };
  auto beta_2 = [PI](const Eigen::Vector2d& /*x*/) -> Eigen::Vector2d {
    return (Eigen::VectorXd(2) << 2.0, 1.0).finished();
  };
  auto f_2 = [PI](const Eigen::Vector2d& x) -> double {
    return 2.0 * PI * PI * std::sin(PI * x[0]) * std::sin(PI * x[1]) +
           PI * (2.0 * std::cos(PI * x[0]) * std::sin(PI * x[1]) +
                 std::sin(PI * x[0]) * std::cos(PI * x[1]));
  };
  auto g_2 = [](const Eigen::Vector2d& /*x*/) -> double { return 0.0; };

  // BVP 1 on a triangular mesh
  // run the convergence tests choosing different approxiamtion orders and
  std::cout << "Running convergence tests for BVP 1 on a triangular mesh \n";
  for (int deg_p = 0; deg_p <= 1; deg_p++) {
    const int delta_p = 2;
    auto errors = projects::dpg::test::
        TestConververgenceUltraWeakDPGConvectionDiffusionDirichletBVP(
            reflev, u_exact, u_grad, epsilon_1, beta_1, f_1, g_1, deg_p,
            delta_p, lf::base::RefEl::kTria());

    // output to a corresponding file
    const std::string file_name = "cs2/triangular/bvp1_" + std::to_string(deg_p);
    std::ofstream file(file_name);
    for (auto [dofs, erroru, errorsigma, errorestimate] : errors) {
      file << dofs << ", " << erroru << "," << errorsigma << ","
           << errorestimate << '\n';
    }
    file.close();
  }

  // BVP 2 on a triangular mesh
  // run the convergence tests choosing different approxiamtion orders
  std::cout << "Running convergence tests for BVP 2 on a triangular mesh \n";
  for (int deg_p = 0; deg_p <= 1; deg_p++) {
    const int delta_p = 2;
    auto errors = projects::dpg::test::
        TestConververgenceUltraWeakDPGConvectionDiffusionDirichletBVP(
            reflev, u_exact, u_grad, epsilon_2, beta_2, f_2, g_2, deg_p,
            delta_p, lf::base::RefEl::kTria());

    // output to a corresponding file
    const std::string file_name = "cs2/triangular/bvp2_" + std::to_string(deg_p);
    std::ofstream file(file_name);
    for (auto [dofs, erroru, errorsigma, errorestimate] : errors) {
      file << dofs << ", " << erroru << "," << errorsigma << ","
           << errorestimate << '\n';
    }
    file.close();
  }

  // BVP 1 on a quadrilateral mesh
  // run the convergence tests choosing different approxiamtion orders
  std::cout << "Running convergence tests for BVP 1 on a quadrilateral mesh \n";
  for (int deg_p = 0; deg_p <= 1; deg_p++) {
    const int delta_p = 2;
    auto errors = projects::dpg::test::
        TestConververgenceUltraWeakDPGConvectionDiffusionDirichletBVP(
            reflev, u_exact, u_grad, epsilon_1, beta_1, f_1, g_1, deg_p,
            delta_p, lf::base::RefEl::kQuad());

    // output to a corresponding file
    const std::string file_name = "cs2/quadrilateral/bvp1_" + std::to_string(deg_p);
    std::ofstream file(file_name);
    for (auto [dofs, erroru, errorsigma, errorestimate] : errors) {
      file << dofs << ", " << erroru << "," << errorsigma << ","
           << errorestimate << '\n';
    }
    file.close();
  }

  // BVP 2 on a quadrilateral mesh
  // run the convergence tests choosing different approxiamtion orders and
  // enrichement degrees:
  std::cout << "Running convergence tests for BVP 2 on a quadrilateral mesh \n";
  for (int deg_p = 0; deg_p <= 1; deg_p++) {
    const int delta_p = 2;
    auto errors = projects::dpg::test::
        TestConververgenceUltraWeakDPGConvectionDiffusionDirichletBVP(
            reflev, u_exact, u_grad, epsilon_2, beta_2, f_2, g_2, deg_p,
            delta_p, lf::base::RefEl::kQuad());

    // output to a corresponding file
    const std::string file_name = "cs2/quadrilateral/bvp2_" + std::to_string(deg_p);
    std::ofstream file(file_name);
    for (auto [dofs, erroru, errorsigma, errorestimate] : errors) {
      file << dofs << ", " << erroru << "," << errorsigma << ","
           << errorestimate << '\n';
    }
    file.close();
  }

  return 0;
}
