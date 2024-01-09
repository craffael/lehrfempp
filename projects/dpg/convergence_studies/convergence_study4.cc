#include <cmath>
#include <fstream>
#include <iomanip>

#include "primal_dpg.h"
#include "ultraweak_dpg.h"

lf::base::size_type reflev =  5;

int main() {
  std::cout << "Convergence study 4: Ultraweak DPG method, solution with a "
               "boundary layer  \n";
  std::cout << "Solution develops a boundary layer at top and right edge. \n";

  const std::vector<double> epsilon_values = {1, 0.5, 0.1, 0.05, 0.02, 0.01};
  const double b1 = 2.0;
  const double b2 = 1.0;

  for (const double eps : epsilon_values) {
    //------------------------
    // Specification of the BVP
    //------------------------
    // exact solution:
    auto u1_exact = [eps, b1](double x) -> double {
      return x + (std::exp(x * b1 / eps) - 1) / (1 - std::exp(b1 / eps));
    };
    auto u1_exact_grad = [eps, b1](double x) -> double {
      return 1 + b1 / eps * std::exp(x * b1 / eps) / (1 - std::exp(b1 / eps));
    };
    auto u1_exact_d2 = [eps, b1](double x) -> double {
      return b1 * b1 / (eps * eps) * std::exp(x * b1 / eps) /
             (1 - std::exp(b1 / eps));
    };

    auto u2_exact = [eps, b2](double y) -> double {
      return y + (std::exp(y * b2 / eps) - 1) / (1 - std::exp(b2 / eps));
    };
    auto u2_exact_grad = [eps, b2](double y) -> double {
      return 1 + b2 / eps * std::exp(y * b2 / eps) / (1 - std::exp(b2 / eps));
    };
    auto u2_exact_d2 = [eps, b2](double y) -> double {
      return b2 * b2 / (eps * eps) * std::exp(y * b2 / eps) /
             (1 - std::exp(b2 / eps));
    };

    auto u_exact = [&u1_exact, &u2_exact](const Eigen::Vector2d& x) -> double {
      return u1_exact(x[0]) * u2_exact(x[1]);
    };
    auto sigma = [eps, &u1_exact, &u2_exact, &u1_exact_grad,
                  &u2_exact_grad](const Eigen::Vector2d& x) -> Eigen::Vector2d {
      return (Eigen::VectorXd(2) << eps * u1_exact_grad(x[0]) * u2_exact(x[1]),
              eps * u1_exact(x[0]) * u2_exact_grad(x[1]))
          .finished();
    };

    auto epsilon = [eps](const Eigen::Vector2d& /*x*/) -> double {
      return eps;
    };
    auto beta = [b1, b2](const Eigen::Vector2d& /*x*/) -> Eigen::Vector2d {
      return (Eigen::VectorXd(2) << b1, b2).finished();
    };

    auto f = [eps, b1, b2, &u1_exact, &u1_exact_grad, &u1_exact_d2, &u2_exact,
              &u2_exact_grad,
              &u2_exact_d2](const Eigen::Vector2d& x) -> double {
      return -eps * (u1_exact_d2(x[0]) * u2_exact(x[1]) +
                     u1_exact(x[0]) * u2_exact_d2(x[1])) +
             b1 * u1_exact_grad(x[0]) * u2_exact(x[1]) +
             b2 * u1_exact(x[0]) * u2_exact_grad(x[1]);
    };
    auto g = [](const Eigen::Vector2d& /*x*/) -> double { return 0.0; };

    // specification of degree and enrichement for the method
    const int deg_p = 1;
    const int delta_p = 2;

    //------------------------
    // run the convergence test
    //------------------------
    std::cout << "Running convergence tests on a triangular mesh  (eps= " << eps
              << ")\n";
    auto errors = projects::dpg::test::
        TestConververgenceUltraWeakDPGConvectionDiffusionDirichletBVP(
            reflev, u_exact, sigma, epsilon, beta, f, g, deg_p, delta_p,
            lf::base::RefEl::kTria());

    // output to a corresponding file
    const std::string file_name = "cs4/" + std::to_string(eps);
    std::ofstream file(file_name);
    for (auto [dofs, error_u, error_s, errorestimator] : errors) {
      file << dofs << ", " << error_u << ", " << error_s << ", "
           << errorestimator << '\n';
    }
    file.close();
  }

  return 0;
}
