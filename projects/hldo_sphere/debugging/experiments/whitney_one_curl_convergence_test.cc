/**
 * @file whitney_one_curl_convergence_test.cc
 */

#include <whitney_one_curl_test.h>

/**
 * @brief Prints the L2 norm
 */
int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " max_refinement_level" << std::endl;
    exit(1);
  }

  unsigned refinement_level = atoi(argv[1]);
  std::cout << "max_refinement_level : " << refinement_level << std::endl;

  // mathematica function output requries the following helpers
  auto Power = [](double a, double b) -> double { return std::pow(a, b); };
  auto Sin = [](double a) -> double { return std::sin(a); };
  auto Cos = [](double a) -> double { return std::cos(a); };
  auto Sqrt = [](double a) -> double { return std::sqrt(a); };

  // Compute the analytic solution of the problem
  auto u_one = [&](const Eigen::Vector3d x_vec) -> Eigen::Vector3d {
    // first scale to the circle
    Eigen::Vector3d x_ = x_vec / x_vec.norm();
    double x = x_(0);
    double y = x_(1);
    double z = x_(2);

    Eigen::VectorXd ret(3);

    // mathematica autocompute
    ret << -y, x, 0;
    return ret;
  };

  // Compute the analytic solution of the problem
  auto u_two = [&](const Eigen::Vector3d x_vec) -> double {
    // first scale to the circle
    Eigen::Vector3d x_ = x_vec / x_vec.norm();
    double x = x_(0);
    double y = x_(1);
    double z = x_(2);

    return sin(z);
  };

  double ana_sol = -7.5691944739878640736;

  projects::hldo_sphere::debugging::WhitneyOneCurlTest test;
  test.SetAnaSol(ana_sol);
  test.SetFunctions(u_two, u_one);
  test.Compute(refinement_level);

  return 0;
}
