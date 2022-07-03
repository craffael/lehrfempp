/**
 * @file whitney_one_basis_expansion_test.cc
 */

#include <whitney_one_basis_expansion_coeffs.h>

/**
 * @brief Prints the L2 norm
 */
int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " max_refinement_level " << std::endl;
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
  auto u = [&](const Eigen::Vector3d x_vec) -> Eigen::Vector3d {
    // first scale to the circle
    Eigen::Vector3d x_ = x_vec / x_vec.norm();
    double x = x_(0);
    double y = x_(1);
    double z = x_(2);

    Eigen::VectorXd ret(3);

    // mathematica autocompute
    ret << (-(x * z * Sin(x)) + (Power(y, 2) + Power(z, 2)) * Sin(y) -
            x * y * Sin(z)) /
               (Power(x, 2) + Power(y, 2) + Power(z, 2)),
        (-(y * z * Sin(x)) - x * y * Sin(y) +
         (Power(x, 2) + Power(z, 2)) * Sin(z)) /
            (Power(x, 2) + Power(y, 2) + Power(z, 2)),
        ((Power(x, 2) + Power(y, 2)) * Sin(x) - z * (x * Sin(y) + y * Sin(z))) /
            (Power(x, 2) + Power(y, 2) + Power(z, 2));
    return ret;
  };

  projects::hldo_sphere::debugging::WhitneyOneBasisExpansionCoeffs::Experiemnt(
      refinement_level, u);

  return 0;
}
