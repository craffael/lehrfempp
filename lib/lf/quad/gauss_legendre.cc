/**
 * @file
 * @brief Implementation of the GaussLegendre function to compute
 *        quadrature Weights and nodes
 * @author Raffael Casagrande
 * @date   2018-08-11 09:08:26
 * @copyright MIT License
 */

#include "gauss_legendre.h"
#include <boost/multiprecision/cpp_bin_float.hpp>

namespace lf::quad {
std::tuple<Eigen::VectorXd, Eigen::VectorXd> GaussLegendre(
    unsigned num_points) {
  LF_ASSERT_MSG(num_points > 0, "num_points must be positive.");

  using scalar_t =
      boost::multiprecision::number<boost::multiprecision::cpp_bin_float<
          57, boost::multiprecision::digit_base_2>>;

  static const scalar_t kPi = boost::math::constants::pi<scalar_t>();

  Eigen::VectorXd points(num_points), weights(num_points);

  // the roots are symmetric in the interval, so we only have to find half of
  // them
  unsigned int m = (num_points + 1) / 2;

  // approximation for the roots:
  for (int i = 0; i < m; ++i) {
    // initial guess for the i-th root:
    scalar_t z = cos(kPi * (i + 0.75) / (num_points + 0.5));
    scalar_t z1, pp;

    // start of newton
    do {
      // calculate value of legendre polynomial at z using recurrence relation:
      scalar_t p1 = 1.0;
      scalar_t p2 = 0.0;
      scalar_t p3;

      for (int j = 0; j < num_points; ++j) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j + 1.) * z * p2 - j * p3) / (j + 1.);
      }

      // p1 is now the value of the legendre polynomial at z. Next we compute
      // its derivative using a standard relation involving also p2, the
      // polynomial of one lower order:
      pp = num_points * (z * p1 - p2) / (z * z - 1.0);

      z1 = z;
      z = z1 - p1 / pp;
    } while (abs(z - z1) > 1e-17);

    points(i) = (0.5 * (1 - z)).convert_to<double>();
    points(num_points - 1 - i) = (0.5 * (1 + z)).convert_to<double>();
    weights(i) = (1. / ((1.0 - z * z) * pp * pp)).convert_to<double>();
    weights(num_points - 1 - i) = weights(i);
  }

  return {points, weights};
}
}  // namespace lf::quad
