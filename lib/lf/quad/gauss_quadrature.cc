/**
 * @file
 * @brief Implementation of the GaussLegendre function to compute
 *        quadrature Weights and nodes
 * @author Raffael Casagrande
 * @date   2018-08-11 09:08:26
 * @copyright MIT License
 */

#include "gauss_quadrature.h"

#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>

namespace lf::quad {
std::tuple<Eigen::VectorXd, Eigen::VectorXd> GaussLegendre(
    unsigned num_points) {
  LF_ASSERT_MSG(num_points > 0, "num_points must be positive.");

  using scalar_t =
      boost::multiprecision::number<boost::multiprecision::cpp_bin_float<
          57, boost::multiprecision::digit_base_2>>;

  static const scalar_t kPi = boost::math::constants::pi<scalar_t>();

  Eigen::VectorXd points(num_points);
  Eigen::VectorXd weights(num_points);

  // the roots are symmetric in the interval, so we only have to find half of
  // them
  unsigned int m = (num_points + 1) / 2;

  // approximation for the roots:
  for (unsigned int i = 0; i < m; ++i) {
    // initial guess for the i-th root:
    scalar_t z = cos(kPi * (i + 0.75) / (num_points + 0.5));
    scalar_t z1;
    scalar_t pp;

    // start of newton
    do {
      // calculate value of legendre polynomial at z using recurrence relation:
      scalar_t p1 = 1.0;
      scalar_t p2 = 0.0;
      scalar_t p3;

      for (unsigned int j = 0; j < num_points; ++j) {
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

  return {std::move(points), std::move(weights)};
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd> GaussJacobi(
    quadDegree_t num_points, double alpha, double beta) {
  LF_ASSERT_MSG(num_points > 0, "num_points must be positive.");
  LF_ASSERT_MSG(alpha > -1, "alpha > -1 required");
  LF_ASSERT_MSG(beta > -1, "beta > -1 required.");

  using boost::math::lgamma;
  using scalar_t =
      boost::multiprecision::number<boost::multiprecision::cpp_bin_float<
          57, boost::multiprecision::digit_base_2>>;

  int MAX_IT = 10;

  scalar_t alfbet;
  scalar_t an;
  scalar_t bn;
  scalar_t r1;
  scalar_t r2;
  scalar_t r3;
  scalar_t a;
  scalar_t b;
  scalar_t c;
  scalar_t p1;
  scalar_t p2;
  scalar_t p3;
  scalar_t pp;
  scalar_t temp;
  scalar_t z;
  scalar_t z1;

  std::vector<scalar_t> points(num_points);
  Eigen::VectorXd weights(num_points);

  // Make an initial guess for the zeros:
  for (quadDegree_t i = 0; i < num_points; ++i) {
    if (i == 0) {
      // initial guess for the largest root
      an = alpha / num_points;
      bn = beta / num_points;
      r1 = (1.0 + alpha) *
           (2.78 / (4.0 + num_points * num_points) + 0.768 * an / num_points);
      r2 = 1.0 + 1.48 * an + 0.96 * bn + 0.452 * an * an + 0.83 * an * bn;
      z = 1.0 - r1 / r2;
    } else if (i == 1) {
      // initial guess for the second largest root
      r1 = (4.1 + alpha) / ((1.0 + alpha) * (1.0 + 0.156 * alpha));
      r2 = 1.0 + 0.06 * (num_points - 8.0) * (1.0 + 0.12 * alpha) / num_points;
      r3 = 1.0 + 0.012 * beta * (1.0 + 0.25 * std::abs(alpha)) / num_points;
      z -= (1.0 - z) * r1 * r2 * r3;
    } else if (i == 2) {
      // initial guess for the third largest root
      r1 = (1.67 + 0.28 * alpha) / (1.0 + 0.37 * alpha);
      r2 = 1.0 + 0.22 * (num_points - 8.0) / num_points;
      r3 = 1.0 + 8.0 * beta / ((6.28 + beta) * num_points * num_points);
      z -= (points[0] - z) * r1 * r2 * r3;
    } else if (i == num_points - 2) {
      // initial guess for the second smallest root
      r1 = (1.0 + 0.235 * beta) / (0.766 + 0.119 * beta);
      r2 = 1.0 / (1.0 + 0.639 * (num_points - 4.0) /
                            (1.0 + 0.71 * (num_points - 4.0)));
      r3 = 1.0 /
           (1.0 + 20.0 * alpha / ((7.5 + alpha) * num_points * num_points));
      z += (z - points[num_points - 4]) * r1 * r2 * r3;
    } else if (i == num_points - 1) {
      // initial guess for the smallest root
      r1 = (1.0 + 0.37 * beta) / (1.67 + 0.28 * beta);
      r2 = 1.0 / (1.0 + 0.22 * (num_points - 8.0) / num_points);
      r3 = 1.0 /
           (1.0 + 8.0 * alpha / ((6.28 + alpha) * num_points * num_points));
      z += (z - points[num_points - 3]) * r1 * r2 * r3;
    } else {
      // initial guess for the other points
      z = 3.0 * points[i - 1] - 3.0 * points[i - 2] + points[i - 3];
    }
    alfbet = alpha + beta;
    quadDegree_t its;
    for (its = 1; its <= MAX_IT; ++its) {
      // refinement by Newton's method
      temp = 2.0 + alfbet;

      // Start the recurrence with P_0 and P1 to avoid a division by zero when
      // alpha * beta = 0 or -1
      p1 = (alpha - beta + temp * z) / 2.0;
      p2 = 1.0;
      for (quadDegree_t j = 2; j <= num_points; ++j) {
        p3 = p2;
        p2 = p1;
        temp = 2 * j + alfbet;
        a = 2 * j * (j + alfbet) * (temp - 2.0);
        b = (temp - 1.0) *
            (alpha * alpha - beta * beta + temp * (temp - 2.0) * z);
        c = 2.0 * (j - 1 + alpha) * (j - 1 + beta) * temp;
        p1 = (b * p2 - c * p3) / a;
      }
      pp = (num_points * (alpha - beta - temp * z) * p1 +
            2.0 * (num_points + alpha) * (num_points + beta) * p2) /
           (temp * (1.0 - z * z));
      // p1 is now the desired jacobia polynomial. We next compute pp, its
      // derivative, by a standard relation involving p2, the polynomial of one
      // lower order
      z1 = z;
      z = z1 - p1 / pp;  // Newtons Formula
      if (abs(z - z1) <= 1e-17) {
        break;
      }
    }
    LF_VERIFY_MSG(its <= MAX_IT, "too many iterations.");

    points[i] = z;
    weights(i) = (exp(lgamma(alpha + num_points) + lgamma(beta + num_points) -
                      lgamma(num_points + 1.) -
                      lgamma(double(num_points) + alfbet + 1.0)) *
                  temp * pow(2.0, alfbet) / (pp * p2))
                     .convert_to<double>();
  }

  Eigen::VectorXd points_result(num_points);
  for (quadDegree_t i = 0; i < num_points; ++i) {
    points_result(i) = points[i].convert_to<double>();
  }

  return {points_result, weights};
}
}  // namespace lf::quad
