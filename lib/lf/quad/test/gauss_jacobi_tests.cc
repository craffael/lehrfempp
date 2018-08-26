/**
 * @file
 * @brief Test the functioning of the method GaussJacobi
 * @author Raffael Casagrande
 * @date   2018-08-26 09:34:12
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/quad/quad.h>
#include "lf/quad/gauss_quadrature.h"

namespace lf::quad::test {

double IntegrateMonomial(quadOrder_t num_points, quadOrder_t exponent,
                         double alpha, double beta) {
  auto [points, weights] = GaussJacobi(num_points, alpha, beta);

  EXPECT_EQ(points.size(), num_points);
  EXPECT_EQ(weights.size(), num_points);

  double result = 0;
  for (quadOrder_t i = 0; i < num_points; ++i) {
    result += std::pow(points(i), exponent) * weights(i);
  }
  return result;
}

TEST(quad_GaussJacobi, Legendre) {
  EXPECT_DOUBLE_EQ(IntegrateMonomial(1, 0, 0, 0), 2);
  EXPECT_DOUBLE_EQ(IntegrateMonomial(1, 1, 0, 0), 0);
  EXPECT_GT(std::abs(IntegrateMonomial(1, 2, 0, 0) - 2. / 3.), 1e-12);
  EXPECT_NEAR(IntegrateMonomial(2, 2, 0, 0), 2. / 3., 1e-12);
  EXPECT_DOUBLE_EQ(IntegrateMonomial(2, 3, 0, 0), 0);
  EXPECT_GT(std::abs(IntegrateMonomial(2, 4, 0, 0) - 2. / 3.), 1e-12);
  EXPECT_NEAR(IntegrateMonomial(3, 4, 0, 0), 2. / 5., 1e-12);
  EXPECT_DOUBLE_EQ(IntegrateMonomial(3, 5, 0, 0), 0);
}

TEST(quad_GaussJacobi, AlphaOneBetaZero) {
  EXPECT_DOUBLE_EQ(IntegrateMonomial(1, 0, 1, 0), 2);
  EXPECT_DOUBLE_EQ(IntegrateMonomial(1, 1, 1, 0), -2. / 3.);
  EXPECT_GT(std::abs(IntegrateMonomial(1, 2, 1, 0) - 2. / 3.), 1e-12);
  EXPECT_NEAR(IntegrateMonomial(2, 2, 1, 0), 2. / 3., 1e-12);
  EXPECT_NEAR(IntegrateMonomial(2, 3, 1, 0), -2. / 5., 1e-12);
  EXPECT_GT(std::abs(IntegrateMonomial(2, 4, 1, 0) - 2. / 5.), 1e-12);
  EXPECT_NEAR(IntegrateMonomial(3, 4, 1, 0), 2. / 5., 1e-12);
  EXPECT_NEAR(IntegrateMonomial(3, 5, 1, 0), -2. / 7., 1e-12);
}

TEST(quad_GaussJacobi, AlphaZeroBetaOne) {
  EXPECT_DOUBLE_EQ(IntegrateMonomial(1, 0, 0, 1), 2);
  EXPECT_DOUBLE_EQ(IntegrateMonomial(1, 1, 0, 1), 2. / 3.);
  EXPECT_GT(std::abs(IntegrateMonomial(1, 2, 0, 1) - 2. / 3.), 1e-12);
  EXPECT_NEAR(IntegrateMonomial(2, 2, 0, 1), 2. / 3., 1e-12);
  EXPECT_NEAR(IntegrateMonomial(2, 3, 0, 1), 2. / 5., 1e-12);
  EXPECT_GT(std::abs(IntegrateMonomial(2, 4, 0, 1) - 2. / 5.), 1e-12);
  EXPECT_NEAR(IntegrateMonomial(3, 4, 0, 1), 2. / 5., 1e-12);
  EXPECT_NEAR(IntegrateMonomial(3, 5, 0, 1), 2. / 7., 1e-12);
}

TEST(quad_GaussJacobi, AlphaOneBetaOne) {
  EXPECT_DOUBLE_EQ(IntegrateMonomial(1, 0, 1, 1), 4. / 3.);
  EXPECT_DOUBLE_EQ(IntegrateMonomial(1, 1, 1, 1), 0.);
  EXPECT_GT(std::abs(IntegrateMonomial(1, 2, 1, 1) - 4. / 15.), 1e-12);
  EXPECT_NEAR(IntegrateMonomial(2, 2, 1, 1), 4. / 15., 1e-12);
  EXPECT_NEAR(IntegrateMonomial(2, 3, 1, 1), 0. / 5., 1e-12);
  EXPECT_GT(std::abs(IntegrateMonomial(2, 4, 1, 1) - 4. / 35.), 1e-12);
  EXPECT_NEAR(IntegrateMonomial(3, 4, 1, 1), 4. / 35., 1e-12);
  EXPECT_NEAR(IntegrateMonomial(3, 5, 1, 1), 0., 1e-12);
  GaussJacobi(30, 1, 0);
}

}  // namespace lf::quad::test
