/**
 * @file
 * @brief Simple tests to check if the GaussLegendre method returns reasonable
 *        values.
 * @author Raffael Casagrande
 * @date   2018-08-11 09:48:06
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/quad/gauss_legendre.h>

namespace lf::quad::test {
TEST(GaussLegendre, compareToTable) {
  // n = 1
  auto [points, weights] = GaussLegendre(1);
  EXPECT_EQ(points.size(), 1);
  EXPECT_EQ(points(0), 0.5);
  EXPECT_EQ(weights.size(), 1);
  EXPECT_EQ(weights(0), 1.);

  // n = 2
  std::tie(points, weights) = GaussLegendre(2);
  EXPECT_EQ(points.size(), 2);
  EXPECT_DOUBLE_EQ(points(0), (1 - std::sqrt(1. / 3.)) / 2.);
  EXPECT_DOUBLE_EQ(points(1), (1 + std::sqrt(1. / 3.)) / 2.);
  EXPECT_EQ(weights.size(), 2);
  EXPECT_DOUBLE_EQ(weights(0), 0.5);
  EXPECT_DOUBLE_EQ(weights(1), 0.5);

  // n = 3
  std::tie(points, weights) = GaussLegendre(3);
  EXPECT_EQ(points.size(), 3);
  EXPECT_DOUBLE_EQ(points(0), (1 - std::sqrt(3. / 5.)) / 2.);
  EXPECT_EQ(points(1), 0.5);
  EXPECT_DOUBLE_EQ(points(2), (1 + std::sqrt(3. / 5.)) / 2.);
  EXPECT_EQ(weights.size(), 3);
  EXPECT_DOUBLE_EQ(weights(0), 5. / 18.);
  EXPECT_DOUBLE_EQ(weights(1), 4. / 9.);
  EXPECT_DOUBLE_EQ(weights(2), 5. / 18.);

  // n = 4
  std::tie(points, weights) = GaussLegendre(4);
  EXPECT_EQ(points.size(), 4);
  EXPECT_DOUBLE_EQ(
      points(0), (1 - std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.))) / 2.);
  EXPECT_DOUBLE_EQ(
      points(1), (1 - std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.))) / 2.);
  EXPECT_DOUBLE_EQ(points(2), 1 - points(1));
  EXPECT_DOUBLE_EQ(points(3), 1 - points(0));
  EXPECT_EQ(weights.size(), 4);
  EXPECT_DOUBLE_EQ(weights(0), (18 - std::sqrt(30.)) / 72.);
  EXPECT_DOUBLE_EQ(weights(1), (18 + std::sqrt(30.)) / 72.);
  EXPECT_DOUBLE_EQ(weights(2), weights(1));
  EXPECT_DOUBLE_EQ(weights(3), weights(0));

  GaussLegendre(100);
}

}  // namespace lf::quad::test
