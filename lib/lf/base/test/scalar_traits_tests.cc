/**
 * @file
 * @brief Test functionality in scalar_traits.h
 * @author Raffael Casagrande
 * @date   2021-01-01 04:01:43
 * @copyright MIT License
 */

#include <gtest/gtest.h>

#include <lf/base/base.h>

namespace lf::base::test {

TEST(ScalarTraits, IsScalar) {
  EXPECT_TRUE(is_scalar<double>);
  EXPECT_TRUE(is_scalar<float>);
  EXPECT_TRUE(is_scalar<int>);
  EXPECT_TRUE(is_scalar<std::complex<double>>);
  EXPECT_TRUE(is_scalar<std::complex<int>>);

  EXPECT_FALSE(is_scalar<std::string>);
  EXPECT_FALSE(is_scalar<Eigen::MatrixXd>);
}

}  // namespace lf::base::test
