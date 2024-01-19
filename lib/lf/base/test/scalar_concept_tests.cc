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
  EXPECT_TRUE(Scalar<double>);
  EXPECT_TRUE(Scalar<float>);
  EXPECT_TRUE(Scalar<int>);
  EXPECT_TRUE(Scalar<std::complex<double>>);
  EXPECT_TRUE(Scalar<std::complex<int>>);

  EXPECT_FALSE(Scalar<std::string>);
  EXPECT_FALSE(Scalar<Eigen::MatrixXd>);
}

}  // namespace lf::base::test
