/**
 * @file
 * @brief Check that the Eigen tools work as expected.
 * @author Raffael Casagrande
 * @date   2019-01-19 07:27:25
 * @copyright MIT License
 */

// clang-format off
#include <Eigen/Eigen>
// clang-format on

#include <gtest/gtest.h>
#include <lf/base/base.h>

#include "lf/base/eigen_tools.h"

namespace lf::base::test {

TEST(eigenTools, IsEigen) {
  EXPECT_FALSE(is_eigen_matrix<double>);
  EXPECT_FALSE(is_eigen_matrix<int>);
  EXPECT_FALSE(is_eigen_matrix<std::vector<Eigen::Vector2d>>);
  EXPECT_TRUE(is_eigen_matrix<Eigen::VectorXd>);
  EXPECT_TRUE(is_eigen_matrix<Eigen::Vector2d>);
  EXPECT_TRUE(is_eigen_matrix<Eigen::Matrix3cd>);
}

}  // namespace lf::base::test
