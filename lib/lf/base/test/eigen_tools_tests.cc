/**
 * @file
 * @brief Check that the Eigen tools work as expected.
 * @author Raffael Casagrande
 * @date   2019-01-19 07:27:25
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/base/base.h>

#include <Eigen/Core>

#include "lf/base/eigen_tools.h"

namespace lf::base::test {

TEST(eigenTools, IsEigen) {
  EXPECT_FALSE(EigenMatrix<double>);
  EXPECT_FALSE(EigenMatrix<int>);
  EXPECT_FALSE(EigenMatrix<std::vector<Eigen::Vector2d>>);
  EXPECT_TRUE(EigenMatrix<Eigen::VectorXd>);
  EXPECT_TRUE(EigenMatrix<Eigen::Vector2d>);
  EXPECT_TRUE(EigenMatrix<Eigen::Matrix3cd>);
  EXPECT_FALSE(EigenMatrix<decltype(Eigen::Vector3d(1,2,3) * 5)>);
  EXPECT_FALSE(EigenMatrix<decltype(Eigen::Vector3d::Zero())>);
}

}  // namespace lf::base::test
