/**
 * @file
 * @brief Implementation of JacobianInverseGramian() test for geometry objects
 * @author Anian Ruoss
 * @date   2019-02-11 17:53:17
 * @copyright MIT License
 */

#include "check_jacobian_inverse_gramian.h"

#include <gtest/gtest.h>

namespace lf::geometry::test_utils {

void checkJacobianInverseGramian(const lf::geometry::Geometry &geom,
                                 const Eigen::MatrixXd &eval_points) {
  const size_t num_points = eval_points.cols();
  const size_t dim_local = geom.DimLocal();
  const size_t dim_global = geom.DimGlobal();

  Eigen::MatrixXd jacobians = geom.Jacobian(eval_points);
  Eigen::MatrixXd jacInvGrams = geom.JacobianInverseGramian(eval_points);

  EXPECT_EQ(jacInvGrams.rows(), dim_global)
      << "JacobianInverseGramian has " << jacInvGrams.rows()
      << " rows instead of " << dim_global;
  EXPECT_EQ(jacInvGrams.cols(), num_points * dim_local)
      << "JacobianInverseGramian has " << jacInvGrams.cols()
      << " cols instead of " << num_points * dim_local;

  for (int j = 0; j < num_points; ++j) {
    Eigen::MatrixXd jacInvGram =
        jacInvGrams.block(0, j * dim_local, dim_global, dim_local);
    Eigen::MatrixXd jacobian =
        jacobians.block(0, j * dim_local, dim_global, dim_local);

    EXPECT_TRUE(jacInvGram.isApprox(
        jacobian * (jacobian.transpose() * jacobian).inverse()))
        << "JacobianInverseGramian incorrect at point " << eval_points.col(j);
  }
}

}  // namespace lf::geometry::test_utils
