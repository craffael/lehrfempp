/**
 * @file
 * @brief Implementation of Jacobian() test for geometry objects
 * @author Anian Ruoss
 * @date   2019-02-11 17:30:17
 * @copyright MIT License
 */

#include "check_jacobian.h"

#include <gtest/gtest.h>

namespace lf::geometry::test_utils {

void checkJacobian(const lf::geometry::Geometry &geom,
                   const Eigen::MatrixXd &eval_points,
                   const double &tolerance) {
  const double h = 1e-6;

  const size_t num_points = eval_points.cols();
  const size_t dim_local = geom.DimLocal();
  const size_t dim_global = geom.DimGlobal();

  Eigen::MatrixXd jacobians = geom.Jacobian(eval_points);

  EXPECT_EQ(jacobians.rows(), dim_global) << "Jacobian has " << jacobians.rows()
                                          << " rows instead of " << dim_global;
  EXPECT_EQ(jacobians.cols(), num_points * dim_local)
      << "Jacobian has " << jacobians.cols() << " cols instead of "
      << num_points * dim_local;

  for (size_t j = 0; j < num_points; ++j) {
    auto point = eval_points.col(j);

    Eigen::MatrixXd jacobian =
        jacobians.block(0, j * dim_local, dim_global, dim_local);
    Eigen::MatrixXd approx_jacobian =
        Eigen::MatrixXd::Zero(dim_global, dim_local);

    for (size_t i = 0; i < dim_local; ++i) {
      Eigen::VectorXd h_vec = Eigen::VectorXd::Zero(dim_local);
      h_vec(i) = h;

      // approximate gradient with symmetric difference quotient
      approx_jacobian.col(i) =
          (geom.Global(point + h_vec) - geom.Global(point - h_vec)) / (2. * h);
    }

    EXPECT_TRUE(jacobian.isApprox(approx_jacobian, tolerance))
        << "Jacobian incorrect at point " << point;
  }
}

}  // namespace lf::geometry::test_utils
