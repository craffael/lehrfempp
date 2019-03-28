/**
 * @file
 * @brief Implementation of IntegrationElement() test for geometry objects
 * @author Anian Ruoss
 * @date   2019-02-11 17:59:17
 * @copyright MIT License
 */

#include "check_integration_element.h"
#include <gtest/gtest.h>

namespace lf::geometry::test_utils {

void checkIntegrationElement(const lf::geometry::Geometry &geom,
                             const Eigen::MatrixXd &eval_points) {
  const size_t num_points = eval_points.cols();
  const size_t dim_local = geom.DimLocal();
  const size_t dim_global = geom.DimGlobal();

  Eigen::MatrixXd jacobians = geom.Jacobian(eval_points);
  Eigen::VectorXd integrationElements = geom.IntegrationElement(eval_points);

  EXPECT_EQ(integrationElements.rows(), num_points)
      << "IntegrationElement has " << integrationElements.rows()
      << " rows instead of " << num_points;
  EXPECT_EQ(integrationElements.cols(), 1)
      << "IntegrationElement has " << integrationElements.cols()
      << " cols instead of " << 1;

  for (int j = 0; j < num_points; ++j) {
    Eigen::MatrixXd jacobian =
        jacobians.block(0, j * dim_local, dim_global, dim_local);

    const double integrationElement = integrationElements(j);
    const double approx_integrationElement =
        std::sqrt((jacobian.transpose() * jacobian).determinant());

    EXPECT_FLOAT_EQ(integrationElement, approx_integrationElement)
        << "IntegrationElement incorrect at point " << eval_points.col(j);
  }
}

}  // namespace lf::geometry::test_utils
