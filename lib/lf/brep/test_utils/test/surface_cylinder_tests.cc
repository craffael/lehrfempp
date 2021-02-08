
/**
 * @file
 * @brief Test implementation of SurfaceCylinder
 * @author Raffael Casagrande
 * @date   2021-02-03 09:13:44
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/brep/test_utils/check_brep_geometry.h>
#include <lf/brep/test_utils/surface_cylinder.h>

#include <unsupported/Eigen/KroneckerProduct>

namespace lf::brep::test_utils::test {
TEST(lf_brep_testUtils, SurfaceCylinderTest) {
  Eigen::Vector3d origin(1, 1, 1);
  double radius = 2;
  Eigen::Vector3d axis(0, 3, 0);

  SurfaceCylinder cylinder(origin, axis, radius);
  Eigen::MatrixXd local(2, 100);
  local.row(0) = Eigen::KroneckerProduct(
      Eigen::RowVectorXd::LinSpaced(10, -base::kPi, base::kPi),
      Eigen::RowVectorXd::Constant(10, 1.));
  local.row(1) =
      Eigen::KroneckerProduct(Eigen::RowVectorXd::Constant(10, 1),
                              Eigen::RowVectorXd::LinSpaced(10, 0, 1.));
  CheckBrepGeometry(cylinder, local);

  // check projection:
  auto [dist, p] = cylinder.Project({1, 3, 3});
  EXPECT_NEAR(dist, 0, 1e-6);
  std::tie(dist, p) = cylinder.Project({1, 3, 4});
  EXPECT_NEAR(dist, 1, 1e-6);
}
}  // namespace lf::brep::test_utils::test
