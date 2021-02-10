/**
 * @file
 * @brief
 * @author Raffael Casagrande
 * @date   2021-01-29 10:33:15
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/brep/test_utils/check_brep_geometry.h>
#include <lf/brep/test_utils/curve_circle.h>

namespace lf::brep::test_utils::test {

TEST(lf_brep_testUtils, CurveCircleTest) {
  Eigen::Vector3d origin(2, 3, 0);
  double radius = 5;
  CurveCircle cc(origin, radius);

  auto local = Eigen::RowVectorXd::LinSpaced(10, -base::kPi, base::kPi);

  CheckBrepGeometry(cc, local);
}
}  // namespace lf::brep::test_utils::test
