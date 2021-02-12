/**
 * @file
 * @brief Test BrepSurfaceSegment class
 * @author Raffael Casagrande
 * @date   2021-02-12 02:05:06
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/brep/geom/geom.h>
#include <lf/brep/test_utils/surface_cylinder.h>

#include "lf/brep/test_utils/check_brep_geometry.h"

namespace lf::brep::geom::test {

TEST(lf_brep_geom, BrepSurfaceSegment) {
  auto unit_cylinder = std::make_shared<test_utils::SurfaceCylinder>(
      Eigen::Vector3d(0, 0, 0), Eigen::Vector3d::UnitZ(), 1.);

  // construct a line that makes a circle:
  auto half_circle = BrepSurfaceSegment(
      unit_cylinder, (Eigen::Matrix2d() << 0, base::kPi, 0, 0).finished());
  auto period = half_circle.Periods();
  ASSERT_NEAR(period(0), 2, 1e-6);
  test_utils::CheckBrepGeometry(half_circle,
                                Eigen::RowVectorXd::LinSpaced(11, 0, 1));

  // construct a line that goes along the z-axis:
  auto z_axis = BrepSurfaceSegment(
      unit_cylinder, (Eigen::Matrix2d() << 0, 0, 0, 1).finished());
  period = z_axis.Periods();
  ASSERT_EQ(period(0), 0);
  test_utils::CheckBrepGeometry(z_axis,
                                Eigen::RowVectorXd::LinSpaced(11, 0, 1));

  // construct a spiral:
  auto spiral = BrepSurfaceSegment(
      unit_cylinder, (Eigen::Matrix2d() << 0, 3 * base::kPi, 0, 1).finished());
  period = spiral.Periods();
  ASSERT_EQ(period(0), 0);
  test_utils::CheckBrepGeometry(spiral,
                                Eigen::RowVectorXd::LinSpaced(11, 0, 1));
}

}  // namespace lf::brep::geom::test
