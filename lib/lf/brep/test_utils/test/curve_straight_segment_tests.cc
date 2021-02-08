/**
 * @file
 * @brief Test implementation of CurveStraightSegment
 * @author Raffael Casagrande
 * @date   2021-02-08 04:13:15
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/brep/test_utils/check_brep_geometry.h>
#include <lf/brep/test_utils/curve_straight_segment.h>

namespace lf::brep::test_utils::test {

TEST(lf_brep_testUtils, CurveStraightSegmentTest) {
  Eigen::Vector3d start(-1, 2, 3);
  Eigen::Vector3d end(2, 1, 0);
  CurveStraightLine line(start, end);

  auto local = Eigen::RowVectorXd::LinSpaced(10, 0, 1.).eval();
  CheckBrepGeometry(line, local);
}

}  // namespace lf::brep::test_utils::test
