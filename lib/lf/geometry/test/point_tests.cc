#include <gtest/gtest.h>
#include <lf/geometry/geometry.h>

namespace lf::geometry::test {

TEST(PointTest, checkCoord) {
  Point p(Eigen::Vector3d(1, 0, 0));
  EXPECT_EQ(p.DimGlobal(), 3);
  EXPECT_EQ(p.DimLocal(), 0);
  EXPECT_EQ(p.Global({0, 1}), Eigen::Vector3d(1, 0, 0));
  EXPECT_EQ(p.RefEl(), base::RefEl::kPoint());
  EXPECT_EQ(p.subGeometry(0, 0)->Global({0, 1}), Eigen::Vector3d(1, 0, 0));
}

}  // namespace lf::geometry::test
