/**
 * @file
 * @brief Test the OcctBrepModel class
 * @author Raffael Casagrande
 * @date   2020-11-11 02:56:08
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include "utils.h"

namespace lf::brep::occt::test {

TEST(occt, fileDoesntExist) {
  EXPECT_DEATH(OcctBrepModel("invalid.brep2"), "Could not open file");
}

TEST(occt, cubeEdgeTest) {
  auto model = LoadModel("cube_rotated.brep");
  EXPECT_EQ(model->NumGeometries(1), 12);
  EXPECT_EQ(model->NumGeometries(2), 6);

  // retrieve edge from [20,0,0] -> [20,10,0]:
  Eigen::MatrixXd global(3, 2);
  // clang-format off
  global << 20,  20,
            0,  10,
            0,  0;
  // clang-format on
  auto geom = model->FindGeometry(1, global);
  EXPECT_TRUE(geom.first);
  EXPECT_EQ(geom.second.cols(), 2);
  EXPECT_TRUE(geom.first->Global(geom.second).isApprox(global));
}

}  // namespace lf::brep::occt::test
