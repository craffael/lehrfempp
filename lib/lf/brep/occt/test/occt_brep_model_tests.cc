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

TEST(occt, rotatedCubeTest) {
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
  auto find_result = model->FindGeometries(1, global);
  EXPECT_EQ(find_result.size(), 1);
  EXPECT_TRUE(find_result[0].first);
  EXPECT_EQ(find_result[0].second.cols(), 2);
  EXPECT_EQ(find_result[0].second.rows(), 1);
  EXPECT_TRUE(
      find_result[0].first->Global(find_result[0].second).isApprox(global));

  // try to retrieve an edge that is outside of the parameter bounds:
  // clang-format off
  global << 20, 20,
             0, 20,
             0,  0;
  // clang-format on
  find_result = model->FindGeometries(1, global);
  EXPECT_EQ(find_result.size(), 0);

  //// retrieve face
  // global.resize(3, 4);
  //// clang-format off
  // global << 20, 20, 20, 20,
  //           0,  0, 10, 10,
  //          10,  0,  0, 10;
  //// clang-format on
  // find_result = model->FindGeometries(2, global);
  // EXPECT_EQ(find_result.size(), 1);
  // EXPECT_TRUE(find_result[0].first);
  // EXPECT_EQ(find_result[0].second.rows(), 2);
  // EXPECT_EQ(find_result[0].second.cols(), 4);
  // EXPECT_TRUE(
  //    find_result[0].first->Global(find_result[0].second).isApprox(global));
}

// Here we check that FindGeometries() gives the correct result also if the
// global point lies on the curve but the parameter is out of bounds.
TEST(occt, bSpline2dTest) {
  auto model = LoadModel("bspline_2d.brep");
  EXPECT_EQ(model->NumGeometries(1), 2);
  EXPECT_EQ(model->NumGeometries(2), 0);

  // get both two curves:
  Eigen::MatrixXd global(3, 2);
  // clang-format off
  global << -10, -13,
            -10, -6,
              0,  0;
  // clang-format on
  auto find_result = model->FindGeometries(1, global);
  EXPECT_EQ(find_result.size(), 2);

  // find the b-spline curve (which is not linear):
  Eigen::Vector3d midpoint = global.rowwise().sum() / 2;
  for (auto& g : find_result) {
    auto [distance, parameter] = g.first->Project(midpoint);
    if (distance(0, 0) > 1e-5) {
      // this is the bspline curve...
      // -> check what happens if we retrieve the geometry that lies slightly
      // outside the parameter bounds of the bspline curve:
      Eigen::MatrixXd local(1, 1);
      local(0, 0) = g.second(0, 0) - (g.second(0, 1) - g.second(0, 0)) / 100.;

      auto temp = g.first->Global(local);
      EXPECT_TRUE(g.first->IsInBoundingBox(temp)[0]);
      auto temp_result = model->FindGeometries(1, temp);
      EXPECT_EQ(temp_result.size(), 0);
    }
  }
}

// This test checks that the number of faces/edges is correct for two cubes
// that share a face.
TEST(occt, doubleCubeTest) {
  auto model = LoadModel("double_cube.brep");
  EXPECT_EQ(model->NumGeometries(1), 20);
  EXPECT_EQ(model->NumGeometries(2), 11);
}

}  // namespace lf::brep::occt::test
