/**
 * @file
 * @brief Test the OcctBrepModel class
 * @author Raffael Casagrande
 * @date   2020-11-11 02:56:08
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/brep/test_utils/check_brep_geometry.h>

#include "utils.h"

namespace lf::brep::occt::test {

/**
 * @brief Checks that all Curves returned by FindCurvesMulti are contained in
 * the result of FindCurves()
 */
std::vector<std::pair<interface::BrepCurve const*, Eigen::RowVectorXd>>
FindCurves(const OcctBrepModel& model, const Eigen::Matrix3Xd& global) {
  auto all = model.FindCurvesMulti(global);
  for (auto& [c, param] : all) {
    EXPECT_TRUE(c->GlobalMulti(param).isApprox(global));
    test_utils::CheckBrepGeometry(*c, param);
  }
  for (int i = 0; i < global.cols(); ++i) {
    auto curves = model.FindCurves(global.col(i));
    EXPECT_TRUE(std::all_of(all.begin(), all.end(), [&](const auto& a) {
      return std::any_of(curves.begin(), curves.end(),
                         [&](const auto& c) { return c.first == a.first; });
    }));
  }
  return all;
}

std::vector<std::pair<interface::BrepSurface const*, Eigen::Matrix2Xd>>
FindSurfaces(const OcctBrepModel& model, const Eigen::Matrix3Xd& global) {
  auto all = model.FindSurfacesMulti(global);
  for (auto& [s, param] : all) {
    EXPECT_TRUE(s->GlobalMulti(param).isApprox(global));
    test_utils::CheckBrepGeometry(*s, param);
  }
  for (int i = 0; i < global.cols(); ++i) {
    auto surfaces = model.FindSurfaces(global.col(i));
    EXPECT_TRUE(std::all_of(all.begin(), all.end(), [&](const auto& a) {
      return std::any_of(surfaces.begin(), surfaces.end(),
                         [&](const auto& c) { return c.first == a.first; });
    }));
  }
  return all;
}

TEST(occtDeathTest, fileDoesntExist) {
  EXPECT_DEATH(OcctBrepModel("invalid.brep2"), "Could not open file");
}

TEST(occt, rotatedCubeTest) {
  auto model = LoadModel("cube_rotated.brep");
  EXPECT_EQ(model->NumCurves(), 12);
  EXPECT_EQ(model->NumSurfaces(), 6);

  // retrieve edge from [20,0,0] -> [20,10,0]:
  Eigen::MatrixXd global(3, 2);
  // clang-format off
  global << 20,  20,
            0,  10,
            0,  0;
  // clang-format on
  auto find_result = FindCurves(*model, global);
  EXPECT_EQ(find_result.size(), 1);
  EXPECT_TRUE(find_result[0].first);

  // try to retrieve an edge that is outside of the parameter bounds:
  // clang-format off
  global << 20, 20,
             0, 20,
             0,  0;
  // clang-format on
  find_result = FindCurves(*model, global);
  EXPECT_EQ(find_result.size(), 0);

  // retrieve face
  global.resize(3, 4);
  // clang-format off
   global << 20, 20, 20, 20,
             0,  0, 10, 10,
            10,  0,  0, 10;
  // clang-format on
  auto find_surfaces = FindSurfaces(*model, global);
  EXPECT_EQ(find_surfaces.size(), 1);
  EXPECT_TRUE(find_surfaces[0].first);
  EXPECT_EQ(find_surfaces[0].second.rows(), 2);
  EXPECT_EQ(find_surfaces[0].second.cols(), 4);
  EXPECT_TRUE(find_surfaces[0]
                  .first->GlobalMulti(find_surfaces[0].second)
                  .isApprox(global));

  // try to retrieve a face that is outside the parameter bounds:
  // clang-format off
  global(1,2) = 20;
  find_surfaces = FindSurfaces(*model, global);
  EXPECT_EQ(find_surfaces.size(), 0);
}

// Here we check that FindGeometries() gives the correct result also if the
// global point lies on the curve but the parameter is out of bounds.
TEST(occt, bSpline2dTest) {
  auto model = LoadModel("bspline_2d.brep");
  EXPECT_EQ(model->NumCurves(), 2);
  EXPECT_EQ(model->NumSurfaces(), 1);

  // get both two curves:
  Eigen::MatrixXd global(3, 2);
  // clang-format off
  global << -10, -13,
            -10, -6,
              0,  0;
  // clang-format on
  auto find_result = FindCurves(*model, global);
  EXPECT_EQ(find_result.size(), 2);

  // find the b-spline curve (which is not linear):
  Eigen::Vector3d midpoint = global.rowwise().sum() / 2;
  interface::BrepCurve const* bspline;
  Eigen::RowVectorXd bspline_param;
  for (auto& g : find_result) {
    auto [distance, parameter] = g.first->Project(midpoint);
    if (distance > 1e-5) {
      // this is the bspline curve...
      bspline = g.first;
      bspline_param = g.second;
    }
  }
  ASSERT_TRUE(bspline);
  // -> check what happens if we retrieve the geometry that lies slightly
  // outside the parameter bounds of the bspline curve:
  Eigen::MatrixXd local(1, 1);
  local(0, 0) =
      bspline_param(0, 0) - (bspline_param(0, 1) - bspline_param(0, 0)) / 100.;

  auto temp = bspline->GlobalMulti(local);
  EXPECT_TRUE(bspline->IsInBoundingBoxMulti(temp)[0]);
  auto temp_result = FindCurves(*model, temp);
  EXPECT_EQ(temp_result.size(), 0);

  // find curves which go through (-13,-6,0) + eps:
  auto find_result_single =
      model->FindCurves(Eigen::Vector3d(-13 + 1e-7, -6, 0));
  ASSERT_EQ(find_result_single.size(), 2);

  // retrieve the surface:
  auto find_surface = FindSurfaces(*model, global);
  EXPECT_EQ(find_surface.size(), 1);

  // try to retrieve the surface with a point that lies outside:
  find_surface = FindSurfaces(*model, Eigen::Vector3d(-4, -10, 0));
  EXPECT_TRUE(find_surface.empty());

  // try to find surface with a point that lies inside and one outside:
  global.col(1) = Eigen::Vector3d(-4, -10, 0);
  find_surface = FindSurfaces(*model, global);
  EXPECT_TRUE(find_surface.empty());
}

// This test checks that the number of faces/edges is correct for two cubes
// that share a face.
TEST(occt, doubleCubeTest) {
  auto model = LoadModel("double_cube.brep");
  EXPECT_EQ(model->NumCurves(), 20);
  EXPECT_EQ(model->NumSurfaces(), 11);

  // try to retrieve the edge that goes through (0,1,0) and (1,1,0)
  Eigen::MatrixXd global(3, 2);
  global.col(0) = Eigen::Vector3d(0, 1, 0);
  global.col(1) = Eigen::Vector3d(1, 1, 0);
  auto find_result = FindCurves(*model, global);
  EXPECT_EQ(find_result.size(), 1);

  // try to retrieve an edge that goes through (0,1,0) and (1,1,1):
  global.col(1) = Eigen::Vector3d(1, 1, 1);
  find_result = FindCurves(*model, global);
  EXPECT_TRUE(find_result.empty());

  // retrieve the face that goes through (0,1,0) and (1,1,1):
  auto find_surf = FindSurfaces(*model, global);
  EXPECT_EQ(find_surf.size(), 1);
}

}  // namespace lf::brep::occt::test
