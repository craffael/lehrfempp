/**
 * @file
 * @brief Check that the Occt
 * @author Raffael Casagrande
 * @date   2020-11-23 09:43:11
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/brep/occt/occt.h>

#include "lf/brep/test_utils/check_brep_geometry.h"
#include "utils.h"

namespace lf::brep::occt::test {

auto make11Matrix = [](double value) {
  return (Eigen::VectorXd(1) << value).finished();
};

TEST(occt, curveProjectionCube) {
  auto model = LoadModel("cube.brep");

  // get curve that goes through (10,10,10) and (10,10,0):
  Eigen::MatrixXd global(3, 2);
  global.col(0) = Eigen::Vector3d(10, 10, 10);
  global.col(1) = Eigen::Vector3d(10, 10, 0);
  auto find_curve = model->FindCurvesMulti(global);
  EXPECT_EQ(find_curve.size(), 1);
  auto curve = find_curve[0].first;
  ASSERT_EQ(curve->Periods()(0), 0);

  // project point (10,10,10) onto curve
  auto [dist, local] = curve->Project(Eigen::Vector3d(10, 10, 10));
  ASSERT_GE(dist, 0);
  ASSERT_LT(dist, 1e-6);
  ASSERT_TRUE(curve->Global(local).isApprox(Eigen::Vector3d(10, 10, 10)));

  // project point (10,10,20) onto curve:
  std::tie(dist, local) = curve->Project(Eigen::Vector3d(10, 10, 20));
  ASSERT_LE(dist, 1e-6);
  ASSERT_FALSE(curve->IsInside(local));

  // project point 12,10,5 onto curve:
  std::tie(dist, local) = curve->Project(Eigen::Vector3d(12, 10, 5));
  ASSERT_LE(std::abs(dist - 2), 1e-6);
  ASSERT_TRUE(curve->IsInside(local));

  // Check rest of the geometry:
  // -> First find the the local coordinates:
  auto [dista, parama] = curve->Project(Eigen::Vector3d(10, 10, 10));
  auto [distb, paramb] = curve->Project(Eigen::Vector3d(10, 10, 0));
  test_utils::CheckBrepGeometry(
      *curve, Eigen::RowVectorXd::LinSpaced(11, parama(0), paramb(0)));
}

TEST(occt, curveProjectionBSpline) {
  auto model = LoadModel("bspline_2d.brep");

  // find b-spline curve:
  auto find_curve = model->FindCurves(Eigen::Vector3d(-10, -10, 0));
  std::shared_ptr<interface::BrepGeometry const> curve;
  double curve_param;
  for (auto [c, param] : find_curve) {
    if (c->Project(Eigen::Vector3d(-11.5, -8, 0)).first > 1e-3) {
      curve = c;
      curve_param = param;
      break;
    }
  }
  ASSERT_TRUE(curve);
  ASSERT_EQ(curve->Periods()(0), 0);

  // project point (-13,-6,0) (other endpoint)
  auto [dist, local] = curve->Project(Eigen::Vector3d(-13, -6, 0));
  ASSERT_LT(dist, 1e-6);
  ASSERT_TRUE(curve->IsInside(local));

  // project point on the outside:
  std::tie(dist, local) = curve->Project(Eigen::Vector3d(-15, -7, 0));
  ASSERT_GT(dist, 1);
  ASSERT_TRUE(
      std::abs(dist -
               (curve->Global(local) - Eigen::Vector3d(-15, -7, 0)).norm()) <
      1e-6);
  ASSERT_TRUE(curve->IsInside(local));
  // std::cout << curve->Global(local) << std::endl;

  // another point:
  Eigen::Vector3d p(-13 - 1e-7, -6, 0);
  std::tie(dist, local) = curve->Project(p);

  ASSERT_TRUE(std::abs(dist - (curve->Global(local) - p).norm()) < 1e-6);
  // std::cout << curve->Global(local) << std::endl;
  ASSERT_LT(dist, 1e-7);
  ASSERT_TRUE(curve->IsInside(local));

  // now map a point that lies slightly inside the bspline and project it back:
  ASSERT_TRUE(curve->IsInside(make11Matrix(curve_param + 1e-3)));
  p = curve->Global(make11Matrix(curve_param + 1e-3));
  std::tie(dist, local) = curve->Project(p);
  ASSERT_LT(dist, 0.1);

  // now map a point that lies slightly outside the bpsline and project it back:
  ASSERT_FALSE(curve->IsInside(make11Matrix(curve_param - 1e-3)));
  p = curve->Global(make11Matrix(curve_param - 1e-3));
  std::tie(dist, local) = curve->Project(p);
  EXPECT_GT(dist, 1);  // in this case we have no guarantee, and here the
                       // distance is actually much bigger
}

TEST(occt, curveCircleTest) {
  // load a circle with origin=(1,2,0) and radius=2
  auto model = LoadModel("circle.brep");

  auto find_curve = model->FindCurves(Eigen::Vector3d(3, 2, 0));
  ASSERT_EQ(find_curve.size(), 1);
  auto circle = find_curve[0].first;
  auto loc0 = find_curve[0].second;
  Eigen::Vector3d origin(1, 2, 0);

  ASSERT_EQ(circle->DimLocal(), 1);
  ASSERT_EQ(circle->DimGlobal(), 3);
  ASSERT_EQ(circle->Periods()(0), 2 * base::kPi);

  // Test global:
  Eigen::RowVectorXd local =
      Eigen::RowVectorXd::LinSpaced(10, 0, 2 * base::kPi);
  auto global = circle->Global(local);
  for (int i = 0; i < global.cols(); ++i) {
    ASSERT_NEAR((global.col(i) - origin).norm(), 2, 1e-7);
  }

  // test bounding box:
  auto in_box = circle->IsInBoundingBox(global);
  for (int i = 0; i < global.cols(); ++i) {
    ASSERT_TRUE(in_box[i]);
  }

  ASSERT_TRUE(circle->IsInBoundingBox(origin)[0]);
  ASSERT_FALSE(circle->IsInBoundingBox(Eigen::Vector3d(4, 2, 0))[0]);

  // test IsInside:
  ASSERT_TRUE(circle->IsInside(make11Matrix(0.)));
  ASSERT_TRUE(circle->IsInside(make11Matrix(2 * base::kPi)));

  // even though circle is periodic, not all values are inside:
  ASSERT_FALSE(circle->IsInside(make11Matrix(2.1 * base::kPi)));
  ASSERT_FALSE(circle->IsInside(make11Matrix(-0.1)));

  // test projection by projecting points on the circle:
  for (int i = 0; i < global.cols() - 1;
       ++i) {  // leave out the last point as it projects onto 0
    auto [dist, p] = circle->Project(global.col(i));
    ASSERT_LT(dist, 1e-7);
    ASSERT_NEAR(p(0, 0), local[i], 1e-7);
    ASSERT_TRUE(circle->IsInside(p));
  }

  // test projection by projecting a point from outside:
  Eigen::Vector3d global2(3, 4, 5);
  Eigen::Vector3d global2_proj;
  global2_proj.topRows(2) =
      (global2.topRows(2) - origin.topRows(2)).normalized() * 2 +
      origin.topRows(2);
  global2_proj.z() = 0.;
  auto [dist, p] = circle->Project(global2);
  ASSERT_TRUE(circle->IsInside(p));
  ASSERT_NEAR(dist, (global2 - global2_proj).norm(), 1e-7);
  ASSERT_TRUE(circle->Global(p).isApprox(global2_proj));
}

}  // namespace lf::brep::occt::test
