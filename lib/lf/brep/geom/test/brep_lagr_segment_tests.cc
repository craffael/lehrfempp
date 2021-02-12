/**
 * @file
 * @brief Test implementation of BrepLagrSegment
 * @author Raffael Casagrande
 * @date   2021-01-29 10:59:52
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/brep/geom/geom.h>
#include <lf/brep/test_utils/curve_circle.h>
#include <lf/quad/quad.h>

#include "utils.h"

namespace lf::brep::geom::test {

template <class GEOM>
void CheckBrepLagrSegment(
    const BrepLagrSegment<GEOM>& geom,
    std::shared_ptr<const interface::BrepGeometry> curve) {
  // check ChildGeometry(codim=1)
  std::vector<std::unique_ptr<geometry::Geometry>> point_children =
      geom.ChildGeometry(SegmentRefPat{}, 1);
  Eigen::Matrix<double, 0, 1> zero_mat;
  auto make_11Mat = [](double value) {
    return (Eigen::Matrix<double, 1, 1>() << value).finished();
  };
  auto checkOnCurve = [&curve](Eigen::MatrixXd x) {
    for (int i = 0; i < x.cols(); ++i) {
      auto [dist, p] = curve->Project(x.col(i));
      ASSERT_LT(dist, 1e-6);
    }
  };
  ASSERT_EQ(point_children.size(), 3);
  ASSERT_TRUE(
      point_children[0]->Global(zero_mat).isApprox(geom.Global(make_11Mat(0))));
  checkOnCurve(point_children[2]->Global(zero_mat));
  ASSERT_TRUE(
      point_children[2]->Global(zero_mat).isApprox(geom.Global(make_11Mat(1))));

  // check ChildGeometry(codim=0)
  auto segment_children = geom.ChildGeometry(SegmentRefPat{}, 0);
  ASSERT_EQ(segment_children.size(), 2);
  ASSERT_TRUE(segment_children[0]
                  ->Global(make_11Mat(0))
                  .isApprox(geom.Global(make_11Mat(0))));
  ASSERT_TRUE(segment_children[1]
                  ->Global(make_11Mat(1))
                  .isApprox(geom.Global(make_11Mat(1))));
  for (int i = 0; i < 2; ++i) {
    checkOnCurve(segment_children[0]->Global(GEOM::LagrangeNodes()));
  }
}

template <class GEOM>
void CheckBRepLagrSegmentOnCircle() {
  // create a segment that spans the quarter of a circle:
  auto circle =
      std::make_shared<test_utils::CurveCircle>(Eigen::Vector3d(1, 2, 0), 1.);
  Eigen::MatrixXd nodes(3, GEOM::LagrangeNodes().cols());
  nodes.row(0).array() =
      (GEOM::LagrangeNodes().array() * base::kPi / 2).cos() + 1.;
  nodes.row(1).array() =
      (GEOM::LagrangeNodes().array() * base::kPi / 2).sin() + 2.;
  nodes.row(2).setZero();
  auto geom = std::make_unique<BrepLagrSegment<GEOM>>(GEOM{nodes}, circle, 0,
                                                      base::kPi / 2);

  // check the geometry
  CheckBrepLagrSegment(*geom, circle);
  // refine a few times to approximate the circle better::
  std::vector<std::unique_ptr<geometry::Geometry>> geoms;
  geoms.emplace_back(std::move(geom));
  auto qr = quad::make_QuadRule(base::RefEl::kSegment(), 10);
  double arc_length;
  for (int level = 1; level < 10; ++level) {
    std::vector<std::unique_ptr<geometry::Geometry>> next_geoms;
    arc_length = 0.;
    for (auto& g : geoms) {
      auto gs = g->ChildGeometry(SegmentRefPat{}, 0);
      for (auto&& g2 : gs) {
        arc_length += g2->IntegrationElement(qr.Points()).dot(qr.Weights());
        CheckBrepLagrSegment(dynamic_cast<const BrepLagrSegment<GEOM>&>(*g2),
                             circle);
        next_geoms.emplace_back(std::move(g2));
      }
    }
    geoms = std::move(next_geoms);

    fmt::print("Arc Length on Level {} : {}\n", level, arc_length);
  }

  ASSERT_NEAR(arc_length, base::kPi / 2, 1e-4);
}

TEST(lf_brep_geom, BrepLagrSegmentO1Test) {
  CheckBRepLagrSegmentOnCircle<geometry::SegmentO1>();
}

TEST(lf_brep_geom, BrepLagrSegmentO2Test) {
  CheckBRepLagrSegmentOnCircle<geometry::SegmentO2>();
}
}  // namespace lf::brep::geom::test
