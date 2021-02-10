/**
 * @file
 * @brief Test BrepSegment class
 * @author Raffael Casagrande
 * @date   2021-02-01 04:01:55
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/brep/geom/geom.h>
#include <lf/brep/test_utils/curve_circle.h>
#include <lf/geometry/test_utils/test_utils.h>
#include <lf/quad/quad.h>

#include "utils.h"

namespace lf::brep::geom::test {

TEST(lf_brep_geom, BrepSegmentCircleApprox) {
  auto circle =
      std::make_shared<test_utils::CurveCircle>(Eigen::Vector3d(0, 0, 0), 1.);

  std::vector<std::unique_ptr<geometry::Geometry>> geometries;
  geometries.push_back(std::make_unique<BrepCurve>(
      circle, Eigen::RowVector2d(0, base::kPi / 2)));
  auto qr = quad::make_QuadRule(base::RefEl::kSegment(), 20);

  quad::QuadRuleCache qr_cache;
  auto qr_provider = [&](auto ref_el) {
    return qr_cache.Get(base::RefEl::kSegment(), 20);
  };

  for (int level = 0; level < 4; ++level) {
    std::vector<std::unique_ptr<geometry::Geometry>> next_geoms;
    double arc_length = 0.;
    for (auto& g : geometries) {
      geometry::test_utils::checkJacobian(*g, qr.Points(), 1e-5);
      geometry::test_utils::checkJacobianInverseGramian(*g, qr.Points());
      geometry::test_utils::checkIntegrationElement(*g, qr.Points());
      geometry::test_utils::checkSubGeometry(*g, qr_provider);
      geometry::test_utils::checkChildGeometry(*g, SegmentRefPat(),
                                               qr_provider);

      // get volume:
      arc_length += g->IntegrationElement(qr.Points()).dot(qr.Weights());

      // check that all points lie on the circle:
      auto mapped_points = g->Global(qr.Points());
      for (int i = 0; i < mapped_points.cols(); ++i) {
        EXPECT_NEAR(1, mapped_points.col(i).norm(), 1e-6);
      }

      // refine the segment:
      auto temp = g->ChildGeometry(SegmentRefPat(), 0);
      for (auto& t : temp) {
        next_geoms.push_back(std::move(t));
      }
    }
    EXPECT_NEAR(arc_length, base::kPi / 2, 1e-6);
    geometries = std::move(next_geoms);
  }
}

}  // namespace lf::brep::geom::test
