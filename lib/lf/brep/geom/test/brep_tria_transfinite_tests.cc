/**
 * @file
 * @brief Test implementation of BrepTriaTransfinite
 * @author Raffael Casagrande
 * @date   2021-02-03 06:02:32
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/brep/geom/geom.h>
#include <lf/brep/test_utils/curve_circle.h>
#include <lf/geometry/test_utils/test_utils.h>
#include <lf/quad/quad.h>

namespace lf::brep::geom::tests {

TEST(lf_brep_geom, BrepTriaTransfiniteCircleTest) {
  Eigen::Vector3d origin(0, 0, 0);
  test_utils::CurveCircle circle(origin, 1);

  std::vector<std::unique_ptr<geometry::Geometry>> geometries;
  using p_t = std::pair<const interface::BrepCurve*, Eigen::RowVector2d>;
  std::array<p_t, 3> curves{
      p_t{&circle, Eigen::RowVector2d(0., base::kPi / 2.)},
      p_t{&circle, Eigen::RowVector2d(base::kPi / 2., 3 * base::kPi / 2.)},
      p_t{&circle, Eigen::RowVector2d(-base::kPi / 2., 0.0)}};
  geometries.push_back(std::make_unique<BrepTriaTransfinite>(
      curves, std::array<bool, 3>{false, false, false}));

  auto qr = quad::make_QuadRule(base::RefEl::kTria(), 20);
  auto qr_segment = quad::make_QuadRule(base::RefEl::kSegment(), 20);
  quad::QuadRuleCache qr_cache;
  auto qr_provider = [&](base::RefEl ref_el) {
    // return qr_cache.Get(ref_el, 20);
    return quad::make_QuadRuleNodal(ref_el);
  };

  refinement::Hybrid2DRefinementPattern ref_pat(base::RefEl::kTria(),
                                                refinement::RefPat::rp_regular);
  refinement::Hybrid2DRefinementPattern ref_pat_segment(
      base::RefEl::kSegment(), refinement::RefPat::rp_split);

  for (int level = 0; level < 2; ++level) {
    double volume = 0.0;
    std::vector<std::unique_ptr<geometry::Geometry>> next_geoms;
    for (auto& g : geometries) {
      geometry::test_utils::checkJacobian(*g, qr.Points(), 1e-5);
      geometry::test_utils::checkJacobianInverseGramian(*g, qr.Points(), 1e-8);
      geometry::test_utils::checkIntegrationElement(*g, qr.Points());
      geometry::test_utils::checkSubGeometry(*g, qr_provider);
      geometry::test_utils::checkChildGeometry(*g, ref_pat, qr_provider);

      auto global = g->Global(qr.Points());
      for (int i = 0; i < global.cols(); ++i) {
        ASSERT_NEAR(global(2, i), 0, 1e-6);
        ASSERT_LT(global.col(i).topRows(2).norm(), 1);
      }

      // check the sub geometries of codim=1:
      for (int i = 0; i < 3; ++i) {
        auto sg = g->SubGeometry(1, i);
        geometry::test_utils::checkJacobian(*sg, qr_segment.Points(), 1e-5);
        geometry::test_utils::checkJacobianInverseGramian(*sg,
                                                          qr_segment.Points());
        geometry::test_utils::checkIntegrationElement(*sg, qr_segment.Points());
        geometry::test_utils::checkSubGeometry(*sg, qr_provider);
        geometry::test_utils::checkChildGeometry(*sg, ref_pat_segment,
                                                 qr_provider);
      }
      volume += g->IntegrationElement(qr.Points()).dot(qr.Weights());

      for (auto&& sg : g->ChildGeometry(ref_pat, 0)) {
        next_geoms.push_back(std::move(sg));
      }
    }
    fmt::print("Area on Level {} : {}\n", level, volume);
    geometries = std::move(next_geoms);
    EXPECT_NEAR(volume, base::kPi, 1e-6);
  }
}

}  // namespace lf::brep::geom::tests
