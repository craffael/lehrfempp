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
#include <lf/brep/test_utils/curve_straight_segment.h>
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

/**
 * @brief Here we make a triangle that occupies the space [0,1]^2\UnitCircle.
 *        This is a quite hard test case because the determinant of the mapping
 *        becomes extremely small/zero in large part of the domain.
 *        The "problem" is that the two straight line segments are tangential to
 * the circle arc and hence form a cusp.
 */
TEST(lf_brep_geom, BrepTriaTransfiniteCircleCutoutTest) {
  Eigen::Vector3d origin(0, 0, 0);
  test_utils::CurveCircle circle(origin, 1);

  test_utils::CurveStraightLine line0({1, 0, 0}, {1, 1, 0});
  test_utils::CurveStraightLine line1({0, 1, 0}, {1, 1, 0});

  using p_t = std::pair<const interface::BrepCurve*, Eigen::RowVector2d>;
  std::array<p_t, 3> curves{
      p_t{&line0, Eigen::RowVector2d(0., 1.)},
      p_t{&line1, Eigen::RowVector2d(1., 0.)},
      p_t{&circle, Eigen::RowVector2d(base::kPi / 2., 0.0)}};
  BrepTriaTransfinite tria(curves, {false, false, false});

  auto qr = quad::make_QuadRule(base::RefEl::kTria(), 40);

  auto jac = tria.Jacobian(qr.Points());
  for (int i = 0; i < qr.NumPoints(); ++i) {
    auto t = (jac.block(0, 2 * i, 3, 2).transpose() * jac.block(0, 2 * i, 3, 2))
                 .determinant();
    EXPECT_GT(t, 0);
  }

  // Make sure all the points lie within the expected area:
  // (Due to roundoff errors, the tolerance is quite low.)
  auto global = tria.Global(qr.Points());
  for (int i = 0; i < global.cols(); ++i) {
    EXPECT_GT(global.col(i).norm(), 1 - 1e-2);
    EXPECT_LT(global(0, i), 1 + 1e-2);
    EXPECT_LT(global(1, i), 1 + 1e-2);
    EXPECT_GT(global(0, i), -1e-2);
    EXPECT_GT(global(1, i), -1e-2);
  }

  double volume = tria.IntegrationElement(qr.Points()).dot(qr.Weights());
  EXPECT_NEAR(volume, 1 - base::kPi / 4., 0.02);

  // Note that the tolerance for the volume check is quite high.
  // Thats because numerical roundoff errors make it very hard to integrate the
  // volume correctly.
  // (Not even Mathematica can do it correctly)
}

/**
 * @brief This is a test case similar to BrepTriaTransfiniteCircleCutoutTest but
 * much more well-conditioned because the two straight line segments are not
 * tangential to the circle
 */
TEST(lf_brep_geom, BrepTriaTransfiniteCircleCutoutTest2) {
  Eigen::Vector3d origin(0, 0, 0);
  test_utils::CurveCircle circle(origin, 1);

  test_utils::CurveStraightLine line0({1, 0, 0}, {2, 2, 0});
  test_utils::CurveStraightLine line1({0, 1, 0}, {2, 2, 0});

  using p_t = std::pair<const interface::BrepCurve*, Eigen::RowVector2d>;
  std::array<p_t, 3> curves{
      p_t{&line0, Eigen::RowVector2d(0., 1.)},
      p_t{&line1, Eigen::RowVector2d(1., 0.)},
      p_t{&circle, Eigen::RowVector2d(base::kPi / 2., 0.0)}};
  BrepTriaTransfinite tria(curves, {false, false, false});

  auto qr = quad::make_QuadRule(base::RefEl::kTria(), 40);
  double volume = tria.IntegrationElement(qr.Points()).dot(qr.Weights());
  EXPECT_NEAR(volume, 4 - 2 - base::kPi / 4., 1e-6);
}

TEST(lf_brep_geom, BrepTriaTransfinitePerronnetCircleTest) {
  Eigen::Vector3d origin(0, 0, 0);
  test_utils::CurveCircle circle(origin, 1);

  std::vector<std::unique_ptr<geometry::Geometry>> geometries;
  using p_t = std::pair<const interface::BrepCurve*, Eigen::RowVector2d>;
  std::array<p_t, 3> curves{
      p_t{&circle, Eigen::RowVector2d(0., base::kPi / 2.)},
      p_t{&circle, Eigen::RowVector2d(base::kPi / 2., 3 * base::kPi / 2.)},
      p_t{&circle, Eigen::RowVector2d(-base::kPi / 2., 0.0)}};
  geometries.push_back(std::make_unique<BrepTriaTransfinitePerronnet>(
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

/**
 * @brief Here we make a triangle that occupies the space [0,1]^2\UnitCircle.
 *        This is a quite hard test case because the determinant of the mapping
 *        becomes extremely small/zero in large part of the domain.
 *        The "problem" is that the two straight line segments are tangential to
 *        the circle arc and hence form a cusp. However
 *        BrepTriaTransfinitePerronnet solves this problem much better than
 *        BrepTriaTransfinite.
 */
TEST(lf_brep_geom, BrepTriaTransfinitePerronnetCircleCutoutTest) {
  Eigen::Vector3d origin(0, 0, 0);
  test_utils::CurveCircle circle(origin, 1);

  test_utils::CurveStraightLine line0({1, 0, 0}, {1, 1, 0});
  test_utils::CurveStraightLine line1({0, 1, 0}, {1, 1, 0});

  using p_t = std::pair<const interface::BrepCurve*, Eigen::RowVector2d>;
  std::array<p_t, 3> curves{
      p_t{&line0, Eigen::RowVector2d(0., 1.)},
      p_t{&line1, Eigen::RowVector2d(1., 0.)},
      p_t{&circle, Eigen::RowVector2d(base::kPi / 2., 0.0)}};
  BrepTriaTransfinitePerronnet tria(curves, {false, false, false});

  auto qr = quad::make_QuadRule(base::RefEl::kTria(), 40);

  auto jac = tria.Jacobian(qr.Points());
  for (int i = 0; i < qr.NumPoints(); ++i) {
    auto t = (jac.block(0, 2 * i, 3, 2).transpose() * jac.block(0, 2 * i, 3, 2))
                 .determinant();
    EXPECT_GT(t, 0);
  }

  // Make sure all the points lie within the expected area:
  // (Due to roundoff errors, the tolerance is quite low.)
  auto global = tria.Global(qr.Points());
  for (int i = 0; i < global.cols(); ++i) {
    EXPECT_GT(global.col(i).norm(), 1 - 1e-2);
    EXPECT_LT(global(0, i), 1 + 1e-2);
    EXPECT_LT(global(1, i), 1 + 1e-2);
    EXPECT_GT(global(0, i), -1e-2);
    EXPECT_GT(global(1, i), -1e-2);
  }

  double volume = tria.IntegrationElement(qr.Points()).dot(qr.Weights());
  EXPECT_NEAR(volume, 1 - base::kPi / 4., 1e-6);
}

}  // namespace lf::brep::geom::tests
