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
      p_t{&circle, Eigen::RowVector2d(base::kPi / 2., -base::kPi / 2.)},
      p_t{&circle, Eigen::RowVector2d(-base::kPi / 2., 0.0)}};
  geometries.push_back(std::make_unique<BrepTriaTransfinite>(
      curves, std::array<bool, 3>{false, false, false}));

  auto qr = quad::make_QuadRule(base::RefEl::kTria(), 20);
  quad::QuadRuleCache qr_cache;
  auto qr_provider = [&](base::RefEl ref_el) {
    // return qr_cache.Get(ref_el, 20);
    if (ref_el == base::RefEl::kTria()) {
      return quad::QuadRule(base::RefEl::kTria(), Eigen::Vector2d(0.3, 0.3),
                            Eigen::RowVectorXd::Constant(1, 1.), 1);
    } else
      return qr_cache.Get(ref_el, 20);
  };

  refinement::Hybrid2DRefinementPattern ref_pat(base::RefEl::kTria(),
                                                refinement::RefPat::rp_regular);

  for (auto& g : geometries) {
    /*geometry::test_utils::checkJacobian(*g, qr.Points(), 1e-5);
    geometry::test_utils::checkJacobianInverseGramian(*g, qr.Points(), 1e-8);
    geometry::test_utils::checkIntegrationElement(*g, qr.Points());
    geometry::test_utils::checkSubGeometry(*g, qr_provider);*/
    geometry::test_utils::checkChildGeometry(*g, ref_pat, qr_provider);
  }
}

}  // namespace lf::brep::geom::tests
