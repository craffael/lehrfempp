/**
 * @file
 * @brief Tests geometry objects
 * @author Anian Ruoss
 * @date   2018-10-27 15:57:17
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/geometry/test_utils/check_child_geometry.h>
#include <lf/geometry/test_utils/check_integration_element.h>
#include <lf/geometry/test_utils/check_jacobian.h>
#include <lf/geometry/test_utils/check_jacobian_inverse_gramian.h>
#include <lf/geometry/test_utils/check_sub_geometry.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>

namespace lf::geometry::test {

void runGeometryChecks(const lf::geometry::Geometry &geom,
                       const Eigen::MatrixXd &eval_points,
                       const double tolerance) {
  lf::geometry::test_utils::checkJacobian(geom, eval_points, tolerance);
  // JacobianInverseGramian is not defined for points
  if (geom.RefEl() != lf::base::RefEl::kPoint()) {
    lf::geometry::test_utils::checkJacobianInverseGramian(geom, eval_points);
  }
  lf::geometry::test_utils::checkIntegrationElement(geom, eval_points);
  lf::geometry::test_utils::checkSubGeometry(
      geom, [](lf::base::RefEl refEl) -> lf::quad::QuadRule {
        return lf::quad::make_QuadRule(refEl, 5);
      });
}

TEST(Geometry, Point) {
  Eigen::MatrixXd global_nodes_2D(2, 1);
  global_nodes_2D << 1, 2;
  Eigen::MatrixXd global_nodes_3D(3, 1);
  global_nodes_3D << 3, 4, 0;

  for (const auto &global_nodes : {global_nodes_2D, global_nodes_3D}) {
    lf::geometry::Point geom(global_nodes);

    // QuadRule is not implemented and coordinate values are irrelevant
    Eigen::MatrixXd points = Eigen::MatrixXd::Random(0, 3);
    runGeometryChecks(geom, points, 1e-9);

    // check that local nodes are mapped to global nodes
    EXPECT_TRUE(
        geom.Global(Eigen::Matrix<double, 0, 1>()).isApprox(global_nodes));
  }
}

TEST(Geometry, SegmentO1) {
  Eigen::MatrixXd global_nodes_2D(2, 2);
  global_nodes_2D << 1, 1, 0, 4;
  Eigen::MatrixXd global_nodes_3D(3, 2);
  global_nodes_3D << -1, 10, 1, 3, 2, 3;

  for (const auto &global_nodes : {global_nodes_2D, global_nodes_3D}) {
    lf::geometry::SegmentO1 geom(global_nodes);
    auto qr = lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), 5);
    runGeometryChecks(geom, qr.Points(), 1e-9);

    std::vector<lf::refinement::RefPat> segSymmetricRefPats = {
        lf::refinement::RefPat::rp_nil, lf::refinement::RefPat::rp_copy,
        lf::refinement::RefPat::rp_split};

    for (const auto &refPat : segSymmetricRefPats) {
      lf::geometry::test_utils::checkChildGeometryVolume(geom, refPat);
      lf::geometry::test_utils::checkChildGeometry(
          geom, lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat),
          [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
    }

    // check that local nodes are mapped to global nodes
    EXPECT_TRUE(geom.Global(lf::base::RefEl::kSegment().NodeCoords())
                    .isApprox(global_nodes));
  }
}

TEST(Geometry, SegmentO2) {
  Eigen::MatrixXd global_nodes_2D(2, 3);
  global_nodes_2D << 1, 2, 3, 1, 4, 2;
  Eigen::MatrixXd global_nodes_3D(3, 3);
  global_nodes_3D << -1, 3, 1, 1, 2, 3, -1, 1, 0;

  for (const auto &global_nodes : {global_nodes_2D, global_nodes_3D}) {
    lf::geometry::SegmentO2 geom(global_nodes);
    auto qr = lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), 5);
    runGeometryChecks(geom, qr.Points(), 1e-9);

    std::vector<lf::refinement::RefPat> segSymmetricRefPats = {
        lf::refinement::RefPat::rp_nil, lf::refinement::RefPat::rp_copy,
        lf::refinement::RefPat::rp_split};

    for (const auto &refPat : segSymmetricRefPats) {
      lf::geometry::test_utils::checkChildGeometryVolume(geom, refPat);
      lf::geometry::test_utils::checkChildGeometry(
          geom, lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat),
          [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
    }

    Eigen::MatrixXd local_nodes(1, 3);
    local_nodes << 0, 1, 0.5;
    EXPECT_TRUE(geom.Global(local_nodes).isApprox(global_nodes));
  }
}

TEST(Geometry, TriaO1) {
  Eigen::MatrixXd global_nodes_2D(2, 3);
  global_nodes_2D << 1, 4, 3, 1, 2, 5;
  Eigen::MatrixXd global_nodes_3D(3, 3);
  global_nodes_3D << -1, 2, 0, 0, -2, -7, 3, 3, 5;

  for (const auto &global_nodes : {global_nodes_2D, global_nodes_3D}) {
    lf::geometry::TriaO1 geom(global_nodes);
    auto qr = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 5);
    runGeometryChecks(geom, qr.Points(), 1e-9);

    std::vector<lf::refinement::RefPat> triaSymmetricRefPats = {
        lf::refinement::RefPat::rp_nil, lf::refinement::RefPat::rp_copy,
        lf::refinement::RefPat::rp_regular,
        lf::refinement::RefPat::rp_barycentric};

    for (const auto &refPat : triaSymmetricRefPats) {
      lf::geometry::test_utils::checkChildGeometryVolume(geom, refPat);
      lf::geometry::test_utils::checkChildGeometry(
          geom, lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat),
          [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
    }

    std::vector<lf::refinement::RefPat> triaAsymmetricRefPats = {
        lf::refinement::RefPat::rp_bisect, lf::refinement::RefPat::rp_trisect,
        lf::refinement::RefPat::rp_trisect_left,
        lf::refinement::RefPat::rp_quadsect};

    for (const auto &refPat : triaAsymmetricRefPats) {
      for (size_t anchor = 0; anchor < 3; ++anchor) {
        lf::geometry::test_utils::checkChildGeometryVolume(geom, refPat,
                                                           anchor);
        lf::geometry::test_utils::checkChildGeometry(
            geom,
            lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat,
                                                      anchor),
            [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
      }
    }

    // check that local nodes are mapped to global nodes
    EXPECT_TRUE(geom.Global(lf::base::RefEl::kTria().NodeCoords())
                    .isApprox(global_nodes));
  }
}

TEST(Geometry, TriaO2) {
  Eigen::MatrixXd global_nodes_2D(2, 6);
  global_nodes_2D << 1, 4, 3, 2, 4, 0, 2, 1, 5, 3, 4, 4;
  Eigen::MatrixXd global_nodes_3D(3, 6);
  global_nodes_3D << 2, 0, -3, 1, -1, 0, -2, 0, -3, -1.5, -1, -2.8, 0, 5, 1, 2,
      3, .2;

  for (const auto &global_nodes : {global_nodes_2D, global_nodes_3D}) {
    lf::geometry::TriaO2 geom(global_nodes);
    auto qr = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 5);
    runGeometryChecks(geom, qr.Points(), 1e-9);

    std::vector<lf::refinement::RefPat> triaSymmetricRefPats = {
        lf::refinement::RefPat::rp_nil, lf::refinement::RefPat::rp_copy,
        lf::refinement::RefPat::rp_regular,
        lf::refinement::RefPat::rp_barycentric};

    for (const auto &refPat : triaSymmetricRefPats) {
      lf::geometry::test_utils::checkChildGeometryVolume(geom, refPat);
      lf::geometry::test_utils::checkChildGeometry(
          geom, lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat),
          [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
    }

    std::vector<lf::refinement::RefPat> triaAsymmetricRefPats = {
        lf::refinement::RefPat::rp_bisect, lf::refinement::RefPat::rp_trisect,
        lf::refinement::RefPat::rp_trisect_left,
        lf::refinement::RefPat::rp_quadsect};

    for (const auto &refPat : triaAsymmetricRefPats) {
      for (size_t anchor = 0; anchor < 3; ++anchor) {
        lf::geometry::test_utils::checkChildGeometryVolume(geom, refPat,
                                                           anchor);
        lf::geometry::test_utils::checkChildGeometry(
            geom,
            lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat,
                                                      anchor),
            [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
      }
    }

    // check that local nodes are mapped to global nodes
    Eigen::MatrixXd local_nodes(2, 6);
    local_nodes << 0, 1, 0, 0.5, 0.5, 0, 0, 0, 1, 0, 0.5, 0.5;
    EXPECT_TRUE(geom.Global(local_nodes).isApprox(global_nodes));
  }
}

TEST(Geometry, QuadO1) {
  Eigen::MatrixXd global_nodes_2D(2, 4);
  global_nodes_2D << -1, 3, 2, 1, -2, 0, 2, 1;
  Eigen::MatrixXd global_nodes_3D(3, 4);
  global_nodes_3D << 4, 5, -3, -5, -2, 1, 3, -3, -2, -3, 1, 3;

  for (const auto &global_nodes : {global_nodes_2D, global_nodes_3D}) {
    lf::geometry::QuadO1 geom(global_nodes);
    auto qr = lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 5);
    runGeometryChecks(geom, qr.Points(), 1e-9);

    std::vector<lf::refinement::RefPat> quadSymmetricRefPats = {
        lf::refinement::RefPat::rp_nil, lf::refinement::RefPat::rp_copy,
        lf::refinement::RefPat::rp_regular,
        lf::refinement::RefPat::rp_barycentric};

    for (const auto &refPat : quadSymmetricRefPats) {
      lf::geometry::test_utils::checkChildGeometryVolume(geom, refPat);
      lf::geometry::test_utils::checkChildGeometry(
          geom, lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat),
          [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
    }

    std::vector<lf::refinement::RefPat> triaAsymmetricRefPats = {
        lf::refinement::RefPat::rp_split, lf::refinement::RefPat::rp_bisect,
        lf::refinement::RefPat::rp_trisect, lf::refinement::RefPat::rp_quadsect,
        lf::refinement::RefPat::rp_threeedge};

    for (const auto &refPat : triaAsymmetricRefPats) {
      for (size_t anchor = 0; anchor < 4; ++anchor) {
        lf::geometry::test_utils::checkChildGeometryVolume(geom, refPat,
                                                           anchor);
        lf::geometry::test_utils::checkChildGeometry(
            geom,
            lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat,
                                                      anchor),
            lf::quad::make_QuadRuleNodal);
      }
    }

    // check that local nodes are mapped to global nodes
    EXPECT_TRUE(geom.Global(lf::base::RefEl::kQuad().NodeCoords())
                    .isApprox(global_nodes));
  }
}

TEST(Geometry, QuadO2) {
  Eigen::MatrixXd global_nodes_2D(2, 8);
  global_nodes_2D << 3, 7, 4, 1, 5, 5.4, 2.5, 1.8, 1, 3, 7, 8, 2.5, 5.9, 6.9,
      4.1;
  Eigen::MatrixXd global_nodes_3D(3, 8);
  global_nodes_3D << 4, 5, -3, -5, 4.499, 0.999, -4.01, -0.499, -2, 1, 3, -3,
      -0.5009, 2.01, 0.009, -2.501, -2, -3, 1, 3, -2.4999, -1.01, 2.01, .499;

  for (const auto &global_nodes : {global_nodes_2D, global_nodes_3D}) {
    lf::geometry::QuadO2 geom(global_nodes);
    auto qr = lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 5);
    runGeometryChecks(geom, qr.Points(), 1e-9);

    std::vector<lf::refinement::RefPat> quadSymmetricRefPats = {
        lf::refinement::RefPat::rp_nil, lf::refinement::RefPat::rp_copy,
        lf::refinement::RefPat::rp_regular,
        lf::refinement::RefPat::rp_barycentric};

    for (const auto &refPat : quadSymmetricRefPats) {
      lf::geometry::test_utils::checkChildGeometryVolume(geom, refPat);
      lf::geometry::test_utils::checkChildGeometry(
          geom, lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat),
          [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
    }

    std::vector<lf::refinement::RefPat> triaAsymmetricRefPats = {
        lf::refinement::RefPat::rp_split, lf::refinement::RefPat::rp_bisect,
        lf::refinement::RefPat::rp_trisect, lf::refinement::RefPat::rp_quadsect,
        lf::refinement::RefPat::rp_threeedge};

    for (const auto &refPat : triaAsymmetricRefPats) {
      for (size_t anchor = 0; anchor < 4; ++anchor) {
        lf::geometry::test_utils::checkChildGeometryVolume(geom, refPat,
                                                           anchor);
        lf::geometry::test_utils::checkChildGeometry(
            geom,
            lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat,
                                                      anchor),
            lf::quad::make_QuadRuleNodal);
      }
    }

    // check that local nodes are mapped to global nodes
    Eigen::MatrixXd local_nodes(2, 8);
    local_nodes << 0, 1, 1, 0, .5, 1, .5, 0, 0, 0, 1, 1, 0, .5, 1, 0.5;
    EXPECT_TRUE(geom.Global(local_nodes).isApprox(global_nodes));
  }
}

}  // namespace lf::geometry::test
