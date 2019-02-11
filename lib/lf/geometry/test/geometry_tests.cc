/**
 * @file
 * @brief tests for geometry objects
 * @author Anian Ruoss
 * @date   2018-10-27 15:57:17
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/geometry/geometry.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>

/**
 * Checks if Jacobian() is implemented correctly by comparing it with the
 * symmetric difference quotient approximation
 */
void checkJacobian(const lf::geometry::Geometry &geom,
                   const Eigen::MatrixXd &eval_points, double tolerance) {
  const double h = 1e-6;

  const size_t num_points = eval_points.cols();
  const size_t dim_local = geom.DimLocal();
  const size_t dim_global = geom.DimGlobal();

  Eigen::MatrixXd jacobians = geom.Jacobian(eval_points);

  EXPECT_EQ(jacobians.rows(), dim_global) << "Jacobian has " << jacobians.rows()
                                          << " rows instead of " << dim_global;
  EXPECT_EQ(jacobians.cols(), num_points * dim_local)
      << "Jacobian has " << jacobians.cols() << " cols instead of "
      << num_points * dim_local;

  for (size_t j = 0; j < num_points; ++j) {
    auto point = eval_points.col(j);

    Eigen::MatrixXd jacobian =
        jacobians.block(0, j * dim_local, dim_global, dim_local);
    Eigen::MatrixXd approx_jacobian =
        Eigen::MatrixXd::Zero(dim_global, dim_local);

    for (size_t i = 0; i < dim_local; ++i) {
      Eigen::VectorXd h_vec = Eigen::VectorXd::Zero(dim_local);
      h_vec(i) = h;

      // approximate gradient with symmetric difference quotient
      approx_jacobian.col(i) =
          (geom.Global(point + h_vec) - geom.Global(point - h_vec)) / (2. * h);
    }

    EXPECT_TRUE(jacobian.isApprox(approx_jacobian, tolerance))
        << "Jacobian incorrect at point " << point;
  }
}

/**
 * Checks if JacobianInverseGramian() is implemented correctly under the
 * assumption that Jacobian() is correct
 */
void checkJacobianInverseGramian(const lf::geometry::Geometry &geom,
                                 const Eigen::MatrixXd &eval_points) {
  const size_t num_points = eval_points.cols();
  const size_t dim_local = geom.DimLocal();
  const size_t dim_global = geom.DimGlobal();

  Eigen::MatrixXd jacobians = geom.Jacobian(eval_points);
  Eigen::MatrixXd jacInvGrams = geom.JacobianInverseGramian(eval_points);

  EXPECT_EQ(jacInvGrams.rows(), dim_global)
      << "JacobianInverseGramian has " << jacInvGrams.rows()
      << " rows instead of " << dim_global;
  EXPECT_EQ(jacInvGrams.cols(), num_points * dim_local)
      << "JacobianInverseGramian has " << jacInvGrams.cols()
      << " cols instead of " << num_points * dim_local;

  for (int j = 0; j < num_points; ++j) {
    Eigen::MatrixXd jacInvGram =
        jacInvGrams.block(0, j * dim_local, dim_global, dim_local);
    Eigen::MatrixXd jacobian =
        jacobians.block(0, j * dim_local, dim_global, dim_local);

    EXPECT_TRUE(jacInvGram.isApprox(
        jacobian * (jacobian.transpose() * jacobian).inverse()))
        << "JacobianInverseGramian incorrect at point " << eval_points.col(j);
  }
}

/**
 * Checks if IntegrationElement() is implemented correctly under the assumption
 * that Jacobian() is correct
 */
void checkIntegrationElement(const lf::geometry::Geometry &geom,
                             const Eigen::MatrixXd &eval_points) {
  const size_t num_points = eval_points.cols();
  const size_t dim_local = geom.DimLocal();
  const size_t dim_global = geom.DimGlobal();

  Eigen::MatrixXd jacobians = geom.Jacobian(eval_points);
  Eigen::VectorXd integrationElements = geom.IntegrationElement(eval_points);

  EXPECT_EQ(integrationElements.rows(), num_points)
      << "IntegrationElement has " << integrationElements.rows()
      << " rows instead of " << num_points;
  EXPECT_EQ(integrationElements.cols(), 1)
      << "IntegrationElement has " << integrationElements.cols()
      << " cols instead of " << 1;

  for (int j = 0; j < num_points; ++j) {
    Eigen::MatrixXd jacobian =
        jacobians.block(0, j * dim_local, dim_global, dim_local);

    const double integrationElement = integrationElements(j);
    const double approx_integrationElement =
        std::sqrt((jacobian.transpose() * jacobian).determinant());

    EXPECT_FLOAT_EQ(integrationElement, approx_integrationElement)
        << "IntegrationElement incorrect at point " << eval_points.col(j);
  }
}

/**
 * Checks that SubGeometry and Geometry map the same nodes to the same points
 */
void checkSubGeometry(
    const lf::geometry::Geometry &geom,
    const std::function<lf::quad::QuadRule(lf::base::RefEl)> &qrProvider) {
  // nodeCoords is a (refEl.Dimension, refEl.NumNodes) matrix
  const auto refEl = geom.RefEl();
  const Eigen::MatrixXd &nodeCoords = refEl.NodeCoords();

  // iterate over all relative codimensions
  for (size_t codim = 0; codim <= geom.RefEl().Dimension(); ++codim) {
    // iterate over all subEntities in given codimension
    const auto numSubEntities = refEl.NumSubEntities(codim);
    for (size_t subEntity = 0; subEntity < numSubEntities; ++subEntity) {
      // subNodeCoords is a (subRefEl.Dimension, subRefEl.NumNodes) matrix
      auto subGeom = geom.SubGeometry(codim, subEntity);
      auto subRefEl = subGeom->RefEl();
      const Eigen::MatrixXd &subNodeCoords = subRefEl.NodeCoords();

      // iterate over all nodes of subEntity
      for (size_t subNode = 0; subNode < subRefEl.NumNodes(); ++subNode) {
        // map coordinates in subRefEl.Dimension to geom.DimGlobal
        auto globalCoordsFromSub = subGeom->Global(subNodeCoords.col(subNode));
        // get index of subSubEntity with respect to refEl
        int subSubIdx = refEl.SubSubEntity2SubEntity(
            codim, subEntity, geom.DimLocal() - codim, subNode);
        // map coordinates in RefEl.Dimension to geom.DimGlobal
        auto globalCoords = geom.Global(nodeCoords.col(subSubIdx));

        EXPECT_TRUE(globalCoordsFromSub.isApprox(globalCoords))
            << "Global mapping of subNode " << subNode << " of subEntity "
            << subEntity << " in relative codim " << codim
            << " differs from global mapping of node " << subSubIdx;
      }

      // check points inside subEntity
      if (subRefEl != lf::base::RefEl::kPoint()) {
        // select nodes of RefEl referenced by subEntity
        Eigen::MatrixXd mapNodeCoords(nodeCoords.rows(), subRefEl.NumNodes());
        for (size_t subNode = 0; subNode < subRefEl.NumNodes(); ++subNode) {
          const auto subSubIdx = refEl.SubSubEntity2SubEntity(
              codim, subEntity, geom.DimLocal() - codim, subNode);
          mapNodeCoords.col(subNode) = nodeCoords.col(subSubIdx);
        }

        auto points = qrProvider(subRefEl).Points();

        // compute mapping: alpha * subNodeCoords + beta = mapNodeCoords
        Eigen::MatrixXd paddedSubNodeCoords = Eigen::MatrixXd::Ones(
            subNodeCoords.rows() + 1, subNodeCoords.cols());
        paddedSubNodeCoords.block(0, 0, subNodeCoords.rows(),
                                  subNodeCoords.cols()) = subNodeCoords;
        Eigen::MatrixXd alphaBeta = paddedSubNodeCoords.transpose()
                                        .fullPivLu()
                                        .solve(mapNodeCoords.transpose())
                                        .transpose();

        // map points onto RefEl
        Eigen::MatrixXd paddedPoints =
            Eigen::MatrixXd::Ones(points.rows() + 1, points.cols());
        paddedPoints.block(0, 0, points.rows(), points.cols()) = points;
        const Eigen::MatrixXd mappedPoints = alphaBeta * paddedPoints;

        // map coordinates in subRefEl.Dimension to geom.DimGlobal
        auto globalPointsFromSub = subGeom->Global(points);
        // map coordinates in RefEl.Dimension to geom.DimGlobal
        auto globalPoints = geom.Global(mappedPoints);

        EXPECT_TRUE(globalPoints.isApprox(globalPointsFromSub))
            << "Global mapping of points " << points << " from subEntity "
            << subEntity << " in relative codim " << codim
            << " differs from global mapping";
      }
    }
  }
}

void runGeometryChecks(const lf::geometry::Geometry &geom,
                       const Eigen::MatrixXd &eval_points,
                       const double tolerance) {
  checkJacobian(geom, eval_points, tolerance);
  // JacobianInverseGramian is not defined for points
  if (geom.RefEl() != lf::base::RefEl::kPoint()) {
    checkJacobianInverseGramian(geom, eval_points);
  }
  checkIntegrationElement(geom, eval_points);
  checkSubGeometry(geom, [](lf::base::RefEl refEl) -> lf::quad::QuadRule {
    return lf::quad::make_QuadRule(refEl, 5);
  });
}

/**
 * Check if the total volume is conserved after call to ChildGeometry()
 */
void checkChildGeometryVolume(
    const lf::geometry::Geometry &geom, const lf::refinement::RefPat &refPat,
    const lf::base::sub_idx_t &anchor = lf::refinement::idx_nil) {
  // compute volume by means of overkill quadrature
  auto computeVolume = [](const lf::geometry::Geometry &geom) {
    const auto qr = lf::quad::make_QuadRule(geom.RefEl(), 23);

    const auto &points = qr.Points();
    const auto &weights = qr.Weights();
    const auto &integrationElements = geom.IntegrationElement(points);

    double vol = 0.;

    for (size_t j = 0; j < points.cols(); ++j) {
      vol += weights(j) * integrationElements(j);
    }

    return vol;
  };

  const double volume = computeVolume(geom);
  double refinedVolume = 0.;
  auto children = geom.ChildGeometry(
      lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat, anchor),
      0);

  for (auto &childGeom : children) {
    refinedVolume += computeVolume(*childGeom);
  }

  switch (refPat) {
    case lf::refinement::RefPat::rp_nil:
      EXPECT_EQ(refinedVolume, 0.)
          << refPat << " should not produce any children";
      break;
    default:
      EXPECT_FLOAT_EQ(volume, refinedVolume)
          << "Parent and children volumes differ for " << refPat;
  }
}

/**
 * @brief Check that the mapping geom.ChildGeometry() is the same as geom
 *    composed with the mapping imposed by ref_pat
 * @param geom The geometry object whose ChildGeometry() method should be
 * checked
 * @param ref_pat The refinement pattern that is used for the test.
 * @param qr_provider Provides the quadrature rules whose points are used to
 *   check whether the two mappings agree.
 *
 * @warning For non-linear mappings, this test can fail, especially for
 *   quadrilaterals which are split into triangles! In this case it makes sense
 *   to use a nodal quadrature rule for the test (see
 *   quad::make_QuadRuleNodal())
 */
void checkChildGeometry(
    const lf::geometry::Geometry &geom,
    const lf::geometry::RefinementPattern &ref_pat,
    const std::function<lf::quad::QuadRule(lf::base::RefEl)> &qr_provider) {
  LF_ASSERT_MSG(geom.RefEl() == ref_pat.RefEl(),
                "This refinement pattern is not made for a " << geom.RefEl());
  double h = 1. / ref_pat.LatticeConst();

  for (auto codim = 0; codim <= geom.RefEl().Dimension(); ++codim) {
    auto children = geom.ChildGeometry(ref_pat, codim);
    auto child_polygons = ref_pat.ChildPolygons(codim);

    EXPECT_EQ(children.size(), ref_pat.noChildren(codim));

    for (int i = 0; i < children.size(); ++i) {
      auto &child = children[i];
      EXPECT_EQ(child->RefEl().Dimension(), geom.DimLocal() - codim);
      EXPECT_EQ(child->DimGlobal(), geom.DimGlobal());
      EXPECT_EQ(child->DimLocal(), geom.DimLocal() - codim);

      if (child->RefEl() == lf::base::RefEl::kPoint()) {
        // check that child == geom.Global(...)
        EXPECT_EQ(child_polygons[i].rows(), geom.DimLocal());
        EXPECT_EQ(child_polygons[i].cols(), 1);
        auto child_coord = geom.Global(child_polygons[i].cast<double>() * h);

        auto zero = Eigen::Matrix<double, 0, 1>::Zero();
        EXPECT_TRUE(child_coord.isApprox(child->Global(zero)));
      } else {
        // Check that  the mapping child is the same as geom \cdot
        // refinementPattern
        auto qr = qr_provider(child->RefEl());
        auto a = child->Global(qr.Points());

        std::unique_ptr<lf::geometry::Geometry> ref_pat_geo = nullptr;

        switch (child->RefEl()) {
          case lf::base::RefEl::kSegment():
            ref_pat_geo = std::make_unique<lf::geometry::SegmentO1>(
                child_polygons[i].cast<double>() * h);
            break;
          case lf::base::RefEl::kTria():
            ref_pat_geo = std::make_unique<lf::geometry::TriaO1>(
                child_polygons[i].cast<double>() * h);
            break;
          case lf::base::RefEl::kQuad():
            ref_pat_geo = std::make_unique<lf::geometry::QuadO1>(
                child_polygons[i].cast<double>() * h);
            break;
          default:
            LF_VERIFY_MSG(
                false,
                "This reference element is not yet supported by the test.");
        }

        auto b = geom.Global(ref_pat_geo->Global(qr.Points()));

        EXPECT_TRUE(a.isApprox(b)) << std::endl
                                   << a << "\n is not \n"
                                   << b << std::endl;
      }
    }
  }
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
      checkChildGeometryVolume(geom, refPat);
      checkChildGeometry(
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
      checkChildGeometryVolume(geom, refPat);
      checkChildGeometry(
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
      checkChildGeometryVolume(geom, refPat);
      checkChildGeometry(
          geom, lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat),
          [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
    }

    std::vector<lf::refinement::RefPat> triaAsymmetricRefPats = {
        lf::refinement::RefPat::rp_bisect, lf::refinement::RefPat::rp_trisect,
        lf::refinement::RefPat::rp_trisect_left,
        lf::refinement::RefPat::rp_quadsect};

    for (const auto &refPat : triaAsymmetricRefPats) {
      for (size_t anchor = 0; anchor < 3; ++anchor) {
        checkChildGeometryVolume(geom, refPat, anchor);
        checkChildGeometry(
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
      checkChildGeometryVolume(geom, refPat);
      checkChildGeometry(
          geom, lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat),
          [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
    }

    std::vector<lf::refinement::RefPat> triaAsymmetricRefPats = {
        lf::refinement::RefPat::rp_bisect, lf::refinement::RefPat::rp_trisect,
        lf::refinement::RefPat::rp_trisect_left,
        lf::refinement::RefPat::rp_quadsect};

    for (const auto &refPat : triaAsymmetricRefPats) {
      for (size_t anchor = 0; anchor < 3; ++anchor) {
        checkChildGeometryVolume(geom, refPat, anchor);
        checkChildGeometry(
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
      checkChildGeometryVolume(geom, refPat);
      checkChildGeometry(
          geom, lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat),
          [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
    }

    std::vector<lf::refinement::RefPat> triaAsymmetricRefPats = {
        lf::refinement::RefPat::rp_split, lf::refinement::RefPat::rp_bisect,
        lf::refinement::RefPat::rp_trisect, lf::refinement::RefPat::rp_quadsect,
        lf::refinement::RefPat::rp_threeedge};

    for (const auto &refPat : triaAsymmetricRefPats) {
      for (size_t anchor = 0; anchor < 4; ++anchor) {
        checkChildGeometryVolume(geom, refPat, anchor);
        checkChildGeometry(geom,
                           lf::refinement::Hybrid2DRefinementPattern(
                               geom.RefEl(), refPat, anchor),
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
      checkChildGeometryVolume(geom, refPat);
      checkChildGeometry(
          geom, lf::refinement::Hybrid2DRefinementPattern(geom.RefEl(), refPat),
          [](auto ref_el) { return lf::quad::make_QuadRule(ref_el, 5); });
    }

    std::vector<lf::refinement::RefPat> triaAsymmetricRefPats = {
        lf::refinement::RefPat::rp_split, lf::refinement::RefPat::rp_bisect,
        lf::refinement::RefPat::rp_trisect, lf::refinement::RefPat::rp_quadsect,
        lf::refinement::RefPat::rp_threeedge};

    for (const auto &refPat : triaAsymmetricRefPats) {
      for (size_t anchor = 0; anchor < 4; ++anchor) {
        checkChildGeometryVolume(geom, refPat, anchor);
        checkChildGeometry(geom,
                           lf::refinement::Hybrid2DRefinementPattern(
                               geom.RefEl(), refPat, anchor),
                           lf::quad::make_QuadRuleNodal);
      }
    }

    // check that local nodes are mapped to global nodes
    Eigen::MatrixXd local_nodes(2, 8);
    local_nodes << 0, 1, 1, 0, .5, 1, .5, 0, 0, 0, 1, 1, 0, .5, 1, 0.5;
    EXPECT_TRUE(geom.Global(local_nodes).isApprox(global_nodes));
  }
}
