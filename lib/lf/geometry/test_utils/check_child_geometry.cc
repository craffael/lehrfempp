/**
 * @file
 * @brief Implementation of ChildGeometry() tests for geometry objects
 * @author Anian Ruoss
 * @date   2019-02-11 18:08:17
 * @copyright MIT License
 */

#include "check_child_geometry.h"

#include <gtest/gtest.h>

namespace lf::geometry::test_utils {

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

    EXPECT_EQ(children.size(), ref_pat.NumChildren(codim));

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

void checkChildGeometryVolume(const lf::geometry::Geometry &geom,
                              const lf::refinement::RefPat &refPat,
                              const lf::base::sub_idx_t &anchor) {
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

}  // namespace lf::geometry::test_utils
