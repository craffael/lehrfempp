/**
 * @file
 * @brief Implementation of SubGeometry() tests for geometry objects
 * @author Anian Ruoss
 * @date   2019-02-11 18:06:17
 * @copyright MIT License
 */

#include "check_sub_geometry.h"

#include <gtest/gtest.h>

namespace lf::geometry::test_utils {

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

}  // namespace lf::geometry::test_utils
