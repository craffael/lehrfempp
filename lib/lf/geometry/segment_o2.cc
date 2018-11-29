/**
 * @file
 * @brief Implementation of second-order segments
 * @author Anian Ruoss
 * @date   2018-11-18 19:02:17
 * @copyright MIT License
 */

#include "segment_o2.h"
#include "point.h"

namespace lf::geometry {

Eigen::MatrixXd SegmentO2::Global(const Eigen::MatrixXd& local) const {
  Eigen::MatrixXd global(DimGlobal(), local.cols());

  const Eigen::VectorXd vtx0 = coords_.col(0);
  const Eigen::VectorXd vtx1 = coords_.col(1);
  const Eigen::VectorXd midp = coords_.col(2);

  // polynomial of degree 2: alpha * x^2 + beta * x + gamma
  const Eigen::VectorXd alpha = 2. * (vtx1 + vtx0) - 4. * midp;
  const Eigen::VectorXd beta = 4. * midp - 3. * vtx0 - vtx1;
  const Eigen::VectorXd& gamma = vtx0;

  for (size_t point = 0; point < local.cols(); ++point) {
    const double x = local(point);

    if (0. <= x && x <= 1.) {
      global.col(point) = alpha * x * x + beta * x + gamma;
    } else {
      LF_VERIFY_MSG(false,
                    "local coordinate out of bounds for reference element");
    }
  }

  return global;
}

Eigen::MatrixXd SegmentO2::Jacobian(const Eigen::MatrixXd& local) const {
  Eigen::MatrixXd jacobian(DimGlobal(), DimLocal() * local.cols());

  const Eigen::VectorXd vtx0 = coords_.col(0);
  const Eigen::VectorXd vtx1 = coords_.col(1);
  const Eigen::VectorXd midp = coords_.col(2);

  // polynomial of degree 2: alpha * x^2 + beta * x + gamma
  const Eigen::VectorXd alpha = 2. * (vtx1 + vtx0) - 4. * midp;
  const Eigen::VectorXd beta = 4. * midp - 3. * vtx0 - vtx1;

  for (size_t point = 0; point < local.cols(); ++point) {
    const double x = local(point);

    if (0. <= x && x <= 1.) {
      jacobian.col(point) = 2. * alpha * x + beta;
    } else {
      LF_VERIFY_MSG(false,
                    "local coordinate out of bounds for reference element");
    }
  }

  return jacobian;
}

Eigen::MatrixXd SegmentO2::JacobianInverseGramian(
    const ::Eigen::MatrixXd& local) const {
  Eigen::MatrixXd jacInvGram(DimGlobal(), DimLocal() * local.cols());
  Eigen::MatrixXd jacobian = Jacobian(local);

  for (size_t point = 0; point < local.cols(); ++point) {
    Eigen::MatrixXd jac =
        jacobian.block(0, point * DimLocal(), DimGlobal(), DimLocal());

    if (DimGlobal() == DimLocal()) {
      jacInvGram.block(0, point * DimLocal(), DimGlobal(), DimLocal()) =
          jac.inverse().transpose();
    } else {
      jacInvGram.block(0, point * DimLocal(), DimGlobal(), DimLocal()) =
          jac * (jac.transpose() * jac).inverse();
    }
  }

  return jacInvGram;
}

Eigen::VectorXd SegmentO2::IntegrationElement(
    const Eigen::MatrixXd& local) const {
  Eigen::VectorXd intEl(local.cols());
  Eigen::MatrixXd jacobian = Jacobian(local);

  for (size_t point = 0; point < local.cols(); ++point) {
    Eigen::MatrixXd jac =
        jacobian.block(0, point * DimLocal(), DimGlobal(), DimLocal());

    if (DimGlobal() == DimLocal()) {
      intEl(point) = std::abs(jac.determinant());
    } else {
      intEl(point) = std::sqrt((jac.transpose() * jac).determinant());
    }
  }

  return intEl;
}

std::unique_ptr<Geometry> SegmentO2::SubGeometry(dim_t codim, dim_t i) const {
  if (codim == 0) {
    LF_ASSERT_MSG(i == 0, "i is out of bounds");
    return std::make_unique<SegmentO2>(coords_);
  }
  if (codim == 1) {
    LF_ASSERT_MSG(0 <= i && i <= 2, "i is out of bounds");
    return std::make_unique<Point>(coords_.col(i));
  }
  LF_VERIFY_MSG(false, "codim is out of bounds");
}

std::vector<std::unique_ptr<Geometry>> SegmentO2::ChildGeometry(
    const RefinementPattern& refPat, base::dim_t codim) const {
  LF_VERIFY_MSG(refPat.RefEl() == lf::base::RefEl::kSegment(),
                "Refinement pattern for " << refPat.RefEl().ToString());
  LF_VERIFY_MSG(codim < 2, "Illegal codim = " << codim);

  const double hLattice = 1. / static_cast<double>(refPat.LatticeConst());
  std::vector<std::unique_ptr<Geometry>> childGeoPtrs = {};

  // get coordinates of childGeometries
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>> childPolygons(
      refPat.ChildPolygons(codim));

  const size_t noChildren = childPolygons.size();
  LF_VERIFY_MSG(
      noChildren == refPat.noChildren(codim),
      "noChildren " << noChildren << " <-> " << refPat.noChildren(codim));

  // create a geometry object for each child
  for (size_t child = 0; child < noChildren; ++child) {
    // codim == 0: A single child must be described by an interval, that is
    // two different lattice coordinates
    // codim == 1: A point's location is just one number
    LF_VERIFY_MSG(
        childPolygons[child].rows() == 1,
        "childPolygons[child].rows() = " << childPolygons[child].rows());
    LF_VERIFY_MSG(
        childPolygons[child].cols() == (2 - codim),
        "childPolygons[child].cols() = " << childPolygons[child].rows());

    Eigen::MatrixXd locCoords(1, 3 - 2 * codim);

    switch (codim) {
      case 0: {
        // SegmentO2 requires the coordinate of the midpoint
        locCoords << hLattice * childPolygons[child].cast<double>(),
            (hLattice * childPolygons[child].cast<double>()).sum() / 2;
        childGeoPtrs.push_back(std::make_unique<SegmentO2>(Global(locCoords)));

        break;
      }
      case 1: {
        locCoords << hLattice * childPolygons[child].cast<double>();
        childGeoPtrs.push_back(std::make_unique<Point>(Global(locCoords)));

        break;
      }
      default: { LF_VERIFY_MSG(false, "Unreachable code"); }
    }
  }

  return childGeoPtrs;
}

}  // namespace lf::geometry
