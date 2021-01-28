/**
 * @file
 * @brief Implementation of second-order parametric segments
 * @author Anian Ruoss
 * @date   2018-11-18 19:02:17
 * @copyright MIT License
 */

#include "segment_o2.h"

#include "point.h"

namespace lf::geometry {

SegmentO2::SegmentO2(Eigen::Matrix<double, Eigen::Dynamic, 3> coords)
    : coords_(std::move(coords)),
      alpha_squared_(0),
      alpha_beta_(0),
      beta_squared_(0) {
  // polynomial of degree 2: alpha * x^2 + beta * x + gamma
  alpha_ = 2. * (coords_.col(1) + coords_.col(0)) - 4. * coords_.col(2);
  beta_ = 4. * coords_.col(2) - 3. * coords_.col(0) - coords_.col(1);
  gamma_ = coords_.col(0);

  // coefficients for JacobianInverseGramian and IntegrationElement
  alpha_squared_ = alpha_.squaredNorm();
  alpha_beta_ = alpha_.dot(beta_);
  beta_squared_ = beta_.squaredNorm();
}

Eigen::MatrixXd SegmentO2::Global(const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  // evaluate polynomial with Horner scheme
  Eigen::VectorXd local_vec = local.transpose();
  auto tmp = ((alpha_ * local).colwise() + beta_).array().rowwise();

  return (tmp * local_vec.transpose().array()).matrix().colwise() + gamma_;
}

Eigen::MatrixXd SegmentO2::Jacobian(const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  return (2. * alpha_ * local).colwise() + beta_;
}

Eigen::MatrixXd SegmentO2::JacobianInverseGramian(
    const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  auto jacobian = (2. * alpha_ * local).colwise() + beta_;

  if (DimGlobal() == 1) {
    return jacobian.cwiseInverse();
  }

  // evaluate polynomial with Horner scheme
  const auto jTj =
      4. * local.array() * (local.array() * alpha_squared_ + alpha_beta_) +
      beta_squared_;
  const Eigen::VectorXd jTj_inv = jTj.cwiseInverse().transpose();

  return jacobian.array().rowwise() * jTj_inv.transpose().array();
}

Eigen::VectorXd SegmentO2::IntegrationElement(
    const Eigen::MatrixXd& local) const {
  LF_VERIFY_MSG((0. <= local.array()).all() && (local.array() <= 1.).all(),
                "local coordinates out of bounds for reference element");

  const auto jTj =
      4. * local.array() * (local.array() * alpha_squared_ + alpha_beta_) +
      beta_squared_;

  return jTj.cwiseSqrt().transpose();
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
    const RefinementPattern& ref_pat, base::dim_t codim) const {
  LF_VERIFY_MSG(ref_pat.RefEl() == lf::base::RefEl::kSegment(),
                "Refinement pattern for " << ref_pat.RefEl().ToString());
  LF_VERIFY_MSG(codim < 2, "Illegal codim = " << codim);

  const double hLattice = 1. / static_cast<double>(ref_pat.LatticeConst());
  std::vector<std::unique_ptr<Geometry>> childGeoPtrs = {};

  // get coordinates of childGeometries
  std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>> childPolygons(
      ref_pat.ChildPolygons(codim));

  const size_t noChildren = childPolygons.size();
  LF_VERIFY_MSG(
      noChildren == ref_pat.NumChildren(codim),
      "NumChildren " << noChildren << " <-> " << ref_pat.NumChildren(codim));

  // create a geometry object for each child
  for (size_t child = 0; child < noChildren; ++child) {
    // codim == 0: a child segment is described by a polygon with three vertices
    // codim == 1: a child point by a single point ("polygon with one corner")
    LF_VERIFY_MSG(
        childPolygons[child].rows() == 1,
        "childPolygons[child].rows() = " << childPolygons[child].rows());
    LF_VERIFY_MSG(
        childPolygons[child].cols() == (2 - codim),
        "childPolygons[child].cols() = " << childPolygons[child].cols());

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
      default: {
        LF_VERIFY_MSG(false, "Illegal co-dimension");
      }
    }
  }

  return childGeoPtrs;
}

}  // namespace lf::geometry
