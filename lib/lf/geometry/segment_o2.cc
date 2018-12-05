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

SegmentO2::SegmentO2(Eigen::Matrix<double, Eigen::Dynamic, 3> coords)
    : coords_(std::move(coords)),
      alpha_(coords_.rows()),
      beta_(coords_.rows()),
      gamma_(coords_.rows()),
      alpha_squared_(0),
      alpha_beta_(0),
      beta_squared_(0) {
  const Eigen::VectorXd& vtx0 = coords_.col(0);
  const Eigen::VectorXd& vtx1 = coords_.col(1);
  const Eigen::VectorXd& midp = coords_.col(2);

  // polynomial of degree 2: alpha * x^2 + beta * x + gamma
  alpha_ = 2. * (vtx1 + vtx0) - 4. * midp;
  beta_ = 4. * midp - 3. * vtx0 - vtx1;
  gamma_ = vtx0;

  // coefficients for JacobianInverseGramian and IntegrationElement
  alpha_squared_ = alpha_.squaredNorm();
  alpha_beta_ = alpha_.dot(beta_);
  beta_squared_ = beta_.squaredNorm();
}

Eigen::MatrixXd SegmentO2::Global(const Eigen::MatrixXd& local) const {
  Eigen::VectorXd local_vec = local.transpose();

  if ((0. <= local.array()).all() && (local.array() <= 1.).all()) {
    // evaluate polynomial with Horner scheme
    auto tmp = ((alpha_ * local).colwise() + beta_).array().rowwise();
    return (tmp * local_vec.transpose().array()).matrix().colwise() + gamma_;
  } else {
    LF_VERIFY_MSG(false,
                  "local coordinates out of bounds for reference element");
  }
}

Eigen::MatrixXd SegmentO2::Jacobian(const Eigen::MatrixXd& local) const {
  if ((0. <= local.array()).all() && (local.array() <= 1.).all()) {
    return (2. * alpha_ * local).colwise() + beta_;
  } else {
    LF_VERIFY_MSG(false,
                  "local coordinates out of bounds for reference element");
  }
}

Eigen::MatrixXd SegmentO2::JacobianInverseGramian(
    const ::Eigen::MatrixXd& local) const {
  if ((0. <= local.array()).all() && (local.array() <= 1.).all()) {
    auto jacobian = (2. * alpha_ * local).colwise() + beta_;

    if (DimGlobal() == 1) {
      return jacobian.cwiseInverse();
    } else {
      // evaluate polynomial with Horner scheme
      const auto jTj =
          4. * local.array() * (local.array() * alpha_squared_ + alpha_beta_) +
          beta_squared_;
      const Eigen::VectorXd jTj_inv = jTj.cwiseInverse().transpose();
      return jacobian.array().rowwise() * jTj_inv.transpose().array();
    }
  } else {
    LF_VERIFY_MSG(false,
                  "local coordinates out of bounds for reference element");
  }
}

Eigen::VectorXd SegmentO2::IntegrationElement(
    const Eigen::MatrixXd& local) const {
  if ((0. <= local.array()).all() && (local.array() <= 1.).all()) {
    const auto jTj =
        4. * local.array() * (local.array() * alpha_squared_ + alpha_beta_) +
        beta_squared_;
    return jTj.cwiseSqrt().transpose();
  } else {
    LF_VERIFY_MSG(false,
                  "local coordinates out of bounds for reference element");
  }
}  // namespace lf::geometry

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
      noChildren == ref_pat.noChildren(codim),
      "noChildren " << noChildren << " <-> " << ref_pat.noChildren(codim));

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
