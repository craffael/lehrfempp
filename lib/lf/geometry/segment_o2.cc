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

  for (size_t point = 0; point < local.cols(); ++point) {
    const double x = local(point);

    if (0. <= x and x <= 0.5) {
      global.col(point) = 2. * (midp - vtx0) * x + vtx0;
    } else if (0.5 < x and x <= 1.) {
      global.col(point) = 2. * (vtx1 - midp) * x + (2. * midp - vtx1);
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

  for (size_t point = 0; point < local.cols(); ++point) {
    const double x = local(point);

    if (0. <= x and x < 0.5) {
      jacobian.col(point) = 2. * (midp - vtx0);
    } else if (x == 0.5) {
      jacobian.col(point) = (vtx1 - vtx0);
    } else if (0.5 < x and x <= 1.) {
      jacobian.col(point) = 2. * (vtx1 - midp);
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
    LF_ASSERT_MSG(0 <= i and i <= 2, "i is out of bounds");
    return std::make_unique<Point>(coords_.col(i));
  }
  LF_VERIFY_MSG(false, "codim is out of bounds");
}

std::vector<std::unique_ptr<Geometry>> SegmentO2::ChildGeometry(
    const RefinementPattern& ref_pat,  // NOLINT(misc-unused-parameters)
    base::dim_t codim) const {         // NOLINT(misc-unused-parameters)
  LF_VERIFY_MSG(false, "refinement for SegmentO2 not implemented yet");
}

}  // namespace lf::geometry
