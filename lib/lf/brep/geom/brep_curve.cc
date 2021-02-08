/**
 * @file
 * @brief Implementation of BrepSegment
 * @author Raffael Casagrande
 * @date   2021-02-01 03:59:37
 * @copyright MIT License
 */

#include "brep_curve.h"

namespace lf::brep::geom {
BrepCurve::BrepCurve(const interface::BrepCurve* curve,
                     Eigen::RowVector2d curve_param, bool delete_curve)
    : curve_(curve),
      offset_(curve_param(0)),
      slope_(curve_param(1) - curve_param(0)),
      delete_curve_(delete_curve) {}

Eigen::MatrixXd BrepCurve::Global(const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(local.rows() == 1, "unexpected # of rows.");
  return curve_->GlobalMulti((slope_ * local.array() + offset_).matrix());
}

Eigen::MatrixXd BrepCurve::Jacobian(const Eigen::MatrixXd& local) const {
  return curve_->JacobianMulti((slope_ * local.array() + offset_).matrix()) *
         slope_;
}

Eigen::MatrixXd BrepCurve::JacobianInverseGramian(
    const Eigen::MatrixXd& local) const {
  auto jac = Jacobian(local);
  return jac * (jac.colwise().squaredNorm().cwiseInverse().asDiagonal());
}

Eigen::VectorXd BrepCurve::IntegrationElement(
    const Eigen::MatrixXd& local) const {
  auto jac = Jacobian(local);
  return jac.colwise().norm().transpose();
}

std::unique_ptr<geometry::Geometry> BrepCurve::SubGeometry(dim_t codim,
                                                           dim_t i) const {
  LF_ASSERT_MSG(codim >= 0 && codim <= 1, "codim out of bounds.");
  if (codim == 0) {
    LF_ASSERT_MSG(i == 0, "subidx out of bounds.");
    return std::make_unique<BrepCurve>(*this);
  }
  LF_ASSERT_MSG(i >= 0 && i <= 1, "subidx out of bounds.");
  if (i == 0) {
    return std::make_unique<geometry::Point>(curve_->GlobalSingle(offset_));
  }
  return std::make_unique<geometry::Point>(
      curve_->GlobalSingle(offset_ + slope_));
}

std::vector<std::unique_ptr<geometry::Geometry>> BrepCurve::ChildGeometry(
    const geometry::RefinementPattern& ref_pat, lf::base::dim_t codim) const {
  LF_ASSERT_MSG(codim >= 0 && codim <= 1, "codim out of bounds.");
  std::vector<std::unique_ptr<Geometry>> result;
  result.reserve(ref_pat.NumChildren(codim));
  double lattic_const = static_cast<double>(ref_pat.LatticeConst());
  if (codim == 1) {
    for (auto& n : ref_pat.ChildPolygons(1)) {
      LF_ASSERT_MSG(n.cols() == 1, "Unexpected # of columns.");
      LF_ASSERT_MSG(n.rows() == 1, "Unexpected # of rows.");
      result.push_back(std::make_unique<geometry::Point>(
          Global(n.cast<double>() / lattic_const)));
    }
  } else {
    // codim == 0:
    for (auto& e : ref_pat.ChildPolygons(0)) {
      LF_ASSERT_MSG(e.cols() == 2, "Unexpected # of columns.");
      LF_ASSERT_MSG(e.rows() == 1, "Unexpected # of rows.");
      auto local = (e.cast<double>().array() / lattic_const * slope_ + offset_)
                       .matrix()
                       .eval();
      result.push_back(
          std::make_unique<BrepCurve>(curve_, local, delete_curve_));
    }
  }
  return result;
}

BrepCurve::~BrepCurve() {
  if (delete_curve_) delete curve_;
}
}  // namespace lf::brep::geom
