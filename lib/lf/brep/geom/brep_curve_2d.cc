/**
 * @file
 * @brief Implementation of BrepCurve2d
 * @author Raffael Casagrande
 * @date   2021-02-10 04:55:28
 * @copyright MIT License
 */

#include "brep_curve_2d.h"

namespace lf::brep::geom {
Eigen::MatrixXd BrepCurve2d::Global(const Eigen::MatrixXd& local) const {
  auto result = curve_->Global(local);
  LF_ASSERT_MSG(result.row(2).isZero(), "Expected the z-component to be zero.");
  return result.topRows(2);
}

Eigen::MatrixXd BrepCurve2d::Jacobian(const Eigen::MatrixXd& local) const {
  auto result = curve_->Jacobian(local);
  LF_ASSERT_MSG(result.row(2).isZero(), "Expected the z-component to be zero.");
  return result.topRows(2);
}

std::vector<bool> BrepCurve2d::IsInBoundingBox(
    const Eigen::MatrixXd& global) const {
  LF_ASSERT_MSG(global.rows() == 2, "expected 2d global.");
  Eigen::MatrixXd global3d(3, global.cols());
  global3d.topRows(2) = global;
  global3d.row(2).setZero();
  return curve_->IsInBoundingBox(global3d);
}

std::pair<double, Eigen::VectorXd> BrepCurve2d::Project(
    const Eigen::VectorXd& global) const {
  LF_ASSERT_MSG(global.rows() == 2, "expected 2d global.");
  Eigen::MatrixXd global3d(3, global.cols());
  global3d.topRows(2) = global;
  global3d.row(2).setZero();
  return curve_->Project(global3d);
}
}  // namespace lf::brep::geom
