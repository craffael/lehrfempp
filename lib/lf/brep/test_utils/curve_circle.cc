/**
 * @file
 * @brief Implementation of CurveCircle2D
 * @author Raffael Casagrande
 * @date   2021-01-29 08:34:53
 * @copyright MIT License
 */

#include "curve_circle.h"

namespace lf::brep::test_utils {
CurveCircle::CurveCircle(const Eigen::Vector3d& origin, double radius)
    : origin_(origin), radius_(radius) {
  LF_ASSERT_MSG(radius_ > 0, "Radius must be positive.");
}

Eigen::MatrixXd CurveCircle::GlobalMulti(const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(local.rows() == 1, "expected exactly one row.");
  Eigen::MatrixXd result(3, local.cols());
  result.row(0).array() = local.array().cos() * radius_ + origin_.x();
  result.row(1).array() = local.array().sin() * radius_ + origin_.y();
  result.row(2).setConstant(origin_.z());
  return result;
}

Eigen::MatrixXd CurveCircle::JacobianMulti(const Eigen::MatrixXd& local) const {
  Eigen::MatrixXd result(3, local.cols());
  result.row(0).array() = -local.array().sin() * radius_;
  result.row(1).array() = local.array().cos() * radius_;
  result.row(2).setZero();
  return result;
}

std::vector<bool> CurveCircle::IsInBoundingBoxMulti(
    const Eigen::MatrixXd& global) const {
  std::vector<bool> result(global.cols());
  constexpr double eps = 1e-7;
  for (int i = 0; i < global.cols(); ++i) {
    result[i] = global(0, i) >= origin_.x() - radius_ - eps &&
                global(0, i) <= origin_.x() + radius_ + eps &&
                global(1, i) >= origin_.y() - radius_ - eps &&
                global(1, i) <= origin_.y() + radius_ + eps &&
                std::abs(global(2, i) - origin_.z()) <= eps;
  }
  return result;
}

Eigen::Vector3d CurveCircle::GlobalSingle(double local) const {
  return GlobalMulti((Eigen::Matrix<double, 1, 1>() << local).finished())
      .col(0);
}

Eigen::Vector3d CurveCircle::JacobianSingle(double local) const {
  return JacobianMulti((Eigen::Matrix<double, 1, 1>() << local).finished())
      .col(0);
}

std::pair<double, double> CurveCircle::Project(
    const Eigen::Vector3d& global) const {
  auto temp = (global - origin_).topRows(2).eval();
  Eigen::Vector3d globalp;
  globalp.topRows(2) = temp.normalized() * radius_ + origin_.topRows(2);
  globalp(2) = origin_.z();

  return {(global - globalp).norm(), std::atan2(temp.y(), temp.x())};
}

bool CurveCircle::IsInBoundingBoxSingle(const Eigen::Vector3d& global) const {
  return IsInBoundingBoxMulti(global)[0];
}

bool CurveCircle::IsInside(double local) const {
  constexpr double eps = 1e-7;
  return local >= -base::kPi - eps && local <= base::kPi + eps;
}

}  // namespace lf::brep::test_utils
