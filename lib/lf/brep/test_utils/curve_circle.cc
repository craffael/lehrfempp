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

Eigen::MatrixXd CurveCircle::Global(const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(local.rows() == 1, "expected exactly one row.");
  Eigen::MatrixXd result(3, local.cols());
  result.row(0).array() = local.array().cos() * radius_ + origin_.x();
  result.row(1).array() = local.array().sin() * radius_ + origin_.y();
  result.row(2).setConstant(origin_.z());
  return result;
}

Eigen::MatrixXd CurveCircle::Jacobian(const Eigen::MatrixXd& local) const {
  Eigen::MatrixXd result(3, local.cols());
  result.row(0).array() = -local.array().sin() * radius_;
  result.row(1).array() = local.array().cos() * radius_;
  result.row(2).setZero();
  return result;
}

std::vector<bool> CurveCircle::IsInBoundingBox(
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

std::pair<double, Eigen::VectorXd> CurveCircle::Project(
    const Eigen::VectorXd& global) const {
  LF_ASSERT_MSG(global.rows() == 3, "expected global to have 3 rows.");
  auto temp = (global - origin_).topRows(2).eval();
  Eigen::Vector3d globalp;
  if (temp.squaredNorm() > 0) {
    globalp.topRows(2) = temp.normalized() * radius_ + origin_.topRows(2);
  } else {
    globalp.topRows(2) =
        Eigen::Vector2d::UnitX() * radius_ + origin_.topRows(2);
  }
  globalp(2) = origin_.z();

  return {(global - globalp).norm(),
          (Eigen::VectorXd(1) << std::atan2(temp.y(), temp.x())).finished()};
}

bool CurveCircle::IsInside(const Eigen::VectorXd& local) const {
  LF_ASSERT_MSG(local.rows() == 1, "local should have exactly 1 row.");
  constexpr double eps = 1e-7;
  return local(0) >= -base::kPi - eps && local(0) <= base::kPi + eps;
}

}  // namespace lf::brep::test_utils
