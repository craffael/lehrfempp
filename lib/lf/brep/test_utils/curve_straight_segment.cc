/**
 * @file
 * @brief implementation of CurveStraightSegment...
 * @author Raffael Casagrande
 * @date   2021-02-08 03:57:07
 * @copyright MIT License
 */

#include "curve_straight_segment.h"

namespace lf::brep::test_utils {
Eigen::MatrixXd CurveStraightLine::GlobalMulti(
    const Eigen::MatrixXd& local) const {
  return transform_ * local + offset_.replicate(1, local.cols());
}

Eigen::MatrixXd CurveStraightLine::JacobianMulti(
    const Eigen::MatrixXd& local) const {
  return transform_.replicate(1, local.cols());
}

std::vector<bool> CurveStraightLine::IsInBoundingBoxMulti(
    const Eigen::MatrixXd& global) const {
  LF_ASSERT_MSG(global.rows() == 3, "global must be 3d-coordinates.");
  std::vector<bool> result(global.cols());
  Eigen::Vector3d min, max;
  min = offset_.cwiseMin(offset_ + transform_);
  max = offset_.cwiseMax(offset_ + transform_);
  constexpr double eps = 1e-6;
  for (int i = 0; i < global.cols(); ++i) {
    result[i] = global(0, i) >= min(0) - eps && global(0, i) <= max(0) + eps &&
                global(1, i) >= min(1) - eps && global(1, i) <= max(1) + eps &&
                global(2, i) >= min(2) - eps && global(2, i) <= max(2) + eps;
  }
  return result;
}

std::pair<double, double> CurveStraightLine::Project(
    const Eigen::Vector3d& global) const {
  std::pair<double, double> result;
  result.second = (global - offset_).dot(transform_) / transform_.squaredNorm();
  result.first = (GlobalSingle(result.second) - global).norm();
  return result;
}

bool CurveStraightLine::IsInBoundingBoxSingle(
    const Eigen::Vector3d& global) const {
  Eigen::Vector3d min, max;
  min = offset_.cwiseMin(offset_ + transform_);
  max = offset_.cwiseMax(offset_ + transform_);
  constexpr double eps = 1e-6;
  return global(0, 0) >= min(0) - eps && global(0, 0) <= max(0) + eps &&
         global(1, 0) >= min(1) - eps && global(1, 0) <= max(1) + eps &&
         global(2, 0) >= min(2) - eps && global(2, 0) <= max(2) + eps;
}

bool CurveStraightLine::IsInside(double local) const {
  constexpr double eps = 1e-6;
  return local >= -eps && local <= 1 + eps;
}
}  // namespace lf::brep::test_utils
