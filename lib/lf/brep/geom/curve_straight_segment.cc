/**
 * @file
 * @brief implementation of CurveStraightSegment...
 * @author Raffael Casagrande
 * @date   2021-02-08 03:57:07
 * @copyright MIT License
 */

#include "curve_straight_segment.h"

namespace lf::brep::geom {
template <int DIM_WORLD>
Eigen::MatrixXd CurveStraightLine<DIM_WORLD>::Global(
    const Eigen::MatrixXd& local) const {
  return transform_ * local + offset_.replicate(1, local.cols());
}

template <int DIM_WORLD>
Eigen::MatrixXd CurveStraightLine<DIM_WORLD>::Jacobian(
    const Eigen::MatrixXd& local) const {
  return transform_.replicate(1, local.cols());
}

template <int DIM_WORLD>
std::vector<bool> CurveStraightLine<DIM_WORLD>::IsInBoundingBox(
    const Eigen::MatrixXd& global) const {
  LF_ASSERT_MSG(global.rows() == DIM_WORLD, "global has wrong number of rows.");
  std::vector<bool> result(global.cols());
  Eigen::Matrix<double, DIM_WORLD, 1> min, max;
  min = offset_.cwiseMin(offset_ + transform_);
  max = offset_.cwiseMax(offset_ + transform_);
  constexpr double eps = 1e-6;
  for (int i = 0; i < global.cols(); ++i) {
    result[i] = global(0, i) >= min(0) - eps && global(0, i) <= max(0) + eps &&
                global(1, i) >= min(1) - eps && global(1, i) <= max(1) + eps;
    if constexpr (DIM_WORLD == 3) {
      result[i] = result[i] && global(2, i) >= min(2) - eps &&
                  global(2, i) <= max(2) + eps;
    }
  }
  return result;
}

template <int DIM_WORLD>
std::pair<double, Eigen::VectorXd> CurveStraightLine<DIM_WORLD>::Project(
    const Eigen::VectorXd& global) const {
  LF_ASSERT_MSG(global.rows() == DIM_WORLD, "Global has wrong number of rows.");
  std::pair<double, Eigen::VectorXd> result;
  result.second = Eigen::VectorXd(1);
  result.second(0) =
      (global - offset_).dot(transform_) / transform_.squaredNorm();
  result.first = (Global(result.second) - global).norm();
  return result;
}

template <int DIM_WORLD>
bool CurveStraightLine<DIM_WORLD>::IsInside(
    const Eigen::VectorXd& local) const {
  LF_ASSERT_MSG(local.rows() == 1, "Local should have 1 row.");
  constexpr double eps = 1e-6;
  return local(0) >= -eps && local(0) <= 1 + eps;
}

// explicit template instantiation defintion.
template class CurveStraightLine<2>;
template class CurveStraightLine<3>;
}  // namespace lf::brep::geom
