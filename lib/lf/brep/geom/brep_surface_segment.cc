/**
 * @file
 * @brief implementation of BrepSurfaceCurve
 * @author Raffael Casagrande
 * @date   2021-01-29 04:13:42
 * @copyright MIT License
 */

#include "brep_surface_segment.h"

namespace lf::brep::geom {
Eigen::MatrixXd BrepSurfaceSegment::GlobalMulti(
    const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(local.rows() == 1, "local can have at most 1 row.");
  return surface_->GlobalMulti(local_direction_ * local + local_offset_);
}

Eigen::MatrixXd BrepSurfaceSegment::JacobianMulti(
    const Eigen::MatrixXd& local) const {
  auto jac = surface_->JacobianMulti(local_direction_ * local + local_offset_);
  Eigen::MatrixXd result(3, local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    result.col(i) = jac.block(0, 2 * i, 3, 2) * local_direction_;
  }
  return result;
}

std::pair<double, double> BrepSurfaceSegment::Project(
    const Eigen::Vector3d& global) const {
  // Use Gauss-Newtons method to find the local minima. Start with 10 initial,
  // linearly spaced guesses:
  constexpr int num_guesses = 10;
  Eigen::RowVectorXd x = Eigen::RowVectorXd::LinSpaced(num_guesses, 0, 1.);
  Eigen::Matrix<bool, 1, Eigen::Dynamic> converged(num_guesses);
  converged.setConstant(false);

  for (int i = 0; i < 100; ++i) {
    // at most 100 steps
    auto jac = JacobianMulti(x);
    auto delta = (jac.colwise().squaredNorm().inverse() * jac.transpose() *
                  (global - GlobalMulti(x)))
                     .eval();
    x = x + converged.select(Eigen::RowVectorXd::Zero(num_guesses), delta);
    converged = converged.cwiseMax((delta.array() < 1e-6).matrix());
    if (converged.all()) {
      break;
    }
  }

  std::pair<double, double> result{1e6, 0};
  auto distances = (GlobalMulti(x) - global).colwise().norm().eval();
  for (int i = 0; i < num_guesses; ++i) {
    if (distances(i) < result.first) {
      result.first = distances(i);
      result.second = x(i);
    }
  }
  LF_ASSERT_MSG(result.first < 1e5, "Projection failed.");
  return result;
}

bool BrepSurfaceSegment::IsInBoundingBoxSingle(
    const Eigen::Vector3d& global) const {
  return surface_->IsInBoundingBox(global);
}

bool BrepSurfaceSegment::IsInside(double local) const {
  constexpr double eps = 1e-5;
  return local >= -eps && local <= 1 + eps;
}
}  // namespace lf::brep::geom
