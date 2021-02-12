/**
 * @file
 * @brief implementation of BrepSurfaceSegment
 * @author Raffael Casagrande
 * @date   2021-01-29 04:13:42
 * @copyright MIT License
 */

#include "brep_surface_segment.h"

#include <Eigen/Eigen>
namespace lf::brep::geom {
Eigen::MatrixXd BrepSurfaceSegment::Global(const Eigen::MatrixXd& local) const {
  LF_ASSERT_MSG(local.rows() == 1, "local can have at most 1 row.");
  return surface_->Global(local_direction_ * local +
                          local_offset_.replicate(1, local.cols()));
}

Eigen::MatrixXd BrepSurfaceSegment::Jacobian(
    const Eigen::MatrixXd& local) const {
  auto jac = surface_->Jacobian(local_direction_ * local +
                                local_offset_.replicate(1, local.cols()));
  Eigen::MatrixXd result(3, local.cols());
  for (int i = 0; i < local.cols(); ++i) {
    result.col(i) = jac.block(0, 2 * i, 3, 2) * local_direction_;
  }
  return result;
}

std::pair<double, Eigen::VectorXd> BrepSurfaceSegment::Project(
    const Eigen::VectorXd& global) const {
  LF_ASSERT_MSG(global.rows() == DimGlobal(),
                "illegal numer of rows in parameter 'global'");
  // Use Gauss-Newtons method to find the local minima. Start with 10 initial,
  // linearly spaced guesses:
  constexpr int num_guesses = 10;
  Eigen::RowVectorXd x = Eigen::RowVectorXd::LinSpaced(num_guesses, 0, 1.);
  Eigen::Matrix<bool, 1, Eigen::Dynamic> converged(num_guesses);
  converged.setConstant(false);

  for (int i = 0; i < 100; ++i) {
    // at most 100 steps
    auto jac = Jacobian(x);
    auto global_temp = Global(x);
    Eigen::RowVectorXd delta(num_guesses);
    for (int j = 0; j < num_guesses; ++j) {
      delta(j) = jac.col(j).dot(global - global_temp.col(j)) /
                 jac.col(j).squaredNorm();
    }

    x = x + converged.select(Eigen::RowVectorXd::Zero(num_guesses), delta);
    converged = converged.cwiseMax((delta.array() < 1e-6).matrix());
    if (converged.all()) {
      break;
    }
  }

  std::pair<double, Eigen::VectorXd> result{1e6, Eigen::VectorXd(1)};
  result.second(0) = 0.;
  auto distances =
      (Global(x) - global.replicate(1, num_guesses)).colwise().norm().eval();
  for (int i = 0; i < num_guesses; ++i) {
    if (distances(i) < result.first) {
      result.first = distances(i);
      result.second(0) = x(i);
    }
  }
  LF_ASSERT_MSG(result.first < 1e5, "Projection failed.");
  return result;
}

bool BrepSurfaceSegment::IsInside(const Eigen::VectorXd& local) const {
  constexpr double eps = 1e-5;
  return local(0) >= -eps && local(0) <= 1 + eps;
}

Eigen::VectorXd BrepSurfaceSegment::Periods() const {
  auto surface_period = surface_->Periods();
  Eigen::VectorXd result(1);
  if (surface_period[0] != 0 && std::abs(local_direction_.x()) > 1e-6 &&
      std::abs(local_direction_.y()) < 1e-6) {
    result(0) = surface_period[0] / local_direction_.x();
  } else if (surface_period[1] != 0 && std::abs(local_direction_[0]) < 1e-6 &&
             std::abs(std::abs(local_direction_[1]))) {
    result(0) = surface_period[1] / local_direction_.y();
  } else {
    // in all other cases there is no periodicity. (In the case that the
    // surface is periodic in both directions, and the direction is not along
    // one of the axis, it's almost impossible to get periodicty because of
    // floating point rounding errors.
    result(0) = 0;
  }
  return result;
}
}  // namespace lf::brep::geom
