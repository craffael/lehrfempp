/**
 * @file
 * @brief implementation of SurfaceCylinder class.
 * @author Raffael Casagrande
 * @date   2021-02-01 05:08:18
 * @copyright MIT License
 */

#include "surface_cylinder.h"

#include <Eigen/Geometry>
#include <iostream>

namespace lf::brep::test_utils {
SurfaceCylinder::SurfaceCylinder(Eigen::Vector3d base,
                                 const Eigen::Vector3d& axis, double radius)
    : base_(std::move(base)) {
  double height = axis.norm();
  LF_ASSERT_MSG(height > 1e-5, "Cylinder height must be > 1e-5");
  for (int i = 0; i < 3; ++i) {
    if (axis(i) > 0.1 * height) {
      mat_((i + 1) % 3, 0) = axis(i);
      mat_(i, 0) = -axis((i + 1) % 3);
      mat_((i + 2) % 3, 0) = 0.;
      break;
    }
  }
  mat_.col(0).normalize();
  mat_.col(1) = (axis / height).cross(mat_.col(0));
  mat_.col(0) *= radius;
  mat_.col(1) *= radius;
  mat_.col(2) = axis;

#ifndef NDEBUG
  Eigen::Matrix3d temp =
      Eigen::Vector3d(radius * radius, radius * radius, axis.squaredNorm())
          .asDiagonal();
  LF_ASSERT_MSG((mat_.transpose() * mat_).isApprox(temp), "Something's wrong.");
#endif
}

Eigen::MatrixXd SurfaceCylinder::GlobalMulti(
    const Eigen::MatrixXd& local) const {
  return mat_.col(0) * local.row(0).array().cos().matrix() +
         mat_.col(1) * local.row(0).array().sin().matrix() +
         mat_.col(2) * local.row(1) + base_.replicate(1, local.cols());
}

Eigen::MatrixXd SurfaceCylinder::JacobianMulti(
    const Eigen::MatrixXd& local) const {
  Eigen::MatrixXd result(3, 2 * local.cols());
  for (std::size_t i = 0; i < local.cols(); ++i) {
    result.col(i * 2) = -mat_.col(0) * std::sin(local(0, i)) +
                        mat_.col(1) * std::cos(local(0, i));
    result.col(i * 2 + 1) = mat_.col(2);
  }
  return result;
}

std::vector<bool> SurfaceCylinder::IsInBoundingBoxMulti(
    const Eigen::MatrixXd& global) const {
  double rsquared = mat_.col(0).squaredNorm();
  Eigen::MatrixXd local = Eigen::Vector3d(1. / rsquared, 1. / rsquared,
                                          1. / mat_.col(2).squaredNorm())
                              .asDiagonal() *
                          mat_.transpose() *
                          (global - base_.replicate(1, global.cols()));
  std::vector<bool> result(global.cols());
  constexpr double eps = 1e-6;
  for (int i = 0; i < global.cols(); ++i) {
    result[i] = local(0, i) >= -1 - eps && local(0, i) <= 1 + eps &&
                local(1, i) >= -1 - eps && local(1, i) <= 1 + eps &&
                local(2, i) >= -eps && local(2, i) <= 1 + eps;
  }
  return result;
}

Eigen::Vector3d SurfaceCylinder::GlobalSingle(
    const Eigen::Vector2d& local) const {
  return mat_.col(0) * std::cos(local(0)) + mat_.col(1) * std::sin(local(0)) +
         mat_.col(2) * local(1) + base_;
}

Eigen::Matrix<double, 3, 2> SurfaceCylinder::JacobianSingle(
    const Eigen::Vector2d& local) const {
  Eigen::Matrix<double, 3, 2> result;
  result.col(0) =
      -mat_.col(0) * std::sin(local(0)) + mat_.col(1) * std::cos(local(0));
  result.col(1) = mat_.col(2);
  return result;
}

std::pair<double, Eigen::Vector2d> SurfaceCylinder::Project(
    const Eigen::Vector3d& global) const {
  double rsquared = mat_.col(0).squaredNorm();
  Eigen::Vector3d local = Eigen::Vector3d(1. / rsquared, 1. / rsquared,
                                          1. / mat_.col(2).squaredNorm())
                              .asDiagonal() *
                          mat_.transpose() * (global - base_);
  std::pair<double, Eigen::Vector2d> result;
  result.second(0) = std::atan2(local.y(), local.x());
  result.second(1) = local.z();
  result.first = (GlobalSingle(result.second) - global).norm();
  return result;
}

bool SurfaceCylinder::IsInBoundingBox(const Eigen::Vector3d& global) const {
  double rsquared = mat_.col(0).squaredNorm();
  Eigen::Vector3d local = Eigen::Vector3d(1. / rsquared, 1. / rsquared,
                                          1. / mat_.col(2).squaredNorm())
                              .asDiagonal() *
                          mat_.transpose() * (global - base_);
  constexpr double eps = 1e-6;
  return local(0) >= -1 - eps && local(0) <= 1 + eps && local(1) >= -1 - eps &&
         local(1) <= 1 + eps && local(2) >= -eps && local(2) <= 1 + eps;
}

bool SurfaceCylinder::IsInside(const Eigen::Vector2d& local) const {
  constexpr double eps = 1e-6;
  return local(0) >= -base::kPi - eps && local(0) <= base::kPi &&
         local(1) >= -eps && local(1) <= 1 + eps;
}
}  // namespace lf::brep::test_utils
