/**
 * @file
 * @brief Represents a cylindrical surface in 3d space.
 * @author Raffael Casagrande
 * @date   2021-02-01 05:07:14
 * @copyright MIT License
 */

#ifndef __7ddb7ad3a39245928724e60211a096b7
#define __7ddb7ad3a39245928724e60211a096b7

#include <lf/brep/interface/interface.h>

#include <utility>

namespace lf::brep::test_utils {

class SurfaceCylinder : public interface::BrepSurface {
 public:
  SurfaceCylinder(Eigen::Vector3d base, const Eigen::Vector3d& axis,
                  double radius);

  [[nodiscard]] base::dim_t DimGlobal() const override { return 3; }
  [[nodiscard]] base::dim_t DimLocal() const override { return 2; }
  [[nodiscard]] Eigen::MatrixXd GlobalMulti(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::MatrixXd JacobianMulti(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] std::vector<bool> IsInBoundingBoxMulti(
      const Eigen::MatrixXd& global) const override;

  [[nodiscard]] Eigen::Vector3d GlobalSingle(
      const Eigen::Vector2d& local) const override;

  [[nodiscard]] Eigen::Matrix<double, 3, 2> JacobianSingle(
      const Eigen::Vector2d& local) const override;

  [[nodiscard]] std::pair<double, Eigen::Vector2d> Project(
      const Eigen::Vector3d& global) const override;

  [[nodiscard]] bool IsInBoundingBox(
      const Eigen::Vector3d& global) const override;

  [[nodiscard]] bool IsInside(const Eigen::Vector2d& local) const override;

 private:
  Eigen::Matrix3d mat_;
  Eigen::Vector3d base_;
};

}  // namespace lf::brep::test_utils

#endif  // __7ddb7ad3a39245928724e60211a096b7
