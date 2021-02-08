/**
 * @file
 * @brief Definition of BrepSurfaceCurve
 * @author Raffael Casagrande
 * @date   2021-01-29 04:10:53
 * @copyright MIT License
 */

#ifndef __32d19b9405d342b5a804e93c64fd2e1e
#define __32d19b9405d342b5a804e93c64fd2e1e

#include <lf/brep/interface/interface.h>

namespace lf::brep::geom {

/**
 * @brief Represents a curve embedded in a interface::BrepSurface (as a straight
 * line in parameter space).
 */
class BrepSurfaceSegment : public interface::BrepCurve {
 public:
  BrepSurfaceSegment(const interface::BrepSurface* surface,
                     const Eigen::Matrix<double, 2, 2>& surface_local_coords)
      : surface_(surface),
        local_direction_(surface_local_coords.col(1) -
                         surface_local_coords.col(0)),
        local_offset_(surface_local_coords.col(0)) {}

  [[nodiscard]] inline base::dim_t DimGlobal() const override {
    return surface_->DimGlobal();
  }
  [[nodiscard]] inline base::dim_t DimLocal() const override { return 1; }
  [[nodiscard]] Eigen::MatrixXd GlobalMulti(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::MatrixXd JacobianMulti(
      const Eigen::MatrixXd& local) const override;
  [[nodiscard]] std::vector<bool> IsInBoundingBoxMulti(
      const Eigen::MatrixXd& global) const override {
    // suboptimal, but everything else is very complicated.
    return surface_->IsInBoundingBoxMulti(global);
  }
  [[nodiscard]] Eigen::Vector3d GlobalSingle(double local) const override {
    return surface_->GlobalSingle(local_direction_ * local + local_offset_);
  }
  [[nodiscard]] Eigen::Vector3d JacobianSingle(double local) const override {
    return surface_->JacobianSingle(local_direction_ * local + local_offset_) *
           local_direction_;
  }

  [[nodiscard]] std::pair<double, double> Project(
      const Eigen::Vector3d& global) const override;
  [[nodiscard]] bool IsInBoundingBoxSingle(
      const Eigen::Vector3d& global) const override;
  [[nodiscard]] bool IsInside(double local) const override;

 private:
  const interface::BrepSurface* surface_;
  Eigen::Vector2d local_direction_;
  Eigen::Vector2d local_offset_;
};

}  // namespace lf::brep::geom

#endif  // __32d19b9405d342b5a804e93c64fd2e1e
