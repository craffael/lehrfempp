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
 * @brief Represents a curve embedded in a interface::BrepGeometry surface (as a
 * straight line in parameter space).
 */
class BrepSurfaceSegment : public interface::BrepGeometry {
 public:
  BrepSurfaceSegment(std::shared_ptr<const interface::BrepGeometry> surface,

                     const Eigen::Matrix<double, 2, 2>& surface_local_coords)
      : surface_(surface),
        local_direction_(surface_local_coords.col(1) -
                         surface_local_coords.col(0)),
        local_offset_(surface_local_coords.col(0)) {
    LF_ASSERT_MSG(surface_->DimLocal() == 2, "surface must have 2d DimLocal.");
  }

  [[nodiscard]] inline base::dim_t DimGlobal() const override {
    return surface_->DimGlobal();
  }
  [[nodiscard]] inline base::dim_t DimLocal() const override { return 1; }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override;
  [[nodiscard]] std::vector<bool> IsInBoundingBox(
      const Eigen::MatrixXd& global) const override {
    // suboptimal, but everything else is very complicated.
    return surface_->IsInBoundingBox(global);
  }

  [[nodiscard]] std::pair<double, Eigen::VectorXd> Project(
      const Eigen::VectorXd& global) const override;
  [[nodiscard]] bool IsInside(const Eigen::VectorXd& local) const override;

  [[nodiscard]] Eigen::VectorXd Periods() const override;

 private:
  std::shared_ptr<const interface::BrepGeometry> surface_;
  Eigen::Vector2d local_direction_;
  Eigen::Vector2d local_offset_;
};

}  // namespace lf::brep::geom

#endif  // __32d19b9405d342b5a804e93c64fd2e1e
