/**
 * @file
 * @brief A wrapper around a BRepCurve that is defined only in the x-y plane
 * @author Raffael Casagrande
 * @date   2021-02-10 03:06:54
 * @copyright MIT License
 */

#ifndef __747c21629d4242dfa617f75bc783fe53
#define __747c21629d4242dfa617f75bc783fe53

#include <lf/brep/interface/interface.h>

namespace lf::brep::geom {

class BrepCurve2d final : public interface::BrepGeometry {
 public:
  explicit BrepCurve2d(std::shared_ptr<const BrepGeometry> curve)
      : curve_(std::move(curve)) {}

  [[nodiscard]] base::dim_t DimGlobal() const override { return 2; }
  [[nodiscard]] base::dim_t DimLocal() const override { return 1; }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] std::vector<bool> IsInBoundingBox(
      const Eigen::MatrixXd& global) const override;

  [[nodiscard]] std::pair<double, Eigen::VectorXd> Project(
      const Eigen::VectorXd& global) const override;

  [[nodiscard]] bool IsInside(const Eigen::VectorXd& local) const override {
    return curve_->IsInside(local);
  }

  [[nodiscard]] Eigen::VectorXd Periods() const override {
    return curve_->Periods();
  }

 private:
  std::shared_ptr<const BrepGeometry> curve_;
};

}  // namespace lf::brep::geom

#endif  // __747c21629d4242dfa617f75bc783fe53
