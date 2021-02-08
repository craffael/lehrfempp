/**
 * @file
 * @brief A curve that represents a straight segment
 * @author Raffael Casagrande
 * @date   2021-02-08 03:56:38
 * @copyright MIT License
 */

#ifndef __b3e7e0a7a63a4981bcd5d130716484f2
#define __b3e7e0a7a63a4981bcd5d130716484f2

#include <lf/brep/interface/interface.h>

namespace lf::brep::test_utils {

class CurveStraightLine : public interface::BrepCurve {
 public:
  CurveStraightLine(const Eigen::Vector3d& start, const Eigen::Vector3d& end) {
    transform_ = end - start;
    offset_ = start;
  }

  [[nodiscard]] base::dim_t DimGlobal() const override { return 3; }
  [[nodiscard]] base::dim_t DimLocal() const override { return 1; }
  [[nodiscard]] Eigen::MatrixXd GlobalMulti(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::MatrixXd JacobianMulti(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] std::vector<bool> IsInBoundingBoxMulti(
      const Eigen::MatrixXd& global) const override;

  [[nodiscard]] Eigen::Vector3d GlobalSingle(double local) const override {
    return transform_ * local + offset_;
  }
  [[nodiscard]] Eigen::Vector3d JacobianSingle(double local) const override {
    return transform_;
  }

  [[nodiscard]] std::pair<double, double> Project(
      const Eigen::Vector3d& global) const override;

  [[nodiscard]] bool IsInBoundingBoxSingle(
      const Eigen::Vector3d& global) const override;

  [[nodiscard]] bool IsInside(double local) const override;

 private:
  Eigen::Vector3d transform_;
  Eigen::Vector3d offset_;
};

}  // namespace lf::brep::test_utils

#endif  // __b3e7e0a7a63a4981bcd5d130716484f2
