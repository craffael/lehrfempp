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

namespace lf::brep::geom {

template <int DIM_WORLD>
class CurveStraightLine : public interface::BrepGeometry {
 public:
  CurveStraightLine(const Eigen::Matrix<double, DIM_WORLD, 1>& start,
                    const Eigen::Matrix<double, DIM_WORLD, 1>& end) {
    transform_ = end - start;
    offset_ = start;
  }

  [[nodiscard]] base::dim_t DimGlobal() const override { return DIM_WORLD; }
  [[nodiscard]] base::dim_t DimLocal() const override { return 1; }
  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] std::vector<bool> IsInBoundingBox(
      const Eigen::MatrixXd& global) const override;

  [[nodiscard]] std::pair<double, Eigen::VectorXd> Project(
      const Eigen::VectorXd& global) const override;

  [[nodiscard]] bool IsInside(const Eigen::VectorXd& local) const override;

  [[nodiscard]] Eigen::VectorXd Periods() const override {
    return Eigen::VectorXd::Zero(1);
  }

 private:
  Eigen::Matrix<double, DIM_WORLD, 1> transform_;
  Eigen::Matrix<double, DIM_WORLD, 1> offset_;
};

// explicit template instantiation declaration.
extern template class CurveStraightLine<2>;
extern template class CurveStraightLine<3>;

}  // namespace lf::brep::geom

#endif  // __b3e7e0a7a63a4981bcd5d130716484f2
