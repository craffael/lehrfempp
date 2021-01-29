/**
 * @file
 * @brief Definition of a BrepCurve that describes a circle
 * @author Raffael Casagrande
 * @date   2021-01-29 08:29:24
 * @copyright MIT License
 */

#ifndef __85d9c8572d5541a48402af491100e932
#define __85d9c8572d5541a48402af491100e932

#include <lf/brep/interface/interface.h>

namespace lf::brep::test_utils {

/**
 * @brief Represents a circle in the x-y plane with given `origin` and `radius`.
 * The parameter range is [-pi,pi]
 */
class CurveCircle : public interface::BrepCurve {
 public:
  CurveCircle(const Eigen::Vector3d& origin, double radius);

  [[nodiscard]] base::dim_t DimGlobal() const override { return 3; }
  [[nodiscard]] base::dim_t DimLocal() const override { return 1; }
  [[nodiscard]] Eigen::MatrixXd GlobalMulti(
      const Eigen::MatrixXd& local) const override;
  [[nodiscard]] Eigen::MatrixXd JacobianMulti(
      const Eigen::MatrixXd& local) const override;
  [[nodiscard]] std::vector<bool> IsInBoundingBoxMulti(
      const Eigen::MatrixXd& global) const override;
  [[nodiscard]] Eigen::Vector3d GlobalSingle(double local) const override;
  [[nodiscard]] Eigen::Vector3d JacobianSingle(double local) const override;
  [[nodiscard]] std::pair<double, double> Project(
      const Eigen::Vector3d& global) const override;
  [[nodiscard]] bool IsInBoundingBoxSingle(
      const Eigen::Vector3d& global) const override;
  [[nodiscard]] bool IsInside(double local) const override;

 private:
  Eigen::Vector3d origin_;
  double radius_;
};

}  // namespace lf::brep::test_utils

#endif  // __85d9c8572d5541a48402af491100e932
