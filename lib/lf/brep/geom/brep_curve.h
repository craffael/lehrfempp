/**
 * @file
 * @brief Definition of BrepSegment
 * @author Raffael Casagrande
 * @date   2021-02-01 03:23:58
 * @copyright MIT License
 */

#ifndef __0269d75a63484458b2461314ed81c2ff
#define __0269d75a63484458b2461314ed81c2ff

#include <lf/brep/interface/interface.h>
#include <lf/geometry/geometry.h>

namespace lf::brep::geom {

class BrepCurve : public geometry::Geometry {
 public:
  BrepCurve(std::shared_ptr<const interface::BrepCurve> curve,
            Eigen::RowVector2d curve_param);

  BrepCurve(const BrepCurve&) = default;

  [[nodiscard]] dim_t DimLocal() const override { return 1; }
  [[nodiscard]] dim_t DimGlobal() const override { return curve_->DimGlobal(); }
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kSegment();
  }

  [[nodiscard]] Eigen::MatrixXd Global(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::MatrixXd Jacobian(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::MatrixXd JacobianInverseGramian(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override;

  [[nodiscard]] std::unique_ptr<Geometry> SubGeometry(dim_t codim,
                                                      dim_t i) const override;

  [[nodiscard]] std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const geometry::RefinementPattern& ref_pat,
      lf::base::dim_t codim) const override;

 private:
  std::shared_ptr<const interface::BrepCurve> curve_;
  double offset_;
  double slope_;
};

}  // namespace lf::brep::geom

#endif  // __0269d75a63484458b2461314ed81c2ff
