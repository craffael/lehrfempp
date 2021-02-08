/**
 * @file
 * @brief Implementation of Transfinite Triangle Geometry object.
 * @author Raffael Casagrande
 * @date   2021-02-03 03:25:52
 * @copyright MIT License
 */

#ifndef __6f8b2bdbdeeb4535bb4385a4a8e58fcc
#define __6f8b2bdbdeeb4535bb4385a4a8e58fcc

#include <lf/brep/interface/interface.h>
#include <lf/geometry/geometry.h>

namespace lf::brep::geom {

class BrepTriaTransfinite : public geometry::Geometry {
 public:
  BrepTriaTransfinite(const BrepTriaTransfinite&) = default;
  BrepTriaTransfinite(
      std::array<std::pair<const interface::BrepCurve*, Eigen::RowVector2d>, 3>
          curves,
      std::array<bool, 3> delete_curves);

  [[nodiscard]] dim_t DimLocal() const override { return 2; }
  [[nodiscard]] dim_t DimGlobal() const override {
    return curves_[0].first->DimGlobal();
  }
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kTria();
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

  ~BrepTriaTransfinite();

 private:
  std::array<std::pair<const interface::BrepCurve*, Eigen::RowVector2d>, 3>
      curves_;
  std::array<bool, 3> delete_curves_;
  Eigen::VectorXd node0_;
};

}  // namespace lf::brep::geom

#endif  // __6f8b2bdbdeeb4535bb4385a4a8e58fcc
