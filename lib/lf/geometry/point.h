#ifndef __184598b89ca44fe1a1e7a043bc32da06
#define __184598b89ca44fe1a1e7a043bc32da06

#include <lf/base/base.h>

#include "geometry.h"

namespace lf::geometry {

class Point : public Geometry {
 public:
  explicit Point(Eigen::VectorXd coord) : coord_(std::move(coord)) {}

  [[nodiscard]] dim_t DimLocal() const override { return 0; }

  [[nodiscard]] dim_t DimGlobal() const override { return coord_.rows(); }

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kPoint();
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

  /**
   * @brief the child geometry is just a copy of the point geometry
   */
  [[nodiscard]] std::vector<std::unique_ptr<Geometry>> ChildGeometry(
      const RefinementPattern& ref_pattern, base::dim_t codim) const override;

 private:
  Eigen::VectorXd coord_;
};

}  // namespace lf::geometry

#endif  // __184598b89ca44fe1a1e7a043bc32da06
