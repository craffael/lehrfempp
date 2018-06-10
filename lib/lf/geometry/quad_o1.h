#ifndef __b193e889f4e245959d124b0ded8a3b68
#define __b193e889f4e245959d124b0ded8a3b68

#include "geometry_interface.h"

namespace lf::geometry {

class QuadO1 : public Geometry {
 public:
  explicit QuadO1(Eigen::Matrix<double, Eigen::Dynamic, 4> coords)
      : coords_(std::move(coords)) {}

  dim_t DimLocal() const override { return 2; }
  dim_t DimGlobal() const override { return coords_.rows(); }
  base::RefEl RefEl() const override { return base::RefEl::kQuad(); }

  Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd JacobianInverseGramian(
      const ::Eigen::MatrixXd& local) const override;
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override;
  std::unique_ptr<Geometry> subGeometry(dim_t codim, dim_t i) const override;

 private:
  Eigen::Matrix<double, Eigen::Dynamic, 4> coords_;
};

}  // namespace lf::geometry

#endif  // __b193e889f4e245959d124b0ded8a3b68
