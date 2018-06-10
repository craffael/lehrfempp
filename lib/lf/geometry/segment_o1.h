#ifndef __0ec30f6fc5554ec9bea4d9b83e83ea42
#define __0ec30f6fc5554ec9bea4d9b83e83ea42

#include "geometry_interface.h"

namespace lf::geometry {

class SegmentO1 : public Geometry {
 public:
  explicit SegmentO1(Eigen::Matrix<double, Eigen::Dynamic, 2> coords)
      : coords_(std::move(coords)) {}

  dim_t DimLocal() const override { return 1; }
  dim_t DimGlobal() const override { return coords_.rows(); }
  base::RefEl RefEl() const override { return base::RefEl::kSegment(); }
  Eigen::MatrixXd Global(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd Jacobian(const Eigen::MatrixXd& local) const override;
  Eigen::MatrixXd JacobianInverseGramian(
      const ::Eigen::MatrixXd& local) const override;
  Eigen::VectorXd IntegrationElement(
      const Eigen::MatrixXd& local) const override;
  std::unique_ptr<Geometry> subGeometry(dim_t codim, dim_t i) const override;

 private:
  Eigen::Matrix<double, Eigen::Dynamic, 2> coords_;
};

}  // namespace lf::geometry

#endif  // __0ec30f6fc5554ec9bea4d9b83e83ea42
