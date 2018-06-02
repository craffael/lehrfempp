#include "segment_o1.h"
#include "point.h"

namespace lf::geometry {

Eigen::MatrixXd SegmentO1::Global(const Eigen::MatrixXd& local) const {
  return coords_.col(1) * local + coords_.col(0) *
         (Eigen::ArrayXd::Ones(1, local.cols()) - local.array()).matrix();
}

Eigen::MatrixXd SegmentO1::Jacobian(const Eigen::MatrixXd& local) const {
  return (coords_.col(1) - coords_.col(0)).replicate(1, local.cols());
}

Eigen::MatrixXd SegmentO1::JacobianInverseGramian(
  const ::Eigen::MatrixXd& local) const {
  if (DimGlobal() == 1) {
    return (coords_.col(1) - coords_.col(0)).cwiseInverse()
                                            .replicate(1, local.cols());
  } else {
    return ((coords_.col(1) - coords_.col(0))
            / (coords_.col(1) - coords_.col(0))
            .squaredNorm()).replicate(1, local.cols());
  }
}

Eigen::VectorXd SegmentO1::IntegrationElement(
  const Eigen::MatrixXd& local) const {
  return Eigen::VectorXd::Constant(local.cols(),
                                   (coords_.col(1) - coords_.col(0)).norm());
}

std::unique_ptr<Geometry> SegmentO1::subGeometry(dim_t codim, dim_t i) const {
  if (codim == 0) {
    LF_ASSERT_MSG(i == 0, "i is out of bounds.");
    return std::make_unique<SegmentO1>(coords_);
  } else if (codim == 1) {
    LF_ASSERT_MSG(i >= 0 && i<2, "i is out of bounds.");
    return std::make_unique<Point>(coords_.col(i));
  } else {
    LF_VERIFY_MSG(false, "codim is out of bounds.");
  }
}


}
