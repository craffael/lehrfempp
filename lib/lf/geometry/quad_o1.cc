#include "quad_o1.h"
#include "point.h"
#include "segment_o1.h"

namespace lf::geometry {

Eigen::MatrixXd QuadO1::Global(const Eigen::MatrixXd& local) const {
  return coords_.col(0) *
             ((1 - local.array().row(0)) * (1 - local.array().row(1)))
                 .matrix() +
         coords_.col(1) *
             (local.array().row(0) * (1 - local.array().col(1))).matrix() +
         coords_.col(2) *
             (local.array().row(0) * local.array().row(1)).matrix() +
         coords_.col(3) *
             ((1 - local.array().row(0)) * local.array().row(1)).matrix();
}

Eigen::MatrixXd QuadO1::Jacobian(const Eigen::MatrixXd& local) const {
  Eigen::MatrixXd result(DimGlobal(), local.cols() * 2);

  for (int i = 0; i < local.cols(); ++i) {
    result.col(2 * i) = (coords_.col(1) - coords_.col(0)) * (1 - local(1, i)) +
                        (coords_.col(2) - coords_.col(3)) * local(1, i);
    result.col(2 * i + 1) =
        (coords_.col(3) - coords_.col(0)) * (1 - local(0, i)) +
        (coords_.col(2) - coords_.col(1)) * local(0, i);
  }
  return result;
}

Eigen::MatrixXd QuadO1::JacobianInverseGramian(
    const ::Eigen::MatrixXd& local) const {
  Eigen::MatrixXd result(DimGlobal(), local.cols() * 2);
  Eigen::MatrixXd jacobian(DimGlobal(), 2);

  for (int i = 0; i < local.cols(); ++i) {
    jacobian.col(0) = (coords_.col(1) - coords_.col(0)) * (1 - local(1, i)) +
                      (coords_.col(2) - coords_.col(3)) * local(1, i);
    jacobian.col(1) = (coords_.col(3) - coords_.col(0)) * (1 - local(0, i)) +
                      (coords_.col(2) - coords_.col(1)) * local(0, i);

    if (DimGlobal() == 2) {
      result.block(0, 2 * i, DimGlobal(), 2) = jacobian.transpose().inverse();
    } else {
      result.block(0, 2 * i, DimGlobal(), 2) = (jacobian.transpose() * jacobian)
                                                   .colPivHouseholderQr()
                                                   .solve(jacobian.transpose())
                                                   .transpose();
    }
  }

  return result;
}

Eigen::VectorXd QuadO1::IntegrationElement(const Eigen::MatrixXd& local) const {
  Eigen::VectorXd result(local.cols());
  Eigen::MatrixXd jacobian(DimGlobal(), 2);
  for (int i = 0; i < local.cols(); ++i) {
    jacobian.col(0) = (coords_.col(1) - coords_.col(0)) * (1 - local(1, i)) +
                      (coords_.col(2) - coords_.col(3)) * local(1, i);
    jacobian.col(1) = (coords_.col(3) - coords_.col(0)) * (1 - local(0, i)) +
                      (coords_.col(2) - coords_.col(1)) * local(0, i);

    if (DimGlobal() == 2) {
      result(i) = jacobian.determinant();
    } else {
      result(i) = std::sqrt((jacobian.transpose() * jacobian).determinant());
    }
  }
  return result;
}

std::unique_ptr<Geometry> QuadO1::subGeometry(dim_t codim, dim_t i) const {
  using std::make_unique;
  switch (codim) {
    case 0:
      LF_ASSERT_MSG(i == 0, "i is out of bounds.");
      return std::make_unique<QuadO1>(coords_);
    case 1:
      LF_ASSERT_MSG(i >= 0 && i < 4, "i is out of bounds.");
      return make_unique<SegmentO1>(
          (Eigen::Matrix<double, Eigen::Dynamic, 2>(DimGlobal(), 2)
               << coords_.col(RefEl().SubSubEntity2SubEntity(1, i, 1, 0)),
           coords_.col(RefEl().SubSubEntity2SubEntity(1, i, 1, 1)))
              .finished());
    case 2:
      LF_ASSERT_MSG(i >= 0 && i < 4, "i is out of bounds.");
      return make_unique<Point>(coords_.col(i));
    default:
      LF_VERIFY_MSG(false, "codim is out of bounds.");
  }
}

}  // namespace lf::geometry
