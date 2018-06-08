#include "tria_o1.h"
#include "point.h"
#include "segment_o1.h"

namespace lf::geometry {

TriaO1::TriaO1(Eigen::Matrix<double, Eigen::Dynamic, 3> coords)
    : coords_(std::move(coords)),
      jacobian_(coords_.rows(), 2),
      jacobian_inverse_gramian_(coords_.rows(), 2),
      integrationElement_(0) {
  jacobian_ << coords_.col(1) - coords_.col(0), coords_.col(2) - coords_.col(0);

  if (coords_.rows() == 2) {
    jacobian_inverse_gramian_ = jacobian_.transpose().inverse();
    integrationElement_ = jacobian_.determinant();
  } else {
    jacobian_inverse_gramian_ =
        jacobian_ * (jacobian_.transpose() * jacobian_).inverse();
    integrationElement_ =
        std::sqrt((jacobian_.transpose() * jacobian_).determinant());
  }
}

Eigen::MatrixXd TriaO1::Global(const Eigen::MatrixXd& local) const {
  return coords_.col(0) *
             (1 - local.array().row(0) - local.array().row(1)).matrix() +
         coords_.col(1) * local.row(0) + coords_.col(2) * local.row(1);
}

std::unique_ptr<Geometry> TriaO1::subGeometry(dim_t codim, dim_t i) const {
  using std::make_unique;
  switch (codim) {
    case 0:
      LF_ASSERT_MSG(i == 0, "i is out of bounds.");
      return std::make_unique<TriaO1>(coords_);
    case 1:
      LF_ASSERT_MSG(i >= 0 && i < 3, "i is out of bounds.");
      return make_unique<SegmentO1>(
          (Eigen::Matrix<double, Eigen::Dynamic, 2>(DimGlobal(), 2)
               << coords_.col(RefEl().SubSubEntity2SubEntity(1, i, 1, 0)),
           coords_.col(RefEl().SubSubEntity2SubEntity(1, i, 1, 1)))
              .finished());
    case 2:
      LF_ASSERT_MSG(i >= 0 && i < 3, "i is out of bounds.");
      return make_unique<Point>(coords_.col(i));
    default:
      LF_VERIFY_MSG(false, "codim is out of bounds.");
  }
}
}  // namespace lf::geometry
