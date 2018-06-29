#include "ref_el.h"

namespace lf::base {

const Eigen::MatrixXd RefEl::ncoords_point_dynamic_ = Eigen::VectorXd(0);

const Eigen::MatrixXd RefEl::ncoords_segment_dynamic_ =
    Eigen::Vector2d(0, 1).transpose();

const Eigen::MatrixXd RefEl::ncoords_tria_dynamic_ =
    (Eigen::MatrixXd{2, 3} << 0, 1, 0, 0, 0, 1).finished();

const Eigen::MatrixXd RefEl::ncoords_quad_dynamic_ =
    (Eigen::MatrixXd{2, 4} << 0, 1, 1, 0, 0, 0, 1, 1).finished();

Eigen::MatrixXd getRefElCorners(RefElType type) {
  switch (type) {
    case RefElType::kPoint:
      return RefEl::ncoords_point_dynamic_;
    case RefElType::kSegment:
      return RefEl::ncoords_segment_dynamic_;
    case RefElType::kTria:
      return RefEl::ncoords_tria_dynamic_;
    case RefElType::kQuad:
      return RefEl::ncoords_quad_dynamic_;
    default:
      LF_VERIFY_MSG(false,
                    "getRefElCorners() not implemented for this RefElType");
  }
}
}  // namespace lf::base
