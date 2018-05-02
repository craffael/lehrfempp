#include "ref_el.h"


namespace lf::base {
const std::vector<Eigen::VectorXd> RefEl::ncoords_point_dynamic_ =
  std::vector<Eigen::VectorXd>{Eigen::VectorXd::Zero(0)};

const std::vector<Eigen::VectorXd> RefEl::ncoords_segment_dynamic_ =
  std::vector<Eigen::VectorXd>{
    Eigen::VectorXd::Zero(1), Eigen::VectorXd::Constant(1, 1.)
  };

const std::vector<Eigen::VectorXd> RefEl::ncoords_tria_dynamic_ =
  std::vector<Eigen::VectorXd>{
    Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 0), Eigen::Vector2d(0, 1)
  };

const std::vector<Eigen::VectorXd> RefEl::ncoords_quad_dynamic_ =
  std::vector<Eigen::VectorXd>{
    Eigen::Vector2d(0, 0), Eigen::Vector2d(1, 0), Eigen::Vector2d(1, 1),
    Eigen::Vector2d(0, 1)
  };
}
