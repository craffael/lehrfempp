#include "ref_el.h"

namespace lf::base {

CONTROLDECLARECOMMENT(RefEl, output_ctrl_, "output_ctrl_", "Diagnostics control for RefEl");

const Eigen::MatrixXd RefEl::ncoords_point_dynamic_ = Eigen::VectorXd(0);

const Eigen::MatrixXd RefEl::ncoords_segment_dynamic_ =
    Eigen::Vector2d(0, 1).transpose();

const Eigen::MatrixXd RefEl::ncoords_tria_dynamic_ =
    (Eigen::MatrixXd{2, 3} << 0, 1, 0, 0, 0, 1).finished();

const Eigen::MatrixXd RefEl::ncoords_quad_dynamic_ =
    (Eigen::MatrixXd{2, 4} << 0, 1, 1, 0, 0, 0, 1, 1).finished();

// Print function
void PrintInfo(const RefEl &ref_el, std::ostream &o){
    o << "Print function test";
}


}  // namespace lf::base
