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

    int dim_ref_el = ref_el.Dimension();
    int no_nodes = ref_el.NumNodes();
    o << "Type of reference element: " << ref_el.ToString() << std::endl;
    o << "Dimension: " << dim_ref_el << std::endl;
    o << "Number of nodes: " << no_nodes << std::endl;

    if (RefEl::output_ctrl_ > 0){
        // Loop over dimensions
        for (int co_dim = dim_ref_el; co_dim > 0; co_dim--){
            int num_sub_ent = ref_el.NumSubEntities(co_dim);
            o << "Codimension " << co_dim << " has " << num_sub_ent << " entities of type " << ref_el.SubType(co_dim,0).ToString() << std::endl;

            if (RefEl::output_ctrl_ > 10){
                for (num_sub_ent; num_sub_ent > 0; num_sub_ent--){
                    int sub_ent = num_sub_ent - 1;
                    o << " Subentity " << sub_ent << " is of type " << ref_el.SubType(co_dim,0).ToString();

                    if (ref_el.SubType(co_dim,0) == RefEl::kPoint() && RefEl::output_ctrl_ > 20){
                        o << " and has coordinates [" << ref_el.NodeCoords().col(sub_ent)[0] << " " << ref_el.NodeCoords().col(sub_ent)[1] << "]" << std::endl;
                    }
                    o << std::endl;
                }
            }
        }
    }
} // end void PrintInfo

}  // namespace lf::base
