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


    // Keep anyway
    int dim_ref_el = ref_el.Dimension();
    int no_nodes = ref_el.NumNodes();

    // Remove when done
    o << "Dimension: " << dim_ref_el << std::endl;
    o << "Number of nodes: " << no_nodes << std::endl;
    //o << "Coordinates of nodes: " << std::endl;
    //o << ref_el.NodeCoords().col(0)[0] << std::endl;



    // Loop over dimensions
    for (int co_dim = dim_ref_el; co_dim > 0; co_dim--){
        int num_sub_ent = ref_el.NumSubEntities(co_dim);
        o << "Codimension " << co_dim << " has " << num_sub_ent << " entities of type " << ref_el.SubType(co_dim,0).ToString() << std::endl;

        // Loop over entities
        for (num_sub_ent; num_sub_ent > 0; num_sub_ent--){
            int sub_ent = num_sub_ent - 1;
            o << " Subentity " << sub_ent << " of type " << ref_el.SubType(co_dim,0).ToString();
            if (ref_el.SubType(co_dim,0) == RefEl::kPoint()){
                o << " has coordinates [" << ref_el.NodeCoords().col(sub_ent)[0] << " " << ref_el.NodeCoords().col(sub_ent)[1] << "]" << std::endl;
            }
            o << std::endl;
        }

    }


    // Actual implementation
    if (RefEl::output_ctrl_ > 0){
        o << "Reference element: " << ref_el.ToString() << std::endl;
        o << "Dimension: " << dim_ref_el << std::endl;
        o << "Number of nodes: " << no_nodes << std::endl;



    } else if (RefEl::output_ctrl_ > 10){
        o << "output_ctrl_ > 10" << std::endl;
        o << "Reference element: " << ref_el.ToString() << std::endl;
        o << "Dimension: " << dim_ref_el << std::endl;
        o << "Number of nodes: " << no_nodes << std::endl;
        o << "Coordinates of nodes: " << std::endl;
        o << ref_el.NodeCoords() << std::endl;

    } else if (RefEl::output_ctrl_ > 20){
        o << "Reference element: " << ref_el.ToString() << std::endl;
        o << "Dimension: " << dim_ref_el << std::endl;
        o << "Number of nodes: " << no_nodes << std::endl;
        //o << "Coordinates of nodes: " << std::endl;
        //o << ref_el.NodeCoords() << std::endl;

    } else if (RefEl::output_ctrl_ > 30){
        o << "Reference element: " << ref_el.ToString() << std::endl;
        o << "Dimension: " << dim_ref_el << std::endl;
        o << "Number of nodes: " << no_nodes << std::endl;
        //o << "Coordinates of nodes: " << std::endl;
        //o << ref_el.NodeCoords() << std::endl;

    } else {
        // o << "Choose another value for output_ctrl_ > 0";
        o << " ";
    }

}

}  // namespace lf::base
