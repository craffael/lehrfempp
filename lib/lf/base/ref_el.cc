#include "ref_el.h"

namespace lf::base {

const Eigen::MatrixXd RefEl::ncoords_point_dynamic_ = Eigen::VectorXd(0);

const Eigen::MatrixXd RefEl::ncoords_segment_dynamic_ =
    Eigen::Vector2d(0, 1).transpose();

const Eigen::MatrixXd RefEl::ncoords_tria_dynamic_ =
    (Eigen::MatrixXd{2, 3} << 0, 1, 0, 0, 0, 1).finished();

const Eigen::MatrixXd RefEl::ncoords_quad_dynamic_ =
    (Eigen::MatrixXd{2, 4} << 0, 1, 1, 0, 0, 0, 1, 1).finished();

// Print function
void PrintInfo(std::ostream &stream, const RefEl &ref_el, int output_ctrl) {
  base::dim_t dim_ref_el = ref_el.Dimension();
  base::dim_t no_nodes = ref_el.NumNodes();
  stream << "Type of reference element: " << ref_el.ToString() << std::endl;
  stream << "Dimension: " << dim_ref_el << std::endl;
  stream << "Number of nodes: " << no_nodes << std::endl;

  if (output_ctrl > 0) {
    // Loop over dimensions
    for (base::dim_t co_dim = dim_ref_el; co_dim > 0; co_dim--) {
      base::dim_t num_sub_ent = ref_el.NumSubEntities(co_dim);
      stream << "Codimension " << co_dim << " has " << num_sub_ent
             << " entities of type " << ref_el.SubType(co_dim, 0).ToString()
             << std::endl;

      if (output_ctrl > 10) {
        for (; num_sub_ent > 0; num_sub_ent--) {
          std::int32_t sub_ent = static_cast<std::int32_t>(num_sub_ent) - 1;
          stream << " Subentity " << sub_ent << " is of type "
                 << ref_el.SubType(co_dim, 0).ToString();

          if (ref_el.SubType(co_dim, 0) == RefEl::kPoint() &&
              output_ctrl > 20) {
            stream << " and has coordinates ["
                   << ref_el.NodeCoords().col(sub_ent)[0] << " "
                   << ref_el.NodeCoords().col(sub_ent)[1] << "]" << std::endl;
          }
          stream << std::endl;
        }
      }
    }
  }
}  // end void PrintInfo

}  // namespace lf::base
