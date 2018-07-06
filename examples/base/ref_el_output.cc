/** @file ref_el_output.cc
 * Demo for outputting information about reference elements
 */

#include <iostream>
#include "lf/base/base.h"

int main() {
  std::cout << "Output of information on reference elements" << std::endl;


  // Node
  auto re_node = lf::base::RefEl::kPoint();
  std::cout << "\nInformation about node: \n" << re_node << std::endl;
  //std::cout << "Sub entities of node: " << re_node.NumSubEntities(0);
  //int dim_node = re_node.Dimension();
  //std::cout << "Dimension of segment: " << dim_node << std::endl;

  // Segment
  auto re_seg = lf::base::RefEl::kSegment();
  std::cout << "\nInformation about segment: \n" << re_seg << std::endl;

  // Triangle
  auto re_tria = lf::base::RefEl::kTria();
  std::cout << "\nInformation about triangle: \n" << re_tria << std::endl;

  /*
  int dim_tria = re_tria.Dimension();
  std::cout << "Dimension: " << dim_tria << std::endl;
  std::cout << "Coordinates of nodes: " << std::endl;
  std::cout << re_tria.NodeCoords() << std::endl;



  // Loop over dimensions
  for (int co_dim = dim_tria; co_dim >= 0; co_dim--){
      int num_sub_ent = re_tria.NumSubEntities(co_dim);
      std::cout << "Codimension " << co_dim << " has " << num_sub_ent << " entities";
      std::cout << " of type " << re_tria.SubType(co_dim, 0) << std::endl;

      // Loop over entities
      for (num_sub_ent; num_sub_ent >= 0; num_sub_ent--){
          std::cout << " Subentity " << num_sub_ent << " of type " << re_tria.SubType(co_dim,0).ToString();
          std::cout << " has coordinates" << std::endl;

      }
  }
*/

/*
  // Lopp over codimensions
  for (int co_dim = dim_mesh; co_dim >= 0; co_dim--) {
    const size_type no_ent = mesh.Size(co_dim);
    o << "Co-dimension " << co_dim << ": " << no_ent << " entities"
      << std::endl;

    // Loop over entities
    for (const Entity &e : mesh.Entities(co_dim)) {
      size_type e_idx = mesh.Index(e);
      dim_t e_codim = e.Codim();
      const geometry::Geometry *e_geo_ptr = e.Geometry();
      lf::base::RefEl e_refel = e.RefEl();

      LF_VERIFY_MSG(e_geo_ptr,
                    co_dim << "-entity " << e_idx << ": missing geometry");
      LF_VERIFY_MSG(e_geo_ptr->DimLocal() == dim_mesh - co_dim,
                    co_dim << "-entity " << e_idx << ": wrong dimension");
      LF_VERIFY_MSG(e_geo_ptr->RefEl() == e_refel,
                    co_dim << "-entity " << e_idx << ": refEl mismatch");
      LF_VERIFY_MSG(e_codim == co_dim,
                    co_dim << "-entity " << e_idx << " co-dimension mismatch");
      const Eigen::MatrixXd &ref_el_corners(e_refel.NodeCoords());
      o << "entity " << e_idx << " (" << e_refel << "): ";

      // Loop over local co-dimensions
      for (int l = 1; l <= dim_mesh - co_dim; l++) {
        o << "rel codim-" << l << " subent: [";
        // Fetch subentities of co-dimension l
        base::RandomAccessRange<const Entity> sub_ent_range(e.SubEntities(l));
        for (const Entity &sub_ent : sub_ent_range) {
          o << mesh.Index(sub_ent) << ' ';
        }
        o << "], ";
      }
      o << std::endl << e_geo_ptr->Global(ref_el_corners) << std::endl;
    }  // end loop over entities
  }    // end loop over co-dimensions
*/




  // Qaud
  auto re_quad = lf::base::RefEl::kQuad();
  std::cout << "\nInformation about quad: \n" << re_quad << std::endl;


  return 0L;

}
