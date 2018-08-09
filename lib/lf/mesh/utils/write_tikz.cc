/** @file write_tikz.cc */

#include "write_tikz.h"
#include "lf/mesh/mesh.h"
#include "lf/base/base.h"

#include <fstream>

namespace lf::mesh::utils {

void writeTikZ(const Mesh &mesh, std::string filename){
  std::ofstream outfile(filename);

  using size_type = std::size_t; //lf::base::size_type;
  using dim_t = lf::base::RefEl::dim_t; // lf::base::dim_t;
  const Eigen::MatrixXd zero(Eigen::MatrixXd::Zero(0, 1));

  // Obtain topological dimension of the mesh
  const dim_t dim_mesh = mesh.DimMesh();
  LF_VERIFY_MSG(dim_mesh == 2, "writeTikZ() only available for 2D meshes");

  //Run through nodes
  const dim_t node_codim(dim_mesh); // Codimension number for nodes in the mesh
  const size_type no_of_nodes = mesh.Size(node_codim); // No. of nodes in codimension for nodes
  size_type node_count = 0;

  outfile << "% TikZ document graphics \n";
  outfile << "\\begin{tikzpicture}[scale=0.8]\n";

  /*
  // Loop only through nodes
  for (const Entity &node : mesh.Entities(node_codim)){
      const lf::base::glb_idx_t node_index = mesh.Index(node);
      const geometry::Geometry *geo_ptr = node.Geometry();
      Eigen::MatrixXd node_coord(geo_ptr->Global(zero));
      std::cout << node_index << std::endl;
      outfile << "\\draw (" << node_coord(0, 0) << "," << node_coord(1, 0) << ")\n";

      if((node_index+1)%no_of_nodes != 0){
          outfile << " -- ";
      }
      node_count++;
  }
  for (const Entity &e : mesh.Entities(0)){
      size_type e_idx = mesh.Index(e);
      std::cout << e_idx << "--" << std::endl;
  }
*/


  // Loop over codimensions
  for (int co_dim = dim_mesh; co_dim >= 0; co_dim--){

    // Loop over entities
    for (const Entity &e : mesh.Entities(co_dim)){
      size_type e_idx = mesh.Index(e); // Entity number/index
      //const lf::base::glb_idx_t node_index = mesh.Index(node);
      lf::base::RefEl e_refel = e.RefEl();
      std::cout << "Index for codim " << co_dim << ": " << e_idx << "(" << e_refel << ")" << std::endl;
      const geometry::Geometry *e_geo_ptr = e.Geometry();


      const Eigen::MatrixXd &ref_el_corners(e_refel.NodeCoords());
      std::cout << e_geo_ptr->Global(ref_el_corners) << std::endl;

      /*
      if (co_dim == 2){

          outfile << "\\draw (" << e_geo_ptr->Global(ref_el_corners).col(0)[0];
          << "," << e_geo_ptr->Global(ref_el_corners).col(0)[1] << ")"
          << " node[circle, draw](" << e_idx << ") {" << e_idx << "};\n";
          if( e_idx >= 1 && (e_idx)%no_of_nodes == 0){
              outfile << "\\draw [->] (" << e_idx << ") edge (" << 0 << ");\n";
          }
          else if (e_idx >= 1){
              outfile << "\\draw [->] (" << e_idx << ") edge (" << e_idx << ");\n";
          }
      }
*/

    } // loop entities
  } // loop codimensions


  outfile << "\\end{tikzpicture}" << std::endl;

} // writeTikZ mesh




} // namespace lf::mesh::utils
