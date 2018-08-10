/** @file write_tikz.cc */

#include "write_tikz.h"
#include "lf/mesh/mesh.h"
#include "lf/base/base.h"

#include <fstream>

namespace lf::mesh::utils {

void writeTikZ(const Mesh &mesh, std::string filename, int output_ctrl){
  std::ofstream outfile(filename);

  // For the enum flags: TikzOutputCtrl
  bool EdgeNumOn = output_ctrl & TikzOutputCtrl::EdgeNumbering;
  bool NodeNumOn = output_ctrl & TikzOutputCtrl::NodeNumbering;
  bool CellNumOn = output_ctrl & TikzOutputCtrl::CellNumbering;

  if(EdgeNumOn){
      // Display edge numbering
  }
  // Similar goes for NodeNumOn and CellNumOn

  //
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


  // START writing to file
  outfile << "% TikZ document graphics \n";
  outfile << "\\begin{tikzpicture}[scale=0.8]\n";



  /*
  for (const Entity &cell : mesh.Entities(0)){
      size_type cell_idx = mesh.Index(cell);
      lf::base::RefEl cell_refel = cell.RefEl();
      int num_nodes_cell = cell_refel.NumNodes();
      const geometry::Geometry *cell_geo_ptr = cell.Geometry();
      const Eigen::MatrixXd &cell_corners(cell_refel.NodeCoords());
      //Eigen::MatrixXd node_coord(cell_geo_ptr->Global(zero));


      outfile << "\\draw ";
      for (int node = 0; node < num_nodes_cell; node++){
          outfile << "(" << cell_geo_ptr->Global(cell_corners).col(node)[0]
                  << "," << cell_geo_ptr->Global(cell_corners).col(node)[1] << ") -- ";

      } // for nodes
      outfile << "(" << cell_geo_ptr->Global(cell_corners).col(num_nodes_cell)[0]
              << "," << cell_geo_ptr->Global(cell_corners).col(num_nodes_cell)[1] << ");\n";




  } // for cells
*/

  outfile << "\\end{tikzpicture}" << std::endl;


} // writeTikZ mesh




} // namespace lf::mesh::utils
