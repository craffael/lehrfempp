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
  bool VerticeNumOn = output_ctrl & TikzOutputCtrl::VerticeNumbering;

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
  outfile << "\\begin{tikzpicture}[scale=3]\n";

  for (int co_dim = dim_mesh; co_dim >= 0; co_dim--){

      for (const Entity &obj : mesh.Entities(co_dim)){
          size_type obj_idx = mesh.Index(obj);
          lf::base::RefEl obj_refel = obj.RefEl();
          const geometry::Geometry *obj_geo_ptr = obj.Geometry();
          const Eigen::MatrixXd &obj_corners(obj_refel.NodeCoords());

          switch (obj_refel) {
          case lf::base::RefEl::kPoint():
              outfile << "\\draw[red, fill = white] (" << obj_geo_ptr->Global(obj_corners).col(0)[0] << "," << obj_geo_ptr->Global(obj_corners).col(0)[1] << ") node[circle, draw, fill = white] {";

              // Node numbering
              if(NodeNumOn){
                  outfile << obj_idx << "};\n";
              } else {outfile << "};\n"; }
              break;

          case lf::base::RefEl::kSegment():
              if (EdgeNumOn){
                  outfile << "\\draw[->, black] (" << obj_geo_ptr->Global(obj_corners).col(0)[0]
                          << "," << obj_geo_ptr->Global(obj_corners).col(0)[1] << ") -- node[] {" << obj_idx << "} "
                          << "(" << obj_geo_ptr->Global(obj_corners).col(1)[0]
                          << "," << obj_geo_ptr->Global(obj_corners).col(1)[1] << ");\n";
              } else {
                  outfile << "\\draw[->, black] (" << obj_geo_ptr->Global(obj_corners).col(0)[0]
                          << "," << obj_geo_ptr->Global(obj_corners).col(0)[1] << ")"
                          << "(" << obj_geo_ptr->Global(obj_corners).col(1)[0]
                          << "," << obj_geo_ptr->Global(obj_corners).col(1)[1] << ");\n";
              }
              break;

          case lf::base::RefEl::kTria():
              if(VerticeNumOn){
                  outfile << "\\draw[->, black] (" << obj_geo_ptr->Global(obj_corners).col(0)[0]
                          << "," << obj_geo_ptr->Global(obj_corners).col(0)[1] << ") node[very near end, below left] {0} -- ("
                          << obj_geo_ptr->Global(obj_corners).col(1)[0] << "," << obj_geo_ptr->Global(obj_corners).col(1)[1] << ") node[very near end, below left] {1} -- ("
                          << obj_geo_ptr->Global(obj_corners).col(2)[0] << "," << obj_geo_ptr->Global(obj_corners).col(2)[1] << ") node[very near end, below left] {2} -- cycle;\n";
              } else {
                  outfile << "\\draw[->, black] (" << obj_geo_ptr->Global(obj_corners).col(0)[0]
                          << "," << obj_geo_ptr->Global(obj_corners).col(0)[1] << ") -- ("
                          << obj_geo_ptr->Global(obj_corners).col(1)[0] << "," << obj_geo_ptr->Global(obj_corners).col(1)[1] << ") -- ("
                          << obj_geo_ptr->Global(obj_corners).col(2)[0] << "," << obj_geo_ptr->Global(obj_corners).col(2)[1] << ") -- cycle;\n";

              } // if EdgeNumOn
              break;

          case lf::base::RefEl::kQuad():
              if (VerticeNumOn){
                  outfile << "\\draw[->, black] (" << obj_geo_ptr->Global(obj_corners).col(0)[0]
                          << "," << obj_geo_ptr->Global(obj_corners).col(0)[1] << ") node[very near end, below left] {0} -- ("
                          << obj_geo_ptr->Global(obj_corners).col(1)[0] << "," << obj_geo_ptr->Global(obj_corners).col(1)[1] << ") node[very near end, below left] {1} -- ("
                          << obj_geo_ptr->Global(obj_corners).col(2)[0] << "," << obj_geo_ptr->Global(obj_corners).col(2)[1] << ") node[very near end, below left] {2} -- ("
                          << obj_geo_ptr->Global(obj_corners).col(3)[0] << "," << obj_geo_ptr->Global(obj_corners).col(3)[1] << ") node[very near end, below left] {3} -- cycle;\n";
              } else {
                  outfile << "\\draw[->, black] (" << obj_geo_ptr->Global(obj_corners).col(0)[0]
                          << "," << obj_geo_ptr->Global(obj_corners).col(0)[1] << ") -- ("
                          << obj_geo_ptr->Global(obj_corners).col(1)[0] << "," << obj_geo_ptr->Global(obj_corners).col(1)[1] << ") -- ("
                          << obj_geo_ptr->Global(obj_corners).col(2)[0] << "," << obj_geo_ptr->Global(obj_corners).col(2)[1] << ") -- ("
                          << obj_geo_ptr->Global(obj_corners).col(3)[0] << "," << obj_geo_ptr->Global(obj_corners).col(3)[1] << ") -- cycle;\n";
              }

          default:
              std::cout << "Error for object " << obj_idx << " in co-dim " << co_dim << std::endl;
              std::cout << "Object type: " << obj_refel << std::endl;
              break;
          } // switch
      } // for entities

       node_count++;
   } // for codim

  outfile << "\\end{tikzpicture}\n";


} // writeTikZ mesh


} // namespace lf::mesh::utils
