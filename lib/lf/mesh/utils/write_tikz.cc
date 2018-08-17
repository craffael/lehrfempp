/** @file write_tikz.cc */

#include "write_tikz.h"
#include "lf/base/base.h"
#include "lf/mesh/mesh.h"

#include <fstream>

namespace lf::mesh::utils {

void writeTikZ(const Mesh &mesh, std::string filename, int output_ctrl){
  std::ofstream outfile(filename);

  // For the enum flags: TikzOutputCtrl
  bool EdgeNumOn = output_ctrl & TikzOutputCtrl::EdgeNumbering;
  bool NodeNumOn = output_ctrl & TikzOutputCtrl::NodeNumbering;
  bool CellNumOn = output_ctrl & TikzOutputCtrl::CellNumbering;
  bool VerticeNumOn = output_ctrl & TikzOutputCtrl::VerticeNumbering;
  bool RenderCellsOn = output_ctrl & TikzOutputCtrl::RenderCells;


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
  outfile << "\\begin{tikzpicture}[scale=4, >= stealth, inner sep=0pt, minimum size=0.4cm]\n";
  // outfile << "\\draw[grey, very thin] (-4.3,-4.3) grid (4.3, 4.3)\n";

  for (int co_dim = 0; co_dim <= dim_mesh; co_dim++){

      for (const Entity &obj : mesh.Entities(co_dim)){
          size_type obj_idx = mesh.Index(obj);
          lf::base::RefEl obj_refel = obj.RefEl();
          int num_nodes_obj = obj_refel.NumNodes();
          const geometry::Geometry *obj_geo_ptr = obj.Geometry();
          const Eigen::MatrixXd &obj_corners(obj_refel.NodeCoords());
          const Eigen::MatrixXd vertices = obj_geo_ptr->Global(obj_corners);
          const Eigen::MatrixXd center = vertices.rowwise().sum()/vertices.cols();
          Eigen::MatrixXd center_mat(center.rows(), num_nodes_obj);


          switch (obj_refel) {

          case lf::base::RefEl::kPoint():{

              if(NodeNumOn){
                  outfile << "\\draw[red, fill = white] (" << vertices(0,0) << "," << vertices(1,0) << ") "
                          << "node[draw, circle, fill = white] {" << obj_idx << "};\n";
              } else {
                  outfile << "\\draw[red] (" << vertices(0,0) << "," << vertices(1,0) << ") "
                          << "node[] {*};\n";
              } // if NodeNumOn
              break;

          } // case kPoint

          case lf::base::RefEl::kSegment():{
              center_mat << center, center;
              const Eigen::MatrixXd scaled_vertices = vertices * 0.85 + center_mat * 0.15;
              const Eigen::MatrixXd semi_scaled_vertices = vertices * 0.95 + center_mat * 0.05;

              if (EdgeNumOn && NodeNumOn){
                  outfile << "\\draw[->] ("
                          << scaled_vertices(0,0) << "," << scaled_vertices(1,0) << ") -- node[black] {\\huge" << obj_idx << "} "
                          << "(" << scaled_vertices(0,1) << "," << scaled_vertices(1,1) << ");\n";
              } else if (NodeNumOn && !EdgeNumOn){
                  outfile << "\\draw[->] ("
                          << scaled_vertices(0,0) << "," << scaled_vertices(1,0) << ") -- "
                          << "(" << scaled_vertices(0,1) << "," << scaled_vertices(1,1) << ");\n";
              } else if (!NodeNumOn && EdgeNumOn){
                  outfile << "\\draw[->] ("
                          << semi_scaled_vertices(0,0) << "," << semi_scaled_vertices(1,0) << ") -- node[black] {" << obj_idx << "} "
                          << "(" << semi_scaled_vertices(0,1) << "," << semi_scaled_vertices(1,1) << ");\n";
              } else if(!NodeNumOn && !EdgeNumOn){
                  outfile << "\\draw[->] ("
                          << semi_scaled_vertices(0,0) << "," << semi_scaled_vertices(1,0) << ") -- "
                          << "(" << semi_scaled_vertices(0,1) << "," << semi_scaled_vertices(1,1) << ");\n";
              } else {std::cout << "Check EdgeNumOn and NodeNumOn for kSegment " << obj_idx << std::endl;}
              break;
          } // case kSegment

          case lf::base::RefEl::kTria():{
              center_mat << center, center, center;
              const Eigen::MatrixXd scaled_vertices = vertices * 0.8 + center_mat * 0.2;

              if(RenderCellsOn){

                  if(VerticeNumOn){
                      outfile << "\\draw[green] ("
                              << scaled_vertices(0,0) << "," << scaled_vertices(1,0) << ") node[] {0} -- ("
                              << scaled_vertices(0,1) << "," << scaled_vertices(1,1) << ") node[] {1} -- ("
                              << scaled_vertices(0,2) << "," << scaled_vertices(1,2) << ") node[] {2} -- cycle;\n";
                  } else {
                      outfile << "\\draw[green] ("
                              << scaled_vertices(0,0) << "," << scaled_vertices(1,0) << ") -- ("
                              << scaled_vertices(0,1) << "," << scaled_vertices(1,1) << ") -- ("
                              << scaled_vertices(0,2) << "," << scaled_vertices(1,2) << ") -- cycle;\n";
                  } // if EdgeNumOn

                  if(CellNumOn){
                      outfile << "\\draw[green] (" << center(0,0) << "," << center(1,0) << ") node[] {" << obj_idx << "};\n";
                  }

              } // RenderCellsOn
              break;
          } // case kTria

          case lf::base::RefEl::kQuad():{
              center_mat << center, center, center, center;
              const Eigen::MatrixXd scaled_vertices = vertices * 0.8 + center_mat * 0.2;

              if(RenderCellsOn){

                  if (VerticeNumOn){
                      outfile << "\\draw[magenta] ("
                              << scaled_vertices(0,0) << "," << scaled_vertices(1,0) << ") node[] {0} -- ("
                              << scaled_vertices(0,1) << "," << scaled_vertices(1,1) << ") node[] {1} -- ("
                              << scaled_vertices(0,2) << "," << scaled_vertices(1,2) << ") node[] {2} -- ("
                              << scaled_vertices(0,3) << "," << scaled_vertices(1,3) << ") node[] {3} -- cycle;\n";
                  } else {
                      outfile << "\\draw[magenta] ("
                              << scaled_vertices(0,0) << "," << scaled_vertices(1,0) << ") -- ("
                              << scaled_vertices(0,1) << "," << scaled_vertices(1,1) << ") -- ("
                              << scaled_vertices(0,2) << "," << scaled_vertices(1,2) << ") -- ("
                              << scaled_vertices(0,3) << "," << scaled_vertices(1,3) << ") -- cycle;\n";
                  }

                  if(CellNumOn){
                      outfile << "\\draw[magenta] (" << center(0,0) << "," << center(1,0) << ") node[] {" << obj_idx << "};\n";
                  }

              } //RenderCellsOn

              break;
          } // case kQuad

          default:{
              std::cout << "Error for object " << obj_idx << " in co-dim " << co_dim << std::endl;
              std::cout << "Object type: " << obj_refel << std::endl;
              break;
          } // default
          } // switch


      } // for entities

       node_count++;
       //LF_VERIFY_MSG(node_count == no_of_nodes, "Node count mismatch");
   } // for codim

  outfile << "\\end{tikzpicture}\n";


} // writeTikZ mesh

} // namespace lf::mesh::utils
