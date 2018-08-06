/** @file write_tikz.cc */

#include "write_tikz.h"
#include "lf/mesh/mesh.h"
#include "lf/base/base.h"

#include <fstream>

namespace lf::mesh::utils {

void writeTikZ(const Mesh &mesh, std::string filename){

  using size_type = std::size_t; //lf::base::size_type;
  using dim_t = lf::base::RefEl::dim_t; // lf::base::dim_t;
  const Eigen::MatrixXd zero(Eigen::MatrixXd::Zero(0, 1));


  // Obtain topological dimension of the mesh
  const dim_t dim_mesh = mesh.DimMesh();
  //const dim_t dim_world = mesh.DimWorld();
  LF_VERIFY_MSG(dim_mesh == 2, "writeTikZ() only available for 2D meshes");

  const dim_t node_codim(dim_mesh);
  const size_type no_of_nodes = mesh.Size(node_codim);





  std::ofstream outfile(filename);

  outfile << "% TikZ document graphics \n";
  outfile << "\\begin{tikzpicture}[scale=2]\n";

  // Loop over codimensions
  for (int co_dim = dim_mesh; co_dim >= 0; co_dim--){
    // Loop over entities
    for (const Entity &e : mesh.Entities(co_dim)){
      size_type e_idx = mesh.Index(e); // Entity number/index
      dim_t e_codim = e.Codim();
      const geometry::Geometry *e_geo_ptr = e.Geometry();
      lf::base::RefEl e_refel = e.RefEl();

      const Eigen::MatrixXd &ref_el_corners(e_refel.NodeCoords());
      if(co_dim == 2){
          outfile << "\\draw (" << e_geo_ptr->Global(ref_el_corners).col(0)[0]
          << "," << e_geo_ptr->Global(ref_el_corners).col(0)[1] << ")"
          << " node[circle, draw](" << e_idx << ") {" << e_idx << "};\n";
          if( e_idx >= 1 && (e_idx)%no_of_nodes == 0){
              outfile << "\\draw [->] (" << e_idx << ") edge (" << 0 << ");\n";
          }
          else if (e_idx >= 1){
              outfile << "\\draw [->] (" << e_idx << ") edge (" << e_idx << ");\n";
          }
      }

    } // loop entities
  } // loop codimensions


  outfile << "\\end{tikzpicture}" << std::endl;

} // writeTikZ mesh



/*
void writeTikZ(const lf::mesh::Entity &e, std::string filename){
    std::ofstream outfile(filename);

    outfile << "% TikZ document graphics \n";
    outfile << "\\begin{tikzpicture}[scale=2]\n";
    if(e.RefEl() == lf::base::RefEl::kPoint()){
        outfile << "\\draw (0,0) node[circle, draw](0) {0};\n";
    }
    if(e.RefEl() == lf::base::RefEl::kSegment()){
        outfile << "\\draw (0,0) node[circle, draw](0) {0};\n";
        outfile << "\\draw (1,0) node[circle, draw](1) {1};\n";
        outfile << "\\draw [->](0) edge (1);\n";
    }
    if(e.RefEl() == lf::base::RefEl::kTria()){
        outfile << "\\draw (0,0) node[circle, draw](0) {0};\n";
        outfile << "\\draw (1,0) node[circle, draw](1) {1};\n";
        outfile << "\\draw (0,1) node[circle, draw](2) {2};\n";
        outfile << "\\draw [->](0) edge (1) (1) edge (2) (2) edge (0);\n";
    }
    if(e.RefEl() == lf::base::RefEl::kQuad()){
        outfile << "\\draw (0,0) node[circle, draw](0) {0};\n";
        outfile << "\\draw (1,0) node[circle, draw](1) {1};\n";
        outfile << "\\draw (1,1) node[circle, draw](2) {2};\n";
        outfile << "\\draw (0,1) node[circle, draw}(3) {3}:\n";
        outfile << "\\draw [->](0) edge (1) (1) edge (2) (2) edge (3) (3) edge (0);\n";

    }
    outfile << "\\end{tikzpicture}";

  
} // end void writetikz
*/





} // namespace lf::mesh::utils
