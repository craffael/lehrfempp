/** @file write_tikz.cc */

#include "write_tikz.h"
#include "lf/base/base.h"
#include "lf/mesh/mesh.h"

#include <fstream>

namespace lf::mesh::utils {
/*
void writeTikZ(const lf::mesh::Mesh &mesh, std::string filename){

}
*/

void writeTikZ(const lf::mesh::Entity &e, std::string filename) {
  std::ofstream outfile(filename);

  outfile << "% TikZ document graphics \n";
  outfile << "\\begin{tikzpicture}[scale=2]\n";
  if (e.RefEl() == lf::base::RefEl::kPoint()) {
    outfile << "\\draw (0,0) node[circle, draw](0) {0};\n";
  }
  if (e.RefEl() == lf::base::RefEl::kSegment()) {
    outfile << "\\draw (0,0) node[circle, draw](0) {0};\n";
    outfile << "\\draw (1,0) node[circle, draw](1) {1};\n";
    outfile << "\\draw [->](0) edge (1);\n";
  }
  if (e.RefEl() == lf::base::RefEl::kTria()) {
    outfile << "\\draw (0,0) node[circle, draw](0) {0};\n";
    outfile << "\\draw (1,0) node[circle, draw](1) {1};\n";
    outfile << "\\draw (0,1) node[circle, draw](2) {2};\n";
    outfile << "\\draw [->](0) edge (1) (1) edge (2) (2) edge (0);\n";
  }
  if (e.RefEl() == lf::base::RefEl::kQuad()) {
    outfile << "\\draw (0,0) node[circle, draw](0) {0};\n";
    outfile << "\\draw (1,0) node[circle, draw](1) {1};\n";
    outfile << "\\draw (1,1) node[circle, draw](2) {2};\n";
    outfile << "\\draw (0,1) node[circle, draw}(3) {3}:\n";
    outfile
        << "\\draw [->](0) edge (1) (1) edge (2) (2) edge (3) (3) edge (0);\n";
  }
  outfile << "\\end{tikzpicture}";
}

}  // namespace lf::mesh::utils
