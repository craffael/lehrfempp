/** @file write_tikz.cc */

#include "write_tikz.h"
#include "lf/mesh/mesh.h"

#include <fstream>

namespace lf::mesh::utils {
/*
void writeTikZ(const lf::mesh::Mesh &mesh, std::string filename){

}
*/

/*
void writeTikZ(const lf::mesh::Entity &e, std::string filename){
    std::ofstream outfile(filename);

    filename << "% TikZ document graphics \n";
    filename << "\\begin{tikzpicture}[scale=2]\n";
    if(e.RefEl() == lf::RefEl::kPoint()){
        filename << "\\draw (0,0) node[circle, draw](0) {0};\n";
    }
    if(e.RefEl() == lf::RefEl::kSegment()){
        filename << "\\draw (0,0) node[circle, draw](0) {0};\n";
        filename << "\\draw (1,0) node[circle, draw](1) {1};\n";
        filename << "\\draw [->](0) edge (1);\n";
    }
    filename << "\\end{tikzpicture}"
}
*/



} // namespace lf::mesh::utils
