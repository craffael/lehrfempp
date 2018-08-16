
#include "lf/mesh/mesh.h"
#include "lf/mesh/utils/utils.h"
#include <cassert>
#include <iostream>

namespace lf::mesh::utils {
void foo() {

    //![writeTikzUsage]
    // Given a mesh object named mesh
    // With the use of the enum flag for node numbering
    // writeTikZ(*mesh, "filename.txt", TikzOutputCtrl::NodeNumbering);

    // Combining enum flags, enabling more detailed output
    // writeTikZ(*mesh, "filename.txt", TikzOutputCtrl::NodeNumbering|TikzOutputCtrl::EdgeNumbering|TikzOutputCtrl::VerticeNumbering);

    // Without flags

    //![writeTikzUsage]


    //![TikzInLatex]

        // \documentclass{article}
        // \usepackage{tikz}
        // \begin{document}

        // \input{"filename.txt"}

        // \end{document}
    //![TikzInLatex]

} // foo_utils
}  // namespace lf::base


