#ifndef WRITE_TIKZ_H
#define WRITE_TIKZ_H

#include <lf/mesh/mesh.h>
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/hybrid2dp/hybrid2dp.h"

#include <string>
#include <vector>

namespace lf::mesh::utils {

/**
 * @brief Writes mesh to file in TikZ Graphics format. File as input in LaTeX
will draw the mesh.
 *
 *
 * @param mesh the mesh to be stored to file
 * @param filename name of output file.
 * @param output_ctrl enum flags controlling amount of output
 *
 * This function writes a file of code, which included in LaTeX draws a
 * mesh using TikZ Graphics.
 * In particular, edges, cells and nodes in the mesh can be written to file and
visualized.
 * Numbering of the aforementioned and local vertice numbering of cells can be
enabled by using enum flags.
 * Combine the flags by using the binary or (|) operator to get a more detailed
visualization of the mesh.\n
 * Another option is to pass the corresponding integer value directly as an
argument. See the enum definition for correct value.\n
 *
 * #### Output control flags:
 * - TikzOutputCtrl::RenderCells to show the specific cells in the mesh in
addition to the mesh grid
 * - TikzOutputCtrl::CellNumbering to display numbering of cells
 * - TikzOutputCtrl::VerticeNumbering to display local vertice numbering of
cells
 * - TikzOutputCtrl::NodeNumbering to display numbering of nodes
 * - TikzOutputCtrl::EdgeNumbering to display edge numbering
 *
 *
 * @note If no details about nodes, cells or edges are wanted, simply pass 0 as
the 'int output_ctrl' parameter in writeTikZ(). This will draw only the mesh
grid and nodes.
 * @note Omitting the output_ctrl argument completely when calling the function:
Cells, numbering of cells and numbering of vertices will be printed.
(output_ctrl = 7)
 * @note TikzOutputCtrl::RenderCells must be enabled in order to use the flags
for numbering of cells and of local vertices of cells.
 *
 * In the LaTeX document, remember to include "\usepackage{tikz}". Use
"\input{}" to include the code file and visualize the mesh.
 *
 * #### Examples of use
 *
 *
 * ##### Function call
 *
 * \verbatim

    // Enum flag for node numbering
    writeTikZ(*mesh, "filename.txt", TikzOutputCtrl::NodeNumbering);

    // Combining enum flags, enabling more detailed output
    // The two examples are equivalent:
    writeTikZ(*mesh, "filename.txt",
TikzOutputCtrl::RenderCells|TikzOutputCtrl::EdgeNumbering|TikzOutputCtrl::VerticeNumbering);
    writeTikZ(*mesh, "filename.txt", 21);
    // Note that ::VerticeNumbering only works because ::RenderCells is
activated.

    // Without flags
    writeTikZ(*mesh, "filename.txt",0);

    // Without specifying last argument
    writeTikZ(*mesh, "filename.txt"); is equivalent to writeTikZ(*mesh,
"filename.txt", 7);

 \endverbatim
 *
 * ##### LaTeX input
 * \verbatim

     \documentclass{article}
     \usepackage{tikz}
     \begin{document}

     \input{"filename.txt"}

     \end{document}
\endverbatim
 *
 *
 */
void writeTikZ(const lf::mesh::Mesh &mesh, const std::string &filename,
               int output_ctrl = 7);

/**
 * @brief writeTikZ: second version!
 * @param mesh The mesh to be stored to file
 * @param filename Name of output file
 * @param selector
 * @param output_ctrl
 */
void writeTikZ2(const lf::mesh::Mesh &mesh, const std::string &filename,
                std::function<bool(const lf::mesh::Entity &)> selector,
                int output_ctrl = 7);

//

/**
 * @brief Enum flags: TikzOutputCtrl for output control of mesh drawn in TikZ.
 */
enum TikzOutputCtrl {
  RenderCells = 1,
  CellNumbering = 2,
  VerticeNumbering = 4,
  NodeNumbering = 8,
  EdgeNumbering = 16,
  ArrowTips = 32
};

}  // namespace lf::mesh::utils

#endif  // WRITE_TIKZ_H
