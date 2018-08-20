#ifndef WRITE_TIKZ_H
#define WRITE_TIKZ_H

#include <lf/mesh/mesh.h>
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/hybrid2dp/hybrid2dp.h"

#include <string>

namespace lf::mesh::utils {

/**
 * @brief Writes mesh to file in TikZ Graphics format. File as input in LaTeX will draw the mesh.
 *
 *
 * @param mesh the mesh to be stored to file
 * @param filename name of output file.
 * @param output_ctrl enum flags controlling amount of output
 *
 * This function writes a file of code, which included in LaTeX draws a
 * mesh using TikZ Graphics.
 * In particular, edges, cells and nodes in the mesh can be written to file and visualized.
 * Numbering of the aforementioned and local vertice numbering of cells can be enabled by using enum flags.
 * Combine the flags by using the binary or (|) operator to get a more detailed visualization of the mesh.\n
 *
 *
 * #### Output control flags:
 * - TikzOutputCtrl::RenderCells to show the specific cells in the mesh in addition to the mesh grid
 * - TikzOutputCtrl::CellNumbering to display numbering of cells
 * - TikzOutputCtrl::VerticeNumbering to display local vertice numbering of cells
 * - TikzOutputCtrl::NodeNumbering to display numbering of nodes
 * - TikzOutputCtrl::EdgeNumbering to display edge numbering
 *
 *
 * @note If no details about nodes, cells or edges are wanted, simply pass 0 as the 'int output_ctrl' parameter in writeTikZ(). This will draw only the mesh grid and nodes.
 * @note TikzOutputCtrl::RenderCells must be enabled in order to use the flags for numbering of cells and of local vertices of cells.
 *
 * In the LaTeX document, remember to include "\usepackage{tikz}". Use "\input{}" to include the code file and visualize the mesh.
 *
 *
 *
 */
void writeTikZ(const lf::mesh::Mesh &mesh, std::string filename, int output_ctrl);
/* // Things that can be added in write_tikz documentation
 *
 * #### Example of use:
 * - How to use the function, with and without flags
 * @snippet mesh_utils.cc writeTikzUsage
 *
 * - To visualize the mesh, include the produced the file in LaTex in the following way:
 * @snippet mesh_utils.cc TikzInLatex
 *
*/



/*
/**
 * @brief writeTikZ: second version!
 * @param mesh
 * @param filename
 * @param selector
 * @param output_ctrl
 */

//void writeTikZ(const lf::mesh::Mesh &mesh, std::string filename,
//               std::function<bool(const lf::mesh::Entity &)> selector, int output_ctrl = 7);




/**
 * @brief Enum flags: TikzOutputCtrl for output control of mesh drawn in TikZ.
 */
enum TikzOutputCtrl {
    RenderCells = 1,
    CellNumbering = 2,
    VerticeNumbering = 4,
    NodeNumbering = 8,
    EdgeNumbering = 16
};


} // namespace lf::mesh::utils

#endif  // WRITE_TIKZ_H
