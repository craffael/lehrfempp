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
 * In particular, edges, cells and nodes in the mesh can be written to file.
 * Numbering of the aforementioned and local vertice numbering of cells can be enabled by using enum flags.
 * Combine the flags by using the binary or (|) operator to get a more detailed visualization of the mesh.\n
 *
 * @note If no local numbering is wanted, simply pass 0 as the 'int output_ctrl' parameter in writeTikZ().
 *
 * #### Output control flags:
 * - TikzOutputCtrl::EdgeNumbering to display edge numbering
 * - TikzOutputCtrl::NodeNumbering to display numbering of nodes
 * - TikzOutputCtrl::CellNumbering to display numbering of cells / entities
 * - TikzOutputCtrl::VerticeNumbering to display local vertice numbering of cells
 *
 *
 * #### Example of use:
 * - How to use the function, with and without flags
 * @snippet mesh_utils.cc writeTikzUsage
 *
 * - To visualize the mesh, include the produced the file in LaTex in the following way:
 * @snippet mesh_utils.cc TikzInLatex
 *
 *
 *
 *
 *
 *
 *
 * Local vertice numbering, edge numbering and numbering of cells is possible to enable.
 *
 * #### Output control options
 * Output control is determined by enum flags. Combine several flags by the binary or operator (|).
 * Possible enum flags are:
 *
 *
 */
void writeTikZ(const lf::mesh::Mesh &mesh, std::string filename, int output_ctrl);


// enum
enum TikzOutputCtrl {

    EdgeNumbering = 1,
    CellNumbering = 2,
    NodeNumbering = 4,
    VerticeNumbering = 8
};


} // namespace lf::mesh::utils

#endif  // WRITE_TIKZ_H
