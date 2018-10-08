/**
 * @file
 * @brief Declares the writeMatplotlib function which writes a mesh to a .csv
 *        file which can be visualized with matplotlib
 * @author Anian Ruoss
 * @date   2018-10-08 18:27:17
 * @copyright MIT License
 */

#ifndef WRITE_MATPLOTLIB_H
#define WRITE_MATPLOTLIB_H

#include <lf/mesh/mesh.h>

#include <string>

namespace lf::io {
/**
 * @brief Write affine triangulation data to file in matplotlib format
 * @param mesh the mesh to be stored to file
 * @param filename name of output file: .csv appended unless present
 *
 * This function creates a .csv file containing all relevant information
 * about the given mesh in the following format:
 *
 * Points:
 * codim, index, x_coord, y_coord
 *
 * Segments:
 * codim, index, point1_idx, point2_idx
 *
 * Triangles/Quadrilaterals:
 * codim, index, segment1_idx, segment2_idx, ...
 *
 * The .csv file can be read by plot_mesh.py to visualize the mesh
 */
void writeMatplotlib(const lf::mesh::Mesh &mesh, std::string filename);

}  // namespace lf::io

#endif /* WRITE_MATPLOTLIB_H */
