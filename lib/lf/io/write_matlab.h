#ifndef __X633e7bf8dd548de839f74075687e81A
#define __X633e7bf8dd548de839f74075687e81A

#include <lf/mesh/mesh.h>

#include <string>

namespace lf::io {
/**
 * @brief Writes affine triangulation data to file in MATLAB format
 *
 * @param mesh the mesh to be stored to file
 * @param filename name of output file: .m will be appended unless present
 *
 * This function creates a .m-file containing a MATLAB function
 * of the same name
 *
 *    function [x,y,TRI,EDS] = filename()
 *
 * that initializes four variables
 * - _x_ the x-coordinates of the nodes if the affine triangular mesh
 * - _y_ the y-coordinates of the nodes if the affine triangular mesh
 * - _TRI_ an N x 4 -matrix whose rows contain the indices of the nodes
 *   of every triangle of the mesh in positions 1-3, the triangle index in
 *   position 4.
 * - _QUAD_ an N x 5 -matrix; each row contains the node indices of a quad as
 *    first four entries, and the cell index as last entry.
 * - _EDS_ an Mx2 - matrix containing the indices of the endpoints of
 *   the edges of the mesh.
 *
 * The data returned by this function can be visualized by the MATLAB
 * function plot_lf_mesh(), which is also  provided with LehrFEM++.
 *
 */
void writeMatlab(const lf::mesh::Mesh &mesh, std::string filename);

}  // namespace lf::io

#endif  // __X633e7bf8dd548de839f74075687e81A
