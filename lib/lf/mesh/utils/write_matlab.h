#ifndef __X633e7bf8dd548de839f74075687e81A
#define __X633e7bf8dd548de839f74075687e81A

#include <lf/mesh/mesh.h>

#include <string>

namespace lf::mesh::utils {
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
 * - _TRI_ an N x 3 -matrix whose rows contain the indices of the nodes
 *   of every triangle of the mesh
 * - _EDS_ an Mx2 - matrix containing the indices of the endpoints of
 *   the edges of the mesh.
 *
 * The data _x_, _y_, _TRI_ can be processed by MATLAB's `triplot` and
 * 'trisurf' functions.
 *
 */
void writeMatlab(const lf::mesh::Mesh &mesh, std::string filename);

}  // namespace lf::mesh::utils

#endif  // __X633e7bf8dd548de839f74075687e81A
