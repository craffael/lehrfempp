#ifndef __2633e7bf8dd548de839f74075687e81A
#define __2633e7bf8dd548de839f74075687e81A

#include <lf/mesh/mesh.h>

namespace lf::mesh::test_utils {
/**
 * @brief Function testing topological completeness of a mesh
 * @param mesh A reference to the mesh to be checked
 *
 * All entities of dimension dim of a mesh must be subentities
 * of at least one entity of dimension dim+1, if dim is smaller
 * than the topological dimension of the mesh
 */
bool checkMeshCompleteness(const Mesh& mesh);

}  // namespace lf::mesh::test_utils

#endif  // __2633e7bf8dd548de839f74075687e81A
