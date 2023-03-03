#ifndef INCG2633e7bf8dd548de839f74075687e810
#define INCG2633e7bf8dd548de839f74075687e810

#include <lf/mesh/mesh.h>

namespace lf::mesh::test_utils {
/**
 * @brief Function for testing mesh indexing (should be called from google test)
 * @param mesh A reference to the mesh to be checked
 *
 * This function tests whether all entities of a mesh are indexed
 * consecutively
 */
void checkEntityIndexing(const Mesh& mesh);

}  // namespace lf::mesh::test_utils

#endif  // INCG2633e7bf8dd548de839f74075687e810
