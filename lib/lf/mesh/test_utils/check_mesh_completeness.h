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

/**
 * @brief check for match of entity geometries
 * @param mesh A reference to the mesh to be checked
 * @param vertices_only if true, then only the correct vertex positions
 *                 are checked, otherwise edges are included
 * @return a vector of types and indices of those entities for
 *         which the test has failed
 *
 * @note For special geometries the exact matching of cell and edge
 *       shapes may not be intended
 *
 * Additional information is logged to lf::mesh::test_utils::watertight_logger.
 */
std::vector<std::pair<lf::base::RefEl, base::glb_idx_t>> isWatertightMesh(
    const Mesh& mesh, bool vertices_only = true);

/**
 * @brief Logger that is used by isWatertightMesh() to output additional
 * information.
 */
extern std::shared_ptr<spdlog::logger> watertight_logger;

}  // namespace lf::mesh::test_utils

#endif  // __2633e7bf8dd548de839f74075687e81A
