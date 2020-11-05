#ifndef LF_LD_MESH_H
#define LF_LD_MESH_H
/**
 * @file
 * @brief Functions for simple LehrFEM++ demos + sample codes
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include <filesystem>

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include "lecturedemo.h"
#include "lf/mesh/utils/lambda_mesh_data_set.h"

namespace lecturedemo {

/** @brief traverses entities of specified co-dimension and
 *         prints information about them.
 * @param mesh reference to a read-only LehrFEM++ mesh
 * @param codim co-dimension of interest.
 */
int traverseEntities(const lf::mesh::Mesh &mesh, dim_t codim);

/** @brief function counting cell objects of different topological type
 *  @param mesh reference to a read-only LehrFEM++ mesh
 */
std::pair<size_type, size_type> countCellTypes(const lf::mesh::Mesh &mesh);

/** @brief function printing the topological relationships of mesh entities
 *  @param mesh read-only LehrFEM++ mesh object
 *  @param codim co-dimension of interest.
 */
void scanTopology(const lf::mesh::Mesh &mesh, dim_t codim);

/**
 * @brief output of the corners of the cells of a mesh
 * @param mesh const reference to a LehrFEM++ mesh object
 * @codim co-dimension of entities of interest
 *
 * Prints the coordinates of all corners for all entities of
 * the specified co-dimension.
 */
void PrintGeometryInfo(const lf::mesh::Mesh &mesh, dim_t codim);

/** @brief driver routine for LehrFEM++ demos for lecture
 */
void lecturedemomesh();

};  // namespace lecturedemo

#endif  // LF_LD_MESH_H
