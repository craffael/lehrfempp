/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Functions for the initialization of special data sets
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#ifndef _LF_SPEC_DATA_
#define _LF_SPEC_DATA_

#include "all_codim_mesh_data_set.h"
#include "codim_mesh_data_set.h"

namespace lf::mesh::utils {

/**
 * @brief store number of adjacent super-entities
 *
 * @param mesh_p reference to underlying mesh
 * @param codim_sub co-dimension of the queried entities
 * @param codim_super _relative_ co-dimension (with positive sign) of
 *                      super entities.
 * @return a cardinal-valued @ref CodimMeshDataSet (= an array of cardinals
 * indexed by entities of a particular co-dimension) storing adjacency numbers
 *
 * For each entity of a given co-dimension, this function counts the number of
 * adjacent super-entities of some smaller co-dimension.
 *
 * ### Example
 *
 * If, for a 2D mesh we want to count the number of cells adjacent to edges,
 * we have to specify codim_sub = 1, codim_super = 1!
 *
 * If, for a 2D mesh you want to count the number of cells owning a node,
 * specify codim_sub = 2 and codim_super = 2!
 *
 * @note codim_super is _relative_ to codim_sub with flipped sign!
 */
CodimMeshDataSet<lf::base::size_type> CountNumSuperEntities(
    const std::shared_ptr<const Mesh>& mesh_p, lf::base::dim_t codim_sub,
    lf::base::dim_t codim_super);

/**
 * @brief flag entities of a specific co-dimension located on the boundary
 *
 * @param codim co-dimension of entities to be flagged, must be > 0.
 * @return an object of a boolean-valued @ref CodimMeshDataSet (= an array of
 * boolean values index by entities) for the entities of the specified
 * co-dimension
 *
 * An entity of co-dimension 1 is located on the boundary, if it is adjacent
 * to exactly 1 cell (= entity of co-dimension 0).
 *
 * The boundary of a mesh is the set of all entities that are either entities
 * of co-dimension 1 located on the boundary or sub-entities of those.
 *
 * The implementation of this  function relies on @ref CountNumSuperEntities().
 *
 * The following example code shows how to create a flag array for marking the
 * boundary edges of a 2D mesh:
 * ~~~
   std::shared_ptr<lf::mesh::Mesh> mesh_p{
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3)};
   lf::mesh::utils::CodimMeshDataSet<bool> bd_flags{
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
 * ~~~
 */
CodimMeshDataSet<bool> flagEntitiesOnBoundary(
    const std::shared_ptr<const Mesh>& mesh_p, lf::base::dim_t codim);

/**
 * @brief flag entities of _any co-dimension_ located on the boundary
 *
 * @param codim co-dimension of entities to be flagged, must be > 0.
 * @return an array of boolean values (actually an object of type @ref
 * AllCodimMeshDataSet) for the entities of the specified co-dimension
 *
 * An entity of co-dimension 1 is located on the boundary, if it is adjacent
 * to exactly 1 cell (= entity of co-dimension 0).
 *
 * The boundary of a mesh is the set of all entities that are either entities
 * of co-dimension 1 located on the boundary or sub-entities of those.
 *
 */
AllCodimMeshDataSet<bool> flagEntitiesOnBoundary(
    const std::shared_ptr<const Mesh>& mesh_p);

}  // namespace lf::mesh::utils

#endif
