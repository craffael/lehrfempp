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

#include "special_entity_sets.h"

namespace lf::mesh::utils {
CodimMeshDataSet<lf::base::size_type> countNoSuperEntities(
    const std::shared_ptr<const Mesh>& mesh_p, lf::base::dim_t codim_sub,
    lf::base::dim_t codim_super) {
  LF_VERIFY_MSG((mesh_p->DimMesh() >= codim_sub),
                "Illegal codim_sub = " << codim_sub);
  LF_VERIFY_MSG(codim_super <= codim_sub, "Codim_super to large");

  // Declare and initialize the data set
  CodimMeshDataSet<lf::base::size_type> sup_ent_cnt{mesh_p, codim_sub, 0};

  const lf::base::dim_t super_codim = codim_sub - codim_super;
  // Run through all super entities
  for (const lf::mesh::Entity& e : mesh_p->Entities(super_codim)) {
    // Traverse all sub-entities of a specific relative co-dimension
    for (const lf::mesh::Entity& subent : e.SubEntities(codim_super)) {
      sup_ent_cnt(subent) += 1;
    }
  }
  return sup_ent_cnt;
}

CodimMeshDataSet<bool> flagEntitiesOnBoundary(
    const std::shared_ptr<const Mesh>& mesh_p, lf::base::dim_t codim) {
  LF_ASSERT_MSG((codim > 0) && (codim <= mesh_p->DimMesh()),
                "Illegal codim = " << codim);
  // count cells adjacent to entities of co-dimension 1
  CodimMeshDataSet<lf::base::size_type> no_adjacent_cells{
      countNoSuperEntities(mesh_p, 1, 1)};
  // flag array
  CodimMeshDataSet<bool> bd_flags{mesh_p, codim, false};
  // relative codimension with respect to faces (entities of co-dimensio 1)
  const lf::base::dim_t rel_codim = codim - 1;
  // Run through  faces and flag sub-entities
  for (const lf::mesh::Entity& edge : mesh_p->Entities(1)) {
    if (no_adjacent_cells(edge) == 1) {
      // Boundary face detected!
      // Traverse all sub-entities of a specific relative co-dimension
      for (const lf::mesh::Entity& subent : edge.SubEntities(rel_codim)) {
        bd_flags(subent) = true;
      }
    }
  }
  return bd_flags;
}

AllCodimMeshDataSet<bool> flagEntitiesOnBoundary(
    const std::shared_ptr<const Mesh>& mesh_p) {
  // count cells adjacent to entities of co-dimension 1
  CodimMeshDataSet<lf::base::size_type> no_adjacent_cells{
      countNoSuperEntities(mesh_p, 1, 1)};
  // flag array
  AllCodimMeshDataSet<bool> bd_flags{mesh_p, false};

  // Run through faces (entities of co-dimension 1) and flag sub-entities
  for (const lf::mesh::Entity& edge : mesh_p->Entities(1)) {
    if (no_adjacent_cells(edge) == 1) {
      // Boundary face detected!
      // Traverse all sub-entities of all co-dimensions
      const lf::base::dim_t dim_mesh = mesh_p->DimMesh();
      for (lf::base::dim_t rel_codim = 0; rel_codim < dim_mesh; rel_codim++) {
        for (const lf::mesh::Entity& subent : edge.SubEntities(rel_codim)) {
          bd_flags(subent) = true;
        }
      }
    }
  }
  return bd_flags;
}

}  // namespace lf::mesh::utils
