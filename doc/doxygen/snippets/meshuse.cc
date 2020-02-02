/**
 * @file
 * @brief Doxygen snippets to show lf::assemble::COOMatrix usage
 * @author Ralf Hiptmair
 * @date Tue Nov 5 2019
 * @copyright MIT License
 */

#include <iostream>
#include <lf/mesh/mesh.h>

namespace lf::mesh {
// This snippet is partly copied from
// lf::mesh::test_utils::checkEntityIndexing() file: check_entity_indexing.cc
//! [usage]
bool checkEntityIndexing(const lf::mesh::Mesh& mesh) {
  // dimension of meshed manifold
  const lf::base::dim_t dim_mesh = mesh.DimMesh();
  // Now run over all co-dimensions to check indexing
  for (lf::base::dim_t co_dim = 0; co_dim <= dim_mesh; ++co_dim) {
    // Number of entities of current co-dimension. For a 2D hybrid mesh:
    // codim == 0: cells, codim == 1: edges, codim == 2: nodes
    const lf::base::size_type no_of_entities = mesh.NumEntities(co_dim);
    // counting array for occurrences of an index
    std::vector<int> idx_use(no_of_entities, 0);
    // Traverse all entities of a given co-dimension. Note that the loop
    // variable is a pointer!
    for (const lf::mesh::Entity* entity_p : mesh.Entities(co_dim)) {
      // Fetch index of current entity
      const lf::base::glb_idx_t current_idx = mesh.Index(*entity_p);
      // Check whether index is in range and, if so, count it
      if (current_idx < no_of_entities)
        idx_use[current_idx]++;
      else
        return false;
      // Check consistency of indexing
      const lf::mesh::Entity* e_ptr = mesh.EntityByIndex(co_dim, current_idx);
      if ((e_ptr == nullptr) || (entity_p != e_ptr) ||
          (current_idx != mesh.Index(*e_ptr)))
        return false;
    }  // end loop over entities

    // Indices should occur only once. This means that the counting array
    // should contain 1 in every component.
    for (lf::base::size_type idx_cnt = 0; idx_cnt < no_of_entities; ++idx_cnt) {
      if (idx_use[idx_cnt] != 1) return false;
    }
  }  // end for
  return true;
}  // end checkEntityIndexing()
   //! [usage]

  // This snippet is copied from lecturedemomesh.cc
  //! [loopentity]
  int traverseEntities(const lf::mesh::Mesh &mesh, lf::base::dim_t codim) {
    lf::base::size_type cnt = 0;
  // Typical loop for running through all entities of a specific co-dimension
  for (const lf::mesh::Entity *entity : mesh.Entities(codim)) {
    // Print entity information including its unique index
    std::cout << cnt << ": Entity #" << mesh.Index(*entity) << ": " << *entity
              << std::endl;
    cnt++;
  }
  return cnt;
}
  //! [loopentity]
}  // namespace lf::mesh
