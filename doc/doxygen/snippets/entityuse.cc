/**
 * @file
 * @brief Doxygen snippets to show lf::assemble::COOMatrix usage
 * @author Ralf Hiptmair
 * @date Wed Nov 6 2019
 * @copyright MIT License
 */

#include <lf/mesh/mesh.h>

namespace lf::mesh {
// Copied from lf::mesh::test_utils::checkLocalTopolgy
//! [usage]
bool checkLocalTopology(const Entity &e) {
  // Obtain basic information about current Entity
  const lf::base::RefEl ref_el = e.RefEl();
  const lf::base::dim_t dimension = ref_el.Dimension();
  // What this function does in the case of a 2D cell entity (dimension = 2):
  // It runs through all edges (co-dimension = 1), fetches their endnodes
  // (co-dimension again 1 relative to an edge) and test whether they
  // are also sub-entities of co-dimension 2 of the cell.

  // Co-dimensions of sub-entities run from 1 to dimension
  for (lf::base::dim_t sub_codim = 1; sub_codim <= dimension; sub_codim++) {
    // Number of sub-entities with relative co-dimensjon sub_codim
    const lf::base::size_type num_sub_entities =
        ref_el.NumSubEntities(sub_codim);
    // Obtain sequence of pointers to sub-entities of co-dimension sub_codim
    std::span<const Entity *const> sub_ent_range = e.SubEntities(sub_codim);
    // Run through sub-entities
    lf::base::size_type sub_ent_cnt{0};
    for (const lf::mesh::Entity *sub_ent : sub_ent_range) {
      const lf::base::RefEl sub_ref_el = sub_ent->RefEl();
      const lf::base::dim_t sub_dim = sub_ref_el.Dimension();
      if (sub_dim != dimension - sub_codim) return false;
      // The sub-entity has further sub-entities of codimension 1 to sub_dim
      for (lf::base::dim_t sub_sub_codim = 1; sub_sub_codim <= sub_dim;
           sub_sub_codim++) {
        std::span<const Entity *const> sub_sub_ent_range =
            sub_ent->SubEntities(sub_sub_codim);
        for (const lf::mesh::Entity *sub_sub_ent : sub_sub_ent_range) {
          // The entity pointed to by sub_sub_ent has co-dimension
          // sub_codim + sub_sub_codim w.r.t. the entity referenced by e
          // Hence get all corresponding sub-entities of e
          std::span<const Entity *const> e_sub_sub_range =
              e.SubEntities(sub_codim + sub_sub_codim);
          // See whether we find sub_sub_ent in this range
          int found = 0;  // Count how many times we find the sub-entity
          for (const Entity *e_sub_sub : e_sub_sub_range) {
            if (e_sub_sub == sub_sub_ent) found++;
          }
          // Any sub-entity should occur exactly once.
          if (found != 1) return false;
        }
      }
      sub_ent_cnt++;
    }  // end loop over sub-entities
    if (num_sub_entities != sub_ent_cnt) return false;
  }
  return true;
}  // end checklocalTopology
   //! [usage]
}  // namespace lf::mesh
