#include "check_local_topology.h"
#include <gtest/gtest.h>

namespace lf::mesh::test_utils {
void checkLocalTopology(const Entity &e) {
  using dim_t = lf::base::RefEl::dim_t;
  using RefEl = lf::base::RefEl;
  using size_type = lf::mesh::Mesh::size_type;

  // Obtain basic information about current Entity
  const RefEl ref_el = e.RefEl();
  const dim_t dimension = ref_el.Dimension();

  // What this function does in the case of a 2D cell entity (dimension = 2)
  // It runs through all edges (co-dimension = 1), fetches their endnodes
  // (co-dimension again 1 relative to an edge) and test whether they
  // are also sub-entities of co-dimension 2 of the cell.

  // Co-dimensions of sub-entities run from 1 to dimension
  for (dim_t sub_codim = 1; sub_codim <= dimension; sub_codim++) {
    // Number of sub-entities with relative co-dimensjon sub_codim
    const size_type num_sub_entities = ref_el.NumSubEntities(sub_codim);
    // Obtain sequence of sub-entities of co-dimensjon sub_codim
    base::RandomAccessRange<const Entity> sub_ent_range(
        e.SubEntities(sub_codim));
    // Run through sub-entities
    for (const Entity &sub_ent : sub_ent_range) {
      const RefEl sub_ref_el = sub_ent.RefEl();
      const dim_t sub_dim = sub_ref_el.Dimension();
      EXPECT_EQ(sub_dim, dimension - sub_codim) << "Dim/Co-dim mismatch";
      // The sub-entity has further sub-entities of codimension 1 to sub_dim
      for (dim_t sub_sub_codim = 1; sub_sub_codim <= sub_dim; sub_sub_codim++) {
        base::RandomAccessRange<const Entity> sub_sub_ent_range(
            sub_ent.SubEntities(sub_sub_codim));
        for (const Entity &sub_sub_ent : sub_sub_ent_range) {
          // The entity referenced by sub_sub_ent has co-dimension
          // sub_codim + sub_sub_codim w.r.t. the entity referenced by e
          // Hence get all corresponding sub-entities of e
          base::RandomAccessRange<const Entity> e_sub_sub_range(
              e.SubEntities(sub_codim + sub_sub_codim));
          // See whether we find sub_sub_ent in this range
          int found = 0;  // Count how many times we find the sub-entity
          for (const Entity &e_sub_sub : e_sub_sub_range) {
            if (e_sub_sub == sub_sub_ent) {
              found++;
            }
          }
          EXPECT_EQ(found, 1) << "Sub-sub-entity hit " << found << " times";
        }
      }
    }
  }
}

}  // namespace lf::mesh::test_utils
