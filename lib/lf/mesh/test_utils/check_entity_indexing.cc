#include <gtest/gtest.h>
#include <lf/mesh/mesh.h>
#include <iostream>

namespace lf::mesh::test_utils {
void checkEntityIndexing(const Mesh& mesh) {
  using size_type = Mesh::size_type;     // type for indices
  const auto dim_mesh = mesh.DimMesh();  // dimension of meshed manifold
  // Now run over all co-dimensions to check indexing
  for (size_type co_dim = 0; co_dim <= dim_mesh; ++co_dim) {
    // Number of entities of current co-dimension
    const auto no_of_entities = mesh.Size(co_dim);
    // counting array for occurrences of an index
    std::vector<int> idx_use(no_of_entities, 0);
    // Traverse all entities of a given co-dimension
    for (const Entity& e : mesh.Entities(co_dim)) {
      // Fetch index of current entity
      const auto current_idx = mesh.Index(e);
      // Check whether index is in range and, if so, count it
      EXPECT_LT(current_idx, no_of_entities)
          << "Index " << current_idx << " out of range";
      idx_use[current_idx]++;
    }
    // Indices should occur only once. This means that the counting array
    // should contain 1 in every component.
    for (size_type idx_cnt = 0; idx_cnt < no_of_entities; ++idx_cnt) {
      // NOLINTNEXTLINE
      EXPECT_EQ(idx_use[idx_cnt], 1)
          << "Index " << idx_cnt << " occurs " << idx_use[idx_cnt] << " times!";
    }
  }  // end for
}  // end checkEntityIndexing()
}  // namespace lf::mesh::test_utils
