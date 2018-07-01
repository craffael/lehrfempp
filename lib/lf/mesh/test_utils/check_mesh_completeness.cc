#include <gtest/gtest.h>
#include <lf/mesh/mesh.h>
#include <iostream>

namespace lf::mesh::test_utils {
bool checkMeshCompleteness(const Mesh& mesh) {
  bool status = true;                    // Test passed or not ?
  using size_type = Mesh::size_type;     // type for indices
  using dim_t = lf::base::RefEl::dim_t;  // type for dimensions
  // Obtain topological dimension of the mesh
  const dim_t dim_mesh = mesh.DimMesh();
  // Now run over all entities of co-dimension < dim_mesh
  for (size_type co_dim = 0; co_dim < dim_mesh; ++co_dim) {
    // Count occurrences of sub-entities of relative co-dimension 1
    // To that end allocate a vector of counters
    std::vector<size_type> entity_link_cnt(mesh.Size(co_dim + 1), (size_type)0);

    // Diagnostic output
    // std::cout << "co-dim " << co_dim + 1 << ": " << mesh.Size(co_dim + 1)
    //           << " entities" << std::endl;

    // Traverse all entities of a given co-dimension
    for (const Entity& e : mesh.Entities(co_dim)) {
      // Diagnostic output
      // std::cout << "Entity(" << mesh.Index(e) << "): " << std::flush;

      // Fetch subentities of co-dimension 1
      base::RandomAccessRange<const Entity> sub_ent_range(e.SubEntities(1));
      for (const Entity& sub_ent : sub_ent_range) {
        // Diagnostic output
        // std::cout << mesh.Index(sub_ent) << " " << std::flush;
        // Obtain index of the sub-entity to address counter
        entity_link_cnt[mesh.Index(sub_ent)]++;
      }
      // For diagnostic output
      // std::cout << std::endl;
    }  // end loop over entities

    // Diagnostic output
    // for(int j=0; j < entity_link_cnt.size() ; j++)
    //   std::cout << "Entity " << j << " -> " << entity_link_cnt[j] << " links"
    //   << std::endl;

    // Maximal number of occurrences of a subentity, this many bins for counting
    size_type max_subent_cnt =
        *std::max_element(entity_link_cnt.begin(), entity_link_cnt.end());
    // for (size_type cnt : entity_link_cnt) {
    //   if (cnt > max_subent_cnt) max_subent_cnt = cnt;
    // }
    std::vector<size_type> occurrence_cnt(max_subent_cnt + 1, (size_type)0);
    size_type entity_index = 0;
    for (size_type i : entity_link_cnt) {
      occurrence_cnt[i]++;
      EXPECT_GT(i, 0) << "Entity " << entity_index << ", co-dimension "
                      << co_dim + 1 << "not linked";
      if (i == 0) status = false;
      entity_index++;
    }
    // Output of diagnostic information
    // Should depend on some control variable
    std::cout << "Enties of dimension " << dim_mesh - co_dim - 1 << ": "
              << std::endl;
    for (int l = 0; l <= max_subent_cnt; l++) {
      std::cout << l << " times linked: " << occurrence_cnt[l] << " entities"
                << std::endl;
    }
  }  // end loop over co-dimensions
  return status;
}  // end checkMeshCompleteness

}  // namespace lf::mesh::test_utils
