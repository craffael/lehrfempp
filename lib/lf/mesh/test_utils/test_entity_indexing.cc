#include <lf/mesh/mesh.h>
#include <iostream>

namespace lf::mesh::test_utils
{
  bool testEntityIndexing(const Mesh &mesh) {
    using size_type = Mesh::size_type; // type for indices
    const size_type dim_mesh = mesh.DimMesh(); // dimension of meshed manifold
    // Now run over all co-dimensions to check indexing
    for (int co_dim = 0; co_dim <= dim_mesh; ++co_dim) {
      // Number of entities of current co-dimension
      const size_type no_of_entities = mesh.Size(co_dim); 
      // counting array for occurrences of an index
      std::vector<int> idx_use(no_of_entities,0);
      // Traverse all entities of a given co-dimension
      for(const Entity& e : mesh.Entities(co_dim)) {
	// Fetch index of current entity
	const size_type current_idx = mesh.Index(e);
	// Check whether index is in range and, if so, count it
	if (current_idx < no_of_entities) idx_use[current_idx]++;
	else {
	  std::cerr << "Index " << current_idx
		    << " out of range" << std::endl;
	  return false;
	}
      }
      // Indices should occur only once. This means that the counting array
      // should contain 1 in every component.
      bool ok = true;
      for (int idx_cnt = 0; idx_cnt < no_of_entities; ++idx_cnt) {
	if (idx_use[idx_cnt] != 1) {
	  std::cerr << "Index " << idx_cnt << " occurs "
		    << idx_use[idx_cnt] << " times!" << std::endl;
	  ok = false;
	}
      }
      return ok;
    } // end for
  } // end testEntityIndexing()
}
