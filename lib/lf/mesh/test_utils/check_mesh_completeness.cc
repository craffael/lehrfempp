#include "check_mesh_completeness.h"

#include <gtest/gtest.h>
#include <lf/mesh/mesh.h>
#include <iostream>
#include "lf/mesh/hybrid2d/mesh.h"

#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>

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
    std::vector<size_type> entity_link_cnt(mesh.NumEntities(co_dim + 1), 0);

    // Diagnostic output
    // std::cout << "co-dim " << co_dim + 1 << ": " << mesh.Size(co_dim + 1)
    //           << " entities" << std::endl;

    // Traverse all entities of a given co-dimension
    for (const Entity* e : mesh.Entities(co_dim)) {
      // Diagnostic output
      // std::cout << "Entity(" << mesh.Index(e) << "): " << std::flush;

      // Fetch subentities of co-dimension 1
      auto sub_ent_range = e->SubEntities(1);
      for (const Entity* sub_ent : sub_ent_range) {
        // Diagnostic output
        // std::cout << mesh.Index(sub_ent) << " " << std::flush;
        // Obtain index of the sub-entity to address counter
        entity_link_cnt[mesh.Index(*sub_ent)]++;
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
    std::vector<size_type> occurrence_cnt(max_subent_cnt + 1, 0);
    size_type entity_index = 0;
    for (size_type i : entity_link_cnt) {
      occurrence_cnt[i]++;
      EXPECT_GT(i, 0) << "Entity " << entity_index << ", co-dimension "
                      << co_dim + 1 << "not linked";
      if (i == 0) {
        status = false;
      }
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

std::shared_ptr<spdlog::logger>& WatertightLogger() {
  static auto logger =
      base::InitLogger("lf::mesh::test_utils::WatertightLogger");
  return logger;
}

std::vector<std::pair<lf::base::RefEl, base::glb_idx_t>> isWatertightMesh(
    const Mesh& mesh, bool vertices_only) {
  base::dim_t dim_mesh = mesh.DimMesh();
  std::vector<std::pair<lf::base::RefEl, base::glb_idx_t>> ret_vals{};

  // "Reference coordinates" for a point: dummy argument
  Eigen::Matrix<double, 0, 1> pt_ref_coord;

  // Loop over cells and edges
  for (int co_dim = dim_mesh; co_dim > 0; co_dim--) {
    base::dim_t codim_pt = dim_mesh - co_dim;  // co-dim of vertices
    // Loop over entities of the specified co-dimension
    for (const Entity* e : mesh.Entities(co_dim)) {
      const lf::base::RefEl e_refel = e->RefEl();
      const Eigen::MatrixXd& ref_el_corners(e_refel.NodeCoords());
      const Eigen::MatrixXd vertex_coords(
          e->Geometry()->Global(ref_el_corners));
      const double approx_area =
          (e->Geometry()->IntegrationElement(ref_el_corners))[0];
      // Visit all nodes
      auto sub_ents = e->SubEntities(codim_pt);
      for (base::sub_idx_t j = 0; j < e_refel.NumNodes(); ++j) {
        const Eigen::VectorXd node_coords(
            sub_ents[j]->Geometry()->Global(pt_ref_coord));
        // Check agreement of coordinates up to roundoff
        if ((vertex_coords.col(j) - node_coords).squaredNorm() >
            1.0E-8 * approx_area) {
          ret_vals.emplace_back(e_refel, mesh.Index(*e));

          SPDLOG_LOGGER_TRACE(WatertightLogger(),
                              "Node {} of {} ({}): position mismatch", j,
                              e_refel, mesh.Index(*e));

        }  // end geometry test
      }    // end loop over nodes
    }      // end loop over entities
  }        // end loop over co-dimensions
  if (!vertices_only) {
    // Check whether geometry of edges and cells match
    // ASSERT_MSG(false,
    //	    "Geometric compatibility test for edges and cells not yet
    // implemented");
  }

  return ret_vals;
}  // namespace lf::mesh::test_utils

}  // namespace lf::mesh::test_utils
