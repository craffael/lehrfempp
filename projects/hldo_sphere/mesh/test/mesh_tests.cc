#include <gtest/gtest.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/mesh.h>
#include <sphere_triag_mesh_builder.h>

std::shared_ptr<lf::mesh::Mesh> buildSphere(
    const lf::base::size_type refinement_level, const double radius) {
  // prepare needed objects
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

  projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
      projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

  sphere.setRefinementLevel(refinement_level);
  sphere.setRadius(radius);

  return sphere.Build();
}

TEST(projects_hldo_sphere_mesh, mesh_properties) {
  const lf::base::size_type refinement_level = 5;
  const double radius = 2.0;

  // Optain sphere mesh
  const std::shared_ptr<lf::mesh::Mesh> mesh =
      buildSphere(refinement_level, radius);

  // Vertex count must be equal to cell count - 2
  EXPECT_EQ(mesh->NumEntities(0), 2 * (mesh->NumEntities(2) - 2))
      << "The cell count is not equal to the 2 * (vextex count - 2)";

  // TODO rewrite tests with closed form formulas

  // we create four new cells out of one in a refinement step and we start with
  // 8 cells
  lf::base::size_type c_0 = 8;

  // closed form of the recursion c_r = 4 * c_{r-1}
  double expected_cell_count = std::pow(4, refinement_level) * c_0;
  EXPECT_EQ(mesh->NumEntities(0), expected_cell_count)
      << "Expected number of cells doesn't match the outcome";

  // the expected number of edges is based on the recursion e_r = 2e_{r-1} +
  // 3c_{r-1} where r is the refinement_level
  lf::base::size_type e_0 = 12;

  // closed form of the recursion e_r = 2 * e_{r-1} + 3 * c_{r-1} where r is the
  // refinement_level
  double expected_edge_count =
      refinement_level == 0 ? e_0
                            : std::pow(2, refinement_level) * e_0 +
                                  3 * c_0 *
                                      (std::pow(2, 2 * refinement_level - 1) -
                                       std::pow(2, refinement_level - 1));
  EXPECT_EQ(mesh->NumEntities(1), expected_edge_count)
      << "Expected number of edges doesn't match the outcome";

  // closed form of the recursion v_r = v_{r-1} + e_{r-1} where r is the
  // refinement_level
  lf::base::size_type v_0 = 6;
  double expected_vertex_count =
      v_0 + e_0 * (std::pow(2, refinement_level) - 1) + c_0 +
      c_0 * std::pow(2, 2 * refinement_level - 1) -
      3 * c_0 * pow(2, refinement_level - 1);
  EXPECT_EQ(mesh->NumEntities(2), expected_vertex_count)
      << "Expected number of vertices doesn't match the outcome";
}

TEST(projects_hldo_sphere_mesh, vertex_adjacent_edges) {
  const lf::base::size_type refinement_level = 2;
  const double radius = 2.0;

  // Optain sphere mesh
  const std::shared_ptr<lf::mesh::Mesh> mesh =
      buildSphere(refinement_level, radius);

  // Create the coordinates of the six vertices with only four adjacent cells
  std::vector<Eigen::Vector3d> basic_vertices(6);
  basic_vertices[0] << radius * 1, 0, 0;
  basic_vertices[1] << radius * (-1), 0, 0;
  basic_vertices[2] << 0, radius * 1, 0;
  basic_vertices[3] << 0, radius * (-1), 0;
  basic_vertices[4] << 0, 0, radius * 1;
  basic_vertices[5] << 0, 0, radius * (-1);

  // create vector with the adjacent cells count for all vertices
  std::vector<int> adjacent_cell_count(mesh->NumEntities(2));

  // count adjacent cells for each edge
  for (const lf::mesh::Entity* c : mesh->Entities(0)) {
    for (const lf::mesh::Entity* v : c->SubEntities(2)) {
      adjacent_cell_count[mesh->Index(*v)]++;
    }
  }

  // count basic vertices
  const int basic_vertex_count = 0;

  // loop over all vertices and check if the vertex has 6 adjacent edges
  for (const lf::mesh::Entity* v : mesh->Entities(2)) {
    // check if the vertex is one of the four basic ones with four adjacent
    // cells
    const lf::geometry::Geometry* geometry = v->Geometry();
    Eigen::Vector3d v_coords = lf::geometry::Corners(*geometry);

    lf::base::size_type v_idx = mesh->Index(*v);

    if (adjacent_cell_count[v_idx] == 4) {
      double min_rel_error = 1;

      for (int i = 0; i < 6; ++i) {
        if ((basic_vertices[i] - v_coords).norm() < min_rel_error) {
          min_rel_error = (basic_vertices[i] - v_coords).norm();
        }
      }

      EXPECT_NEAR(0.0, min_rel_error, (v_coords.norm()) * 1e-12)
          << "The vertex at " << v_coords << " with global index " << v_idx
          << " only has 4 adjacent cells but is not at one of the basic vertex "
             "positions";

    } else if (adjacent_cell_count[v_idx] == 6) {
      // NOTHING TO DO HERE
    } else {
      EXPECT_TRUE(false) << "The vertex at " << v_coords
                         << " with global index " << v_idx
                         << " only has not 4 not 6 adjacent cells";
    }
  }
}
