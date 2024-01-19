#include <gtest/gtest.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <sphere_triag_mesh_builder.h>

/**
 * Method to create a mesh with desired parameters
 */
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

/**
 * Check the mesh with refinement level 5
 * on the number of cells, edges and vertices.
 */
TEST(projects_hldo_sphere_mesh, mesh_properties) {
  const lf::base::size_type refinement_level = 5;
  const double radius = 2.0;

  // Optain sphere mesh
  const std::shared_ptr<lf::mesh::Mesh> mesh =
      buildSphere(refinement_level, radius);

  // Vertex count must be equal to cell count - 2
  EXPECT_EQ(mesh->NumEntities(0), 2 * (mesh->NumEntities(2) - 2))
      << "The cell count is not equal to the 2 * (vextex count - 2)";

  // see thesis appendix C for reference
  double expected_cell_count = std::pow(2, 2 * refinement_level + 3);
  EXPECT_EQ(mesh->NumEntities(0), expected_cell_count)
      << "Expected number of cells doesn't match the outcome";

  // see thesis appendix C for reference
  double expected_edge_count = 3 * pow(2., 2 * refinement_level + 2);
  EXPECT_EQ(mesh->NumEntities(1), expected_edge_count)
      << "Expected number of edges doesn't match the outcome";

  // see thesis appendix C for reference
  double expected_vertex_count = pow(2., 2 * refinement_level + 2) + 2;
  EXPECT_EQ(mesh->NumEntities(2), expected_vertex_count)
      << "Expected number of vertices doesn't match the outcome";
}

/**
 *
 * Checks if the vertices have the right number of ajacent cells
 *
 */
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
          << " only has 4 adjacent cells but is not one of the basic vertex "
             "positions";

    } else if (adjacent_cell_count[v_idx] == 6) {
      // NOTHING TO DO HERE
    } else {
      EXPECT_TRUE(false) << "The vertex at " << v_coords
                         << " with global index " << v_idx
                         << " only doesn't have 6 adjacent cells";
    }
  }
}

/**
 *
 * Checks if all edges are locally oriented in the same direction.
 * making sure that the sums of relative orientations are zero.
 *
 */
TEST(projects_hldo_sphere_mesh, local_edge_directions) {
  const lf::base::size_type refinement_level = 5;
  const double radius = 1.0;

  for (int l = 0; l < refinement_level; l++) {
    // Optain sphere mesh
    const std::shared_ptr<lf::mesh::Mesh> mesh =
        buildSphere(refinement_level, radius);

    // compute for all edges if the local direction add up to 0
    //
    // loop over all cells and add the values for the edge orientations (they
    // sould cancel out)
    // initialized to 0
    lf::mesh::utils::CodimMeshDataSet<double> edgeSums(mesh, 1, 0);
    for (const lf::mesh::Entity* e : mesh->Entities(0)) {
      auto edgeOrientations = e->RelativeOrientations();
      Eigen::VectorXd s(3);
      auto localEdges = e->SubEntities(1);
      for (int i = 0; i < 3; i++) {
        edgeSums(*localEdges[i]) += lf::mesh::to_sign(edgeOrientations[i]);
      }
    }

    // loop over edges and check
    for (const lf::mesh::Entity* e : mesh->Entities(1)) {
      int index = mesh->Index(*e);
      EXPECT_EQ(0, edgeSums(*e))
          << "The local edge orientation of the edge with index " << index
          << " do not add up to 0";
    }
  }
}
