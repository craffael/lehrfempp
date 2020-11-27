/**
 * @file
 * @brief tests for the TorusMeshBuilder class
 * @author Anian Ruoss
 * @date   2018-10-22 16:17:17
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/torus_mesh_builder.h>
#include <lf/mesh/utils/utils.h>
#include <memory>
#include "lf/mesh/test_utils/check_entity_indexing.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"

namespace lf::mesh::test {

// Test for creation of tensor product grid on torus
TEST(lf_mesh_p, buildTorusMesh) {
  // Enable copious output
  utils::TorusMeshBuilder::Logger()->set_level(spdlog::level::trace);
  // Construct a tensor-product grid with 15 cells
  std::unique_ptr<hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_unique<hybrid2d::MeshFactory>(3);
  utils::TorusMeshBuilder builder(std::move(mesh_factory_ptr));
  // Set mesh parameters following the Builder pattern
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 4})
      .setNumXCells(3)
      .setNumYCells(5);
  auto mesh_p = builder.Build();

  EXPECT_TRUE(mesh_p) << "Oops! no mesh!";
  EXPECT_EQ(mesh_p->DimMesh(), 2) << "Mesh dimension != 2 !";
  EXPECT_EQ(mesh_p->DimWorld(), 3) << "World dimension must be 3";
  EXPECT_EQ(mesh_p->NumEntities(0), 15) << "Mesh should comprise 15 cells";
  EXPECT_EQ(mesh_p->NumEntities(1), 30) << "Mesh should comprise 30 edges";
  EXPECT_EQ(mesh_p->NumEntities(2), 15) << "Mesh should have 15 vertices";

  // check that every edge has two adjacent cells
  auto cells_per_edge = mesh::utils::CountNumSuperEntities(mesh_p, 1, 1);
  for (const auto* edge : mesh_p->Entities(1)) {
    EXPECT_EQ(cells_per_edge(*edge), 2) << "Edge should have 2 adjacent cells";
  }

  // check that every vertex has four adjacent edges and cells
  auto cells_per_vertex = mesh::utils::CountNumSuperEntities(mesh_p, 2, 2);
  auto edges_per_vertex = mesh::utils::CountNumSuperEntities(mesh_p, 2, 1);
  for (const auto* vertex : mesh_p->Entities(2)) {
    EXPECT_EQ(cells_per_vertex(*vertex), 4)
        << "Vertex should have 4 adjacent cells";
    EXPECT_EQ(edges_per_vertex(*vertex), 4)
        << "Vertex should have 4 adjacent edges";
  }

  test_utils::checkEntityIndexing(*mesh_p);
  test_utils::checkMeshCompleteness(*mesh_p);
}

}  // namespace lf::mesh::test
