/**
 * @file
 * @brief tests for the TPTriaMeshBuilder class
 * @author Raffael Casagrande
 * @date   2018-06-22 09:43:11
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/hybrid2dp/hybrid2dp.h>
#include <lf/mesh/tp_triag_mesh_builder.h>
#include <lf/mesh/utils/utils.h>
#include <memory>
#include "lf/mesh/test_utils/check_entity_indexing.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"

namespace lf::mesh::test {
// Test for index-based implementation
TEST(lf_mesh, buildStructuredMesh) {
  // Construct a structured mesh with 8 triangles
  hybrid2d::TPTriagMeshBuilder builder(
      std::make_shared<hybrid2d::MeshFactory>(2));
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNoXCells(2)
      .setNoYCells(2);
  auto mesh_p = builder.Build();

  EXPECT_TRUE(mesh_p) << "Oops! no mesh!";
  EXPECT_EQ(mesh_p->DimMesh(), 2) << "Mesh dimension != 2 !";
  EXPECT_EQ(mesh_p->DimWorld(), 2) << "World dimension must be 2";
  EXPECT_EQ(mesh_p->Size(0), 8) << "Mesh should comprise 8 triangles";
  EXPECT_EQ(mesh_p->Size(1), 16) << "Mesh should have 16 edges";
  EXPECT_EQ(mesh_p->Size(2), 9) << "Mesh should have 9 vertices";

  std::cout << "Checking entity indexing" << std::endl;
  test_utils::checkEntityIndexing(*mesh_p);
  std::cout << "Checking mesh completeness" << std::endl;
  test_utils::checkMeshCompleteness(*mesh_p);
}

// Test for pointer-based implementation
// Note the use of the namespace hybrid2dp
TEST(lf_mesh_p, buildStructuredMesh_p) {
  // Construct a structured mesh with 8 triangles
  std::shared_ptr<hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<hybrid2dp::MeshFactory>(2);
  hybrid2d::TPTriagMeshBuilder builder(mesh_factory_ptr);
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNoXCells(2)
      .setNoYCells(2);
  auto mesh_p = builder.Build();

  EXPECT_TRUE(mesh_p) << "Oops! no mesh!";
  EXPECT_EQ(mesh_p->DimMesh(), 2) << "Mesh dimension != 2 !";
  EXPECT_EQ(mesh_p->DimWorld(), 2) << "World dimension must be 2";
  EXPECT_EQ(mesh_p->Size(0), 8) << "Mesh should comprise 8 triangles";
  EXPECT_EQ(mesh_p->Size(1), 16) << "Mesh should have 16 edges";
  EXPECT_EQ(mesh_p->Size(2), 9) << "Mesh should have 9 vertices";

  std::cout << "Checking entity indexing" << std::endl;
  test_utils::checkEntityIndexing(*mesh_p);
  std::cout << "Checking mesh completeness" << std::endl;
  test_utils::checkMeshCompleteness(*mesh_p);
  std::cout << "Checking geometry compatibility: " << std::flush;
  lf::mesh::test_utils::watertight_mesh_ctrl = 100;
  auto fails = lf::mesh::test_utils::isWatertightMesh(*mesh_p, false);
  if (fails.empty()) {
    std::cout << "consistent!" << std::endl;
  } else {
    std::cout << "INCONSISTENT!" << std::endl;
    for (auto &geo_errs : fails) {
      std::cout << geo_errs.first.ToString() << "(" << geo_errs.second << ")"
                << std::endl;
    }
  }
  std::cout << "Writing MATLAB file" << std::endl;
  utils::writeMatlab(*mesh_p, "tp_triag_test.m");

  // Printing mesh information
  utils::PrintInfo(*mesh_p, std::cout);
}

// Test for pointer-based implementation
// and creation of tensor product grid
TEST(lf_mesh_p, buildTPQuadMesh) {
  // Enable copious output
  hybrid2d::TPQuadMeshBuilder::output_ctrl_ = 100;
  // Construct a tensor-product grid of the unit square
  // with 6 rectangular cells
  std::shared_ptr<hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<hybrid2dp::MeshFactory>(2);
  hybrid2d::TPQuadMeshBuilder builder(mesh_factory_ptr);
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNoXCells(3)
      .setNoYCells(2);
  auto mesh_p = builder.Build();

  EXPECT_TRUE(mesh_p) << "Oops! no mesh!";
  EXPECT_EQ(mesh_p->DimMesh(), 2) << "Mesh dimension != 2 !";
  EXPECT_EQ(mesh_p->DimWorld(), 2) << "World dimension must be 2";
  EXPECT_EQ(mesh_p->Size(0), 6) << "Mesh should comprise 6 squares";
  EXPECT_EQ(mesh_p->Size(1), 17) << "Mesh should comprise 17 edges";
  EXPECT_EQ(mesh_p->Size(2), 12) << "Mesh should have 12 vertices";

  std::cout << "Checking entity indexing" << std::endl;
  test_utils::checkEntityIndexing(*mesh_p);
  std::cout << "Checking mesh completeness" << std::endl;
  test_utils::checkMeshCompleteness(*mesh_p);
  std::cout << "Printing mesh information" << std::endl;
  utils::PrintInfo(*mesh_p, std::cout);
}

// Test for creation of tensor product grid on torus
TEST(lf_mesh_p, buildTorusMesh) {
  // Enable copious output
  hybrid2d::TPQuadMeshBuilder::output_ctrl_ = 100;
  // Construct a tensor-product grid with 15 cells
  std::shared_ptr<hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<hybrid2dp::MeshFactory>(3);
  hybrid2d::TorusMeshBuilder builder(mesh_factory_ptr);
  // Set mesh parameters following the Builder pattern
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 4})
      .setNoXCells(3)
      .setNoYCells(5);
  auto mesh_p = builder.Build();

  EXPECT_TRUE(mesh_p) << "Oops! no mesh!";
  EXPECT_EQ(mesh_p->DimMesh(), 2) << "Mesh dimension != 2 !";
  EXPECT_EQ(mesh_p->DimWorld(), 3) << "World dimension must be 3";
  EXPECT_EQ(mesh_p->Size(0), 15) << "Mesh should comprise 15 cells";
  EXPECT_EQ(mesh_p->Size(1), 30) << "Mesh should comprise 30 edges";
  EXPECT_EQ(mesh_p->Size(2), 15) << "Mesh should have 15 vertices";

  // check that every edge has two adjacent cells
  auto cells_per_edge = mesh::utils::countNoSuperEntities(mesh_p, 1, 1);
  for (const auto &edge : mesh_p->Entities(1)) {
    EXPECT_EQ(cells_per_edge(edge), 2) << "Edge should have 2 adjacent cells";
  }

  // check that every vertex has four adjacent edges and cells
  auto cells_per_vertex = mesh::utils::countNoSuperEntities(mesh_p, 2, 2);
  auto edges_per_vertex = mesh::utils::countNoSuperEntities(mesh_p, 2, 1);
  for (const auto &vertex : mesh_p->Entities(2)) {
    EXPECT_EQ(cells_per_vertex(vertex), 4)
        << "Vertex should have 4 adjacent cells";
    EXPECT_EQ(edges_per_vertex(vertex), 4)
        << "Vertex should have 4 adjacent edges";
  }

  test_utils::checkEntityIndexing(*mesh_p);
  test_utils::checkMeshCompleteness(*mesh_p);
}

}  // namespace lf::mesh::test
