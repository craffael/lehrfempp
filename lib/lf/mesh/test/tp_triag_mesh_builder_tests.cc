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

  EXPECT_NE(mesh_p, nullptr) << "Oops! no mesh!";
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

  EXPECT_NE(mesh_p, nullptr) << "Oops! no mesh!";
  EXPECT_EQ(mesh_p->DimMesh(), 2) << "Mesh dimension != 2 !";
  EXPECT_EQ(mesh_p->DimWorld(), 2) << "World dimension must be 2";
  EXPECT_EQ(mesh_p->Size(0), 8) << "Mesh should comprise 8 triangles";
  EXPECT_EQ(mesh_p->Size(1), 16) << "Mesh should have 16 edges";
  EXPECT_EQ(mesh_p->Size(2), 9) << "Mesh should have 9 vertices";

  std::cout << "Checking entity indexing" << std::endl;
  test_utils::checkEntityIndexing(*mesh_p);
  std::cout << "Checking mesh completeness" << std::endl;
  test_utils::checkMeshCompleteness(*mesh_p);
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

  EXPECT_NE(mesh_p, nullptr) << "Oops! no mesh!";
  EXPECT_EQ(mesh_p->DimMesh(), 2) << "Mesh dimension != 2 !";
  EXPECT_EQ(mesh_p->DimWorld(), 2) << "Wolrd dimension must be 2";
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

}  // namespace lf::mesh::test
