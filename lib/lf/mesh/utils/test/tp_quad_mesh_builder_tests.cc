/**
 * @file
 * @brief tests for the TPQuadMeshBuilder class
 * @author Raffael Casagrande
 * @date   2018-06-22 09:43:11
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/tp_quad_mesh_builder.h>
#include <lf/mesh/utils/utils.h>
#include <memory>
#include "lf/mesh/test_utils/check_entity_indexing.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"

namespace lf::mesh::test {

// Test for pointer-based implementation and creation of tensor product grid
TEST(lf_mesh_p, buildTPQuadMesh) {
  // Enable copious output
  utils::TPQuadMeshBuilder::logger->set_level(spdlog::level::trace);
  // Construct a tensor-product grid of the unit square
  // with 6 rectangular cells
  std::unique_ptr<hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_unique<hybrid2d::MeshFactory>(2);
  utils::TPQuadMeshBuilder builder(std::move(mesh_factory_ptr));
  // Set mesh parameters following the Builder pattern
  // Domain is the unit square
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
      .setTopRightCorner(Eigen::Vector2d{1, 1})
      .setNumXCells(3)
      .setNumYCells(2);
  auto mesh_p = builder.Build();

  EXPECT_TRUE(mesh_p) << "Oops! no mesh!";
  EXPECT_EQ(mesh_p->DimMesh(), 2) << "Mesh dimension != 2 !";
  EXPECT_EQ(mesh_p->DimWorld(), 2) << "World dimension must be 2";
  EXPECT_EQ(mesh_p->NumEntities(0), 6) << "Mesh should comprise 6 squares";
  EXPECT_EQ(mesh_p->NumEntities(1), 17) << "Mesh should comprise 17 edges";
  EXPECT_EQ(mesh_p->NumEntities(2), 12) << "Mesh should have 12 vertices";

  std::cout << "Checking entity indexing" << std::endl;
  test_utils::checkEntityIndexing(*mesh_p);
  std::cout << "Checking mesh completeness" << std::endl;
  test_utils::checkMeshCompleteness(*mesh_p);
  std::cout << "Printing mesh information" << std::endl;
  utils::PrintInfo(std::cout, *mesh_p);
}

}  // namespace lf::mesh::test
