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
#include <lf/mesh/mesh_utils/write_matlab.h>
#include <lf/mesh/tp_triag_mesh_builder.h>
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
  EXPECT_EQ(mesh_p->DimWorld(), 2) << "Wolrd dimension must be 2";
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
  EXPECT_EQ(mesh_p->DimWorld(), 2) << "Wolrd dimension must be 2";
  EXPECT_EQ(mesh_p->Size(0), 8) << "Mesh should comprise 8 triangles";
  EXPECT_EQ(mesh_p->Size(1), 16) << "Mesh should have 16 edges";
  EXPECT_EQ(mesh_p->Size(2), 9) << "Mesh should have 9 vertices";

  std::cout << "Checking entity indexing" << std::endl;
  test_utils::checkEntityIndexing(*mesh_p);
  std::cout << "Checking mesh completeness" << std::endl;
  test_utils::checkMeshCompleteness(*mesh_p);
  std::cout << "Writing MATLAB file" << std::endl;
  mesh_utils::writeMatlab(*mesh_p, "tp_triag_test.m");
}

// Test for generating a mesh with reconstruction of edge information
TEST(lf_edge_create, MeshFactory_p) {
  using coord_t = Eigen::Vector2d;
  using size_type = mesh::Mesh::size_type;
  // Obtain mesh factory
  std::shared_ptr<hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<hybrid2dp::MeshFactory>(2);
  // Setting point coordinate
  mesh_factory_ptr->AddPoint(coord_t({1.5, 2}));
  mesh_factory_ptr->AddPoint(coord_t({1, 1}));
  mesh_factory_ptr->AddPoint(coord_t({2, 1}));
  mesh_factory_ptr->AddPoint(coord_t({0, 0}));
  mesh_factory_ptr->AddPoint(coord_t({1.5, 0}));
  mesh_factory_ptr->AddPoint(coord_t({3, 0}));
  mesh_factory_ptr->AddPoint(coord_t({3, 2}));
  mesh_factory_ptr->AddPoint(coord_t({3, 3}));
  mesh_factory_ptr->AddPoint(coord_t({0, 3}));
  mesh_factory_ptr->AddPoint(coord_t({0, 1}));

  // Setting vertices of cells but not their geometry
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({0, 8, 9}),
      std::unique_ptr<geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({3, 4, 1}),
      std::unique_ptr<geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({4, 2, 1}),
      std::unique_ptr<geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({4, 5, 2}),
      std::unique_ptr<geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({5, 6, 2}),
      std::unique_ptr<geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(),
      lf::base::ForwardRange<const size_type>({3, 1, 0, 9}),
      std::unique_ptr<geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(),
      lf::base::ForwardRange<const size_type>({2, 6, 7, 0}),
      std::unique_ptr<geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({0, 7, 8}),
      std::unique_ptr<geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({0, 1, 2}),
      std::unique_ptr<geometry::Geometry>(nullptr));
  // Inspect data
  std::cout << "**************************************************"
            << std::endl;
  std::cout << "Internal data of MeshFactory" << std::endl;
  mesh_factory_ptr->PrintLists(std::cout);
  std::cout << "**************************************************"
            << std::endl;

  // Building the mesh
  std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
            << std::endl;
  std::cout << "&&&& Building mesh &&&&" << std::endl;
  auto mesh_p = mesh_factory_ptr->Build();
  std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
            << std::endl;

  std::cout << "Checking entity indexing" << std::endl;
  test_utils::checkEntityIndexing(*mesh_p);
  std::cout << "Checking mesh completeness" << std::endl;
  test_utils::checkMeshCompleteness(*mesh_p);

  // Printing mesh information
  printInfo(*mesh_p, std::cout);
}

}  // namespace lf::mesh::test
