/**
 * @file
 * @brief Tests for the lf::mesh::hybrid2dp::MeshFactory
 * @author Raffael Casagrande
 * @date   2018-07-01 01:15:55
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/mesh/hybrid2dp/hybrid2dp.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <Eigen/Eigen>
#include "lf/mesh/test_utils/check_entity_indexing.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"

namespace lf::mesh::hybrid2dp::test {

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
  utils::PrintInfo(*mesh_p, std::cout);
}

}  // namespace lf::mesh::hybrid2dp::test
