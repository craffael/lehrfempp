#include "test_meshes.h"
#include <iostream>

namespace lf::mesh::test_utils {

std::shared_ptr<lf::mesh::Mesh> GenerateHybrid2DTestMesh() {
  using coord_t = Eigen::Vector2d;
  using size_type = lf::mesh::Mesh::size_type;
  // Obtain mesh factory
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
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
      std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({3, 4, 1}),
      std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({4, 2, 1}),
      std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({4, 5, 2}),
      std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({5, 6, 2}),
      std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(),
      lf::base::ForwardRange<const size_type>({3, 1, 0, 9}),
      std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(),
      lf::base::ForwardRange<const size_type>({2, 6, 7, 0}),
      std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({0, 7, 8}),
      std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),
      lf::base::ForwardRange<const size_type>({0, 1, 2}),
      std::unique_ptr<lf::geometry::Geometry>(nullptr));
  // Inspect data
  mesh_factory_ptr->PrintLists(std::cout);
  return mesh_factory_ptr->Build();
}

}  // namespace lf::mesh::test_utils
