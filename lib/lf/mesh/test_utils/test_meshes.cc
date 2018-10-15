/** @file
 * @brief Functions generating 2D hybrid meshes for testing purposes
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#include "test_meshes.h"
#include <iostream>

namespace lf::mesh::test_utils {

std::shared_ptr<lf::mesh::Mesh> GenerateHybrid2DTestMesh(int selector) {
  using coord_t = Eigen::Vector2d;
  using size_type = lf::mesh::Mesh::size_type;
  using quad_coord_t = Eigen::Matrix<double, 2, 4>;
  using tria_coord_t = Eigen::Matrix<double, 2, 3>;

  // Obtain mesh factory
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);

  switch (selector) {
    case 0: {
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
      break;
    }
    case 1: {
      // Set coordinates of nodes
      mesh_factory_ptr->AddPoint(coord_t({0, 0}));  // 0
      mesh_factory_ptr->AddPoint(coord_t({2, 0}));  // 1
      mesh_factory_ptr->AddPoint(coord_t({4, 0}));  // 2
      mesh_factory_ptr->AddPoint(coord_t({0, 2}));  // 3
      mesh_factory_ptr->AddPoint(coord_t({0, 3}));  // 4
      mesh_factory_ptr->AddPoint(coord_t({2, 4}));  // 5
      mesh_factory_ptr->AddPoint(coord_t({0, 4}));  // 6
      mesh_factory_ptr->AddPoint(coord_t({4, 4}));  // 7
      mesh_factory_ptr->AddPoint(coord_t({4, 2}));  // 8
      mesh_factory_ptr->AddPoint(coord_t({1, 1}));  // 9
      mesh_factory_ptr->AddPoint(coord_t({3, 1}));  // 10
      mesh_factory_ptr->AddPoint(coord_t({3, 2}));  // 11
      mesh_factory_ptr->AddPoint(coord_t({2, 3}));  // 12
      mesh_factory_ptr->AddPoint(coord_t({1, 2}));  // 13
      mesh_factory_ptr->AddPoint(coord_t({2, 2}));  // 14

      quad_coord_t quad_coord(2, 4);
      // First cell: a parallelogram
      quad_coord << 0.0, 2.0, 3.0, 1.0, 0.0, 0.0, 1.0, 1.0;
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({0, 1, 10, 9}),
          std::make_unique<lf::geometry::Parallelogram>(quad_coord));
      // Second cell: a triangle
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({1, 2, 10}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      // Third cell: a quadrilateral
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({2, 8, 11, 10}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      // Fourth cell: a triangle
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({11, 8, 7}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      // Fifth cell: a quadrilateral
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({5, 12, 11, 7}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      // Sixth cell: A parallelogram
      quad_coord << 0, 2, 2, 0, 3, 3, 4, 4;
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({4, 12, 5, 6}),
          std::make_unique<lf::geometry::Parallelogram>(quad_coord));
      // Seventh cell: a quadrilateral
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({4, 3, 13, 12}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      // Eigth cell: a quadrilateral
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({3, 0, 9, 13}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      // Ninth cell: a triangle
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({9, 10, 14}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      // Tenth cell: a parallelogram
      quad_coord << 2, 3, 3, 2, 2, 1, 2, 3;
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({14, 10, 11, 12}),
          std::make_unique<lf::geometry::Parallelogram>(quad_coord));
      // 11th cell: a parallelogram
      quad_coord << 1, 1, 2, 2, 2, 1, 2, 3;
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({13, 9, 14, 12}),
          std::make_unique<lf::geometry::Parallelogram>(quad_coord));
      break;
    }
    default: {
      LF_VERIFY_MSG(false, "Illegal selector for test meshes");
      break;
    }
  }  // end switch
  // Inspect data
  mesh_factory_ptr->PrintLists(std::cout);
  return mesh_factory_ptr->Build();
}

}  // namespace lf::mesh::test_utils
