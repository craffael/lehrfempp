/** @file
 * @brief Functions generating 2D hybrid meshes for testing purposes
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#include "test_meshes.h"
#include <iostream>

namespace lf::mesh::test_utils {

std::shared_ptr<lf::mesh::Mesh> GenerateHybrid2DTestMesh(int selector,
                                                         double scale) {
  using coord_t = Eigen::Vector2d;
  using size_type = lf::mesh::Mesh::size_type;
  using quad_coord_t = Eigen::Matrix<double, 2, 4>;
  using tria_coord_t = Eigen::Matrix<double, 2, 3>;

  // Obtain mesh factory
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

  switch (selector) {
    case 0: {
      // Setting point coordinate
      mesh_factory_ptr->AddPoint(coord_t({1.5 * scale, 2 * scale}));
      mesh_factory_ptr->AddPoint(coord_t({1 * scale, 1 * scale}));
      mesh_factory_ptr->AddPoint(coord_t({2 * scale, 1 * scale}));
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 0 * scale}));
      mesh_factory_ptr->AddPoint(coord_t({1.5 * scale, 0 * scale}));
      mesh_factory_ptr->AddPoint(coord_t({3 * scale, 0 * scale}));
      mesh_factory_ptr->AddPoint(coord_t({3 * scale, 2 * scale}));
      mesh_factory_ptr->AddPoint(coord_t({3 * scale, 3 * scale}));
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 3 * scale}));
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 1 * scale}));

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
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 0 * scale}));  // 0
      mesh_factory_ptr->AddPoint(coord_t({2 * scale, 0 * scale}));  // 1
      mesh_factory_ptr->AddPoint(coord_t({4 * scale, 0 * scale}));  // 2
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 2 * scale}));  // 3
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 3 * scale}));  // 4
      mesh_factory_ptr->AddPoint(coord_t({2 * scale, 4 * scale}));  // 5
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 4 * scale}));  // 6
      mesh_factory_ptr->AddPoint(coord_t({4 * scale, 4 * scale}));  // 7
      mesh_factory_ptr->AddPoint(coord_t({4 * scale, 2 * scale}));  // 8
      mesh_factory_ptr->AddPoint(coord_t({1 * scale, 1 * scale}));  // 9
      mesh_factory_ptr->AddPoint(coord_t({3 * scale, 1 * scale}));  // 10
      mesh_factory_ptr->AddPoint(coord_t({3 * scale, 2 * scale}));  // 11
      mesh_factory_ptr->AddPoint(coord_t({2 * scale, 3 * scale}));  // 12
      mesh_factory_ptr->AddPoint(coord_t({1 * scale, 2 * scale}));  // 13
      mesh_factory_ptr->AddPoint(coord_t({2 * scale, 2 * scale}));  // 14

      quad_coord_t quad_coord(2, 4);
      // First cell: a parallelogram
      quad_coord << 0.0, 2.0, 3.0, 1.0, 0.0, 0.0, 1.0, 1.0;
      quad_coord *= scale;
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
      quad_coord *= scale;
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
      quad_coord *= scale;
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({14, 10, 11, 12}),
          std::make_unique<lf::geometry::Parallelogram>(quad_coord));
      // 11th cell: a parallelogram
      quad_coord << 1, 1, 2, 2, 2, 1, 2, 3;
      quad_coord *= scale;
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({13, 9, 14, 12}),
          std::make_unique<lf::geometry::Parallelogram>(quad_coord));
      break;
    }
    case 2: {
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 0 * scale}));      // 0
      mesh_factory_ptr->AddPoint(coord_t({1 * scale, 0 * scale}));      // 1
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 1 * scale}));      // 2
      mesh_factory_ptr->AddPoint(coord_t({1 * scale, 1 * scale}));      // 3
      mesh_factory_ptr->AddPoint(coord_t({1.5 * scale, 0.5 * scale}));  // 4
      quad_coord_t quad_coord(2, 4);
      // First cell: the unit square
      quad_coord << 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0;
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kQuad(),
          lf::base::ForwardRange<const size_type>({0, 1, 3, 2}),
          std::make_unique<lf::geometry::Parallelogram>(quad_coord));
      // Second cell: an affine triangle
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({1, 3, 4}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      break;
    }
    case 3: {
      // Triangular mesh
      // Set coordinates of nodes
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 0 * scale}));    // 0
      mesh_factory_ptr->AddPoint(coord_t({1.5 * scale, 0 * scale}));  // 1
      mesh_factory_ptr->AddPoint(coord_t({3 * scale, 0 * scale}));    // 2
      mesh_factory_ptr->AddPoint(coord_t({1 * scale, 1 * scale}));    // 3
      mesh_factory_ptr->AddPoint(coord_t({2 * scale, 1 * scale}));    // 4
      mesh_factory_ptr->AddPoint(coord_t({3 * scale, 1 * scale}));    // 5
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 1.5 * scale}));  // 6
      mesh_factory_ptr->AddPoint(coord_t({2 * scale, 1.5 * scale}));  // 7
      mesh_factory_ptr->AddPoint(coord_t({1 * scale, 2 * scale}));    // 8
      mesh_factory_ptr->AddPoint(coord_t({3 * scale, 2 * scale}));    // 9
      mesh_factory_ptr->AddPoint(coord_t({0 * scale, 3 * scale}));    // 10
      mesh_factory_ptr->AddPoint(coord_t({1.5 * scale, 3 * scale}));  // 11
      mesh_factory_ptr->AddPoint(coord_t({3 * scale, 3 * scale}));    // 12

      // Setting vertices of triangular cells but not their geometry
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({0, 1, 3}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({1, 2, 4}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({0, 3, 6}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({1, 3, 4}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({4, 2, 5}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({3, 4, 7}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({4, 5, 7}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({3, 6, 8}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({7, 8, 3}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({5, 7, 9}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({6, 8, 10}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({8, 7, 11}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({7, 9, 12}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({10, 8, 11}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      mesh_factory_ptr->AddEntity(
          lf::base::RefEl::kTria(),
          lf::base::ForwardRange<const size_type>({7, 12, 11}),
          std::unique_ptr<lf::geometry::Geometry>(nullptr));
      break;
    }
    case 4: {
      // Triangular tensor product mesh
      // Construct a structured mesh with 18 triangles
      hybrid2d::TPTriagMeshBuilder builder(mesh_factory_ptr);
      // Set mesh parameters following the Builder pattern
      // Domain is the unit square
      builder.setBottomLeftCorner(Eigen::Vector2d{0 * scale, 0 * scale})
          .setTopRightCorner(Eigen::Vector2d{1 * scale, 1 * scale})
          .setNoXCells(3)
          .setNoYCells(3);
      return builder.Build();
      break;
    }
    case 5: {
      // Hybrid mesh of [0,3]^2 with only affine cells
      // Set coordinates of nodes
      std::array<std::array<double, 2>, 16> node_coord{
          std::array<double, 2>({0, 0}), std::array<double, 2>({1, 0}),
          std::array<double, 2>({2, 0}), std::array<double, 2>({3, 0}),
          std::array<double, 2>({0, 1}), std::array<double, 2>({1, 1}),
          std::array<double, 2>({2, 1}), std::array<double, 2>({3, 1}),
          std::array<double, 2>({0, 2}), std::array<double, 2>({1, 2}),
          std::array<double, 2>({2, 2}), std::array<double, 2>({3, 2}),
          std::array<double, 2>({0, 3}), std::array<double, 2>({1, 3}),
          std::array<double, 2>({2, 3}), std::array<double, 2>({3, 3})};

      // Create nodes
      for (const auto &node : node_coord) {
        mesh_factory_ptr->AddPoint(coord_t({node[0] * scale, node[1] * scale}));
      }

      // Specify triangles
      std::array<std::array<size_type, 3>, 4> tria_nodes{
          std::array<size_type, 3>({0, 4, 5}),
          std::array<size_type, 3>({2, 3, 7}),
          std::array<size_type, 3>({8, 9, 12}),
          std::array<size_type, 3>({11, 15, 14}),
      };

      // Specify parallelograms
      std::array<std::array<size_type, 4>, 7> quad_nodes{
          std::array<size_type, 4>({0, 1, 6, 5}),
          std::array<size_type, 4>({1, 2, 7, 6}),
          std::array<size_type, 4>({5, 6, 10, 9}),
          std::array<size_type, 4>({6, 7, 11, 10}),
          std::array<size_type, 4>({9, 10, 13, 12}),
          std::array<size_type, 4>({10, 11, 14, 13}),
          std::array<size_type, 4>({4, 5, 9, 8})};

      // generate triangles
      for (const auto &node : tria_nodes) {
        mesh_factory_ptr->AddEntity(
            lf::base::RefEl::kTria(),
            lf::base::ForwardRange<const size_type>(
                {node[0], node[1], node[2]}),
            std::unique_ptr<lf::geometry::Geometry>(nullptr));
      }

      // generate Parallelograms
      for (const auto &node : quad_nodes) {
        quad_coord_t quad_coord(2, 4);
        for (int n_pt = 0; n_pt < 4; ++n_pt) {
          quad_coord(0, n_pt) = node_coord[node[n_pt]][0];
          quad_coord(1, n_pt) = node_coord[node[n_pt]][1];
        }
        mesh_factory_ptr->AddEntity(
            lf::base::RefEl::kQuad(),
            lf::base::ForwardRange<const size_type>(
                {node[0], node[1], node[2], node[3]}),
            std::make_unique<lf::geometry::Parallelogram>(quad_coord));
      }
      break;
    }
    case 6: {
      // Hybrid mesh for triangle domain with vertices [0,0],[1,0],[0.2,1]
      // Set coordinates of nodes
      std::array<std::array<double, 2>, 8> node_coord{
          std::array<double, 2>({0, 0}),
          std::array<double, 2>({1, 0}),
          std::array<double, 2>({0.2, 1.0}),
          std::array<double, 2>({0.5, 0.0}),
          std::array<double, 2>({0.6, 0.5}),
          std::array<double, 2>({0.1, 0.5}),
          std::array<double, 2>({0.15, 0.75}),
          std::array<double, 2>({0.4, 0.75})};

      // Specify triangles (two)
      std::array<std::array<size_type, 3>, 2> tria_nodes{
          std::array<size_type, 3>({3, 1, 4}),
          std::array<size_type, 3>({7, 6, 2})};

      // Specify general quadrilateral
      std::array<std::array<size_type, 4>, 1> quad_nodes{
          std::array<size_type, 4>({5, 4, 7, 6})};

      // Specify parallelogram
      std::array<std::array<size_type, 4>, 1> parg_nodes{
          std::array<size_type, 4>({0, 3, 4, 5})};

      // Create nodes
      for (const auto &node : node_coord) {
        mesh_factory_ptr->AddPoint(coord_t({node[0] * scale, node[1] * scale}));
      }

      // generate triangles
      for (const auto &node : tria_nodes) {
        mesh_factory_ptr->AddEntity(
            lf::base::RefEl::kTria(),
            lf::base::ForwardRange<const size_type>(
                {node[0], node[1], node[2]}),
            std::unique_ptr<lf::geometry::Geometry>(nullptr));
      }

      // generate quadrilaterals
      for (const auto &node : quad_nodes) {
        quad_coord_t quad_coord(2, 4);
        for (int n_pt = 0; n_pt < 4; ++n_pt) {
          quad_coord(0, n_pt) = node_coord[node[n_pt]][0];
          quad_coord(1, n_pt) = node_coord[node[n_pt]][1];
        }
        mesh_factory_ptr->AddEntity(
            lf::base::RefEl::kQuad(),
            lf::base::ForwardRange<const size_type>(
                {node[0], node[1], node[2], node[3]}),
            std::make_unique<lf::geometry::QuadO1>(quad_coord));
      }

      // generate Parallelograms
      for (const auto &node : parg_nodes) {
        quad_coord_t quad_coord(2, 4);
        for (int n_pt = 0; n_pt < 4; ++n_pt) {
          quad_coord(0, n_pt) = node_coord[node[n_pt]][0];
          quad_coord(1, n_pt) = node_coord[node[n_pt]][1];
        }
        mesh_factory_ptr->AddEntity(
            lf::base::RefEl::kQuad(),
            lf::base::ForwardRange<const size_type>(
                {node[0], node[1], node[2], node[3]}),
            std::make_unique<lf::geometry::Parallelogram>(quad_coord));
      }

      break;
    }
    default: {
      LF_VERIFY_MSG(false, "Illegal selector for test meshes");
      break;
    }
  }  // end switch
  // Optional: Inspect data
  // mesh_factory_ptr->PrintLists(std::cout);
  return mesh_factory_ptr->Build();
}  // namespace lf::mesh::test_utils

}  // namespace lf::mesh::test_utils
