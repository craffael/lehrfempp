/** @file mesh_entity_output.cc
 * Demo for outputting information about mesh entity elements
 */

#include <algorithm>  // for std::find
#include <iostream>
#include <vector>

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

int main() {
  using lf::io::TikzOutputCtrl;
  std::cout << "Output of information for mesh entity elements" << std::endl;

  // Build mesh ----------------------------------
  lf::mesh::hybrid2d::MeshFactory test(2);  // MeshFactory object

  // add nodes
  test.AddPoint(Eigen::Vector2d(0, 0));
  test.AddPoint(Eigen::Vector2d(0, 1));
  test.AddPoint(Eigen::Vector2d(1, 1));
  test.AddPoint(Eigen::Vector2d(1, 0));

  // add an element
  Eigen::MatrixXd node_coord(2, 4);
  node_coord << 0, 1, 1, 0, 0, 0, 1, 1;
  test.AddEntity(lf::base::RefEl::kQuad(),
                 std::array<lf::base::size_type, 4>{0, 1, 2, 3},
                 std::make_unique<lf::geometry::QuadO1>(std::move(node_coord)));

  // explicitly add the right edge:
  node_coord = Eigen::MatrixXd(2, 2);
  node_coord << 1, 1, 0, 1;
  test.AddEntity(lf::base::RefEl::kSegment(),
                 std::array<lf::base::size_type, 2>{1, 2},
                 std::make_unique<lf::geometry::SegmentO1>(node_coord));

  // build the mesh and retrieve a pointer
  auto mesh =
      test.Build();  // mesh is Mesh object // test is MeshFactory object

  // Output information on mesh
  std::cout << "##### Mesh information ######" << std::endl;
  lf::mesh::utils::PrintInfo(std::cout, *mesh);
  std::cout << "#####                   #####" << std::endl;

  std::cout << "****** Output of mesh entities *******" << std::endl;
  // Loop over entities and print associated information
  for (lf::base::dim_t codim = 0; codim <= 2; ++codim) {
    std::cout << "******* Entities of codimension " << static_cast<int>(codim)
              << " ******* " << std::endl;
    for (const lf::mesh::Entity* entity : mesh->Entities(codim)) {
      lf::mesh::utils::PrintInfo(std::cout, *entity, 100);
    }
  }

  // Another example
  // Test mesh ---------------------------------------------------
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // write_tikz version 2.0 --------------------------------------

  auto desiredEntities = [&](const lf::mesh::Entity& entity) -> bool {
    // If list contains entity, return true
    // return std::find(entitiesToPrint.begin(), entitiesToPrint.end(), entity)
    // != entitiesToPrint.end();

    // Nodes only
    return (entity.RefEl() == lf::base::RefEl::kPoint());
    // Segments only
    // return (entity.RefEl() == lf::base::RefEl::kSegment());
    // return true;
  };  // auto

  lf::io::writeTikZ(*mesh_p, "tikz_selector.txt", desiredEntities);

  lf::io::writeTikZ(
      *mesh_p, "tikz_mesh_test.txt", desiredEntities,
      TikzOutputCtrl::RenderCells | TikzOutputCtrl::VerticeNumbering |
          TikzOutputCtrl::EdgeNumbering | TikzOutputCtrl::CellNumbering |
          TikzOutputCtrl::NodeNumbering | TikzOutputCtrl::ArrowTips);

  return 0L;
}
