/** @file mesh_entity_output.cc
 * Demo for outputting information about mesh entity elements
 */

#include <iostream>

#include "lf/mesh/mesh.h"
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/hybrid2dp/hybrid2dp.h"
#include "lf/mesh/utils/utils.h"
#include "lf/mesh/utils/print_info.h"

int main() {
  using namespace lf::mesh::utils;
  
  std::cout << "Output of information for mesh entity elements" << std::endl;

  // NOTE TO SELF: There are several test folders in lf/mesh:
  // - mesh/test
  // - mesh/hybrid2d/test
  // - mesh/hybrid2dp/test

  // Don't use hybrid2d!


  // Entities -------------------------------

  //lf::mesh::utils::writeTikZ(lf::mesh::hybrid2dp::Triangle(), "tikztest.txt");
  // file is saved in /Build/examples/mesh/


  // Test ----------------------------------
  lf::mesh::hybrid2dp::MeshFactory test(2); // MeshFactory object

  // add nodes
  test.AddPoint(Eigen::Vector2d(0,0));
  test.AddPoint(Eigen::Vector2d(0,1));
  test.AddPoint(Eigen::Vector2d(1,1));
  test.AddPoint(Eigen::Vector2d(1,0));

  // add an element
  Eigen::MatrixXd node_coord(2, 4);
  node_coord << 0, 1, 1, 0, 0, 0, 1, 1;
  test.AddEntity(lf::base::RefEl::kQuad(), {0, 1, 2, 3},
		 std::make_unique<lf::geometry::QuadO1>(std::move(node_coord)));


  // explicitly add the right edge:
  node_coord = Eigen::MatrixXd(2, 2);
  node_coord << 1, 1, 0, 1;
  test.AddEntity(lf::base::RefEl::kSegment(), {1, 2},
		 std::make_unique<lf::geometry::SegmentO1>(node_coord));

  // build the mesh and retrieve a pointer
  auto mesh = test.Build(); // mesh is Mesh object

  int dim_mesh_MeshFactory = test.DimWorld(); // test is MeshFactory object
  int dim_mesh_Mesh = mesh->DimWorld();

  // Output information on mesh
  std::cout << "##### Mesh information ######" << std::endl;
  PrintInfo(*mesh, std::cout);
  std::cout << "#####                   #####" << std::endl;

  lf::mesh::Entity::output_ctrl_ = 100; 


  std::cout << "****** Output of mesh entities *******" << std::endl;
  // Loop over entities and print associated information
  for (lf::base::dim_t codim=0; codim <= 2; ++codim) {
    std::cout << "******* Entities of codimension " << codim << " ******* " << std::endl;
    for (const lf::mesh::Entity &entity : mesh->Entities(codim)) {
      std::cout << entity << std::endl;
    }
  }

  return 0L;
}
