/** @file mesh_entity_output.cc
 * Demo for outputting information about mesh entity elements
 */

#include <iostream>
//#include "lf/base/base.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/hybrid2dp/hybrid2dp.h"
#include "lf/mesh/utils/utils.h"
#include "lf/mesh/utils/print_info.h"

/*
#include <gtest/gtest.h>
#include <Eigen/Eigen>
#include "lf/mesh/test_utils/check_entity_indexing.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
*/


int main() {
  std::cout << "Output of information for mesh entity elements" << std::endl;

  // NOTE TO SELF: There are several test folders in lf/mesh:
  // - mesh/test
  // - mesh/hybrid2d/test
  // - mesh/hybrid2dp/test

  // Don't use hybrid2d!


  // Entities -------------------------------

  // !!! Do not use hybrid2d
  lf::mesh::hybrid2d::Entity<2> entity_2;
  lf::mesh::hybrid2d::Entity<1> entity_1;
  lf::mesh::hybrid2d::Entity<0> entity_0;

  std::cout << entity_2 << std::endl;
  std::cout << entity_2.RefEl() << std::endl; // kPoint
  std::cout << std::endl;

  std::cout << entity_0 << std::endl;
  std::cout << entity_0.RefEl() << std::endl; // kQuad
  std::cout << std::endl;

  std::cout << entity_1 << std::endl;
  std::cout << entity_1.RefEl() << std::endl; // kSegment
  // !!!

  lf::mesh::utils::writeTikZ(entity_2, "tikztest.txt");
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
  test.AddEntity(lf::base::RefEl::kQuad(), {0, 1, 2, 3}, std::make_unique<lf::geometry::QuadO1>(std::move(node_coord)));


  // explicitly add the right edge:
  node_coord = Eigen::MatrixXd(2, 2);
  node_coord << 1, 1, 0, 1;
  test.AddEntity(lf::base::RefEl::kSegment(), {1, 2}, std::make_unique<lf::geometry::SegmentO1>(node_coord));


  // build the mesh
  auto mesh = test.Build();

  int dim_mesh_MeshFactory = test.DimWorld();
  int dim_mesh_Mesh = mesh->DimWorld();

  // test is MeshFactory object
  std::cout << "\n" << test << std::endl;
  std::cout << "Dim of MeshFactory mesh: " << dim_mesh_MeshFactory << std::endl;

  // mesh is Mesh object
  std::cout << "\n" << mesh << std::endl;
  std::cout << "Dim of Mesh mesh: " << dim_mesh_Mesh << std::endl;



  lf::mesh::utils::PrintInfo(*mesh, std::cout); // Prints
  lf::mesh::utils::PrintInfo(entity_2, std::cout); // Prints


  return 0L;

}
