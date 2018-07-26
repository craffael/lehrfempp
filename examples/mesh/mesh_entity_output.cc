/** @file mesh_entity_output.cc
 * Demo for outputting information about mesh entity elements
 */

#include <iostream>
//#include "lf/base/base.h"
#include "lf/mesh/mesh.h"
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/hybrid2dp/hybrid2dp.h"
//#include "lf/mesh/utils/utils.h"
//#include "lf/mesh/utils/print_info.h"

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
  test.AddEntity(lf::base::RefEl::kQuad(), {0, 1, 2, 3}, std::make_unique<lf::geometry::QuadO1>(std::move(node_coord)));


  // explicitly add the right edge:
  node_coord = Eigen::MatrixXd(2, 2);
  node_coord << 1, 1, 0, 1;
  test.AddEntity(lf::base::RefEl::kSegment(), {1, 2}, std::make_unique<lf::geometry::SegmentO1>(node_coord));


  // build the mesh
  auto mesh = test.Build(); // mesh is Mesh object

  int dim_mesh_MeshFactory = test.DimWorld(); // test is MeshFactory object
  int dim_mesh_Mesh = mesh->DimWorld();



  //std::cout << "\n" << mesh << std::endl; // Use PrintInfo

  // Mesh
  //lf::mesh::utils::PrintInfo(*mesh, std::cout); // Prints


  // Entities

  lf::mesh::utils::PrintInfo(lf::mesh::hybrid2dp::Point(), std::cout); // Dimension = 2?
  lf::mesh::utils::PrintInfo(lf::mesh::hybrid2dp::Segment(), std::cout); // Dimension = 1
  lf::mesh::utils::PrintInfo(lf::mesh::hybrid2dp::Triangle(), std::cout); // Prints
  lf::mesh::utils::PrintInfo(lf::mesh::hybrid2dp::Quadrilateral(), std::cout); // Prints

  //std::cout << lf::mesh::hybrid2dp::Point() << std::endl;
  //std::cout << lf::mesh::hybrid2dp::Point() << std::endl;



  return 0L;

}
