/** @file mesh_entity_output.cc
 * Demo for outputting information about mesh entity elements
 */

#include <iostream>

#include "lf/mesh/hybrid2d/entity.h"
#include "lf/mesh/hybrid2dp/triangle.h"
#include "lf/mesh/mesh_factory.h"
#include "lf/mesh/entity.h"
#include "lf/base/base.h"
#include "lf/mesh/tp_triag_mesh_builder.h"
#include "lf/mesh/hybrid2d/hybrid2d.h"

/*
#include <gtest/gtest.h>
#include <lf/mesh/hybrid2dp/hybrid2dp.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
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

  lf::mesh::hybrid2d::MeshFactory test(2); // Working, but not entity
  // Need to add nodes and entity, AddPoint, AddEntity

  // add nodes

  //test.AddPoint(Eigen::Vector2cd(0,0));

/*
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, 0)), 0);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(1, 0)), 1);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(1, 1)), 2);
  EXPECT_EQ(mf.AddPoint(Eigen::Vector2d(0, 1)), 3);

  // add an element
  Eigen::MatrixXd node_coord(2, 4);
  node_coord << 0, 1, 1, 0, 0, 0, 1, 1;
  EXPECT_EQ(0, mf.AddEntity(
                   base::RefEl::kQuad(), {0, 1, 2, 3},
                   std::make_unique<geometry::QuadO1>(std::move(node_coord))));
*/


  return 0L;
}
