/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Functions for the initialization of special data sets
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */
#include <gtest/gtest.h>
#include <lf/mesh/utils/utils.h>
#include <iostream>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::mesh::utils::test {
// Test for counting of super-entities
TEST(test_mesh_utils, count_test) {
  std::cout << "Test for counting of super-entities" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  EXPECT_EQ(mesh_p->NumEntities(0), 9) << "Test mesh: 9 cells expected!";
  EXPECT_EQ(mesh_p->NumEntities(1), 18) << "Test mesh: 18 edges expected!";
  EXPECT_EQ(mesh_p->NumEntities(2), 10) << "Test mesh: 10 nodes expected!";

  auto cells_at_edges{countNoSuperEntities(mesh_p, 1, 1)};
  for (const lf::mesh::Entity &edge : mesh_p->Entities(1)) {
    std::cout << cells_at_edges(edge) << " cells @ " << edge << ' '
              << mesh_p->Index(edge) << std::endl;
  }
  auto edges_at_vertices{countNoSuperEntities(mesh_p, 2, 1)};
  for (const lf::mesh::Entity &node : mesh_p->Entities(2)) {
    std::cout << edges_at_vertices(node) << " edges @ " << node << ' '
              << mesh_p->Index(node) << std::endl;
  }
  auto cells_at_nodes{countNoSuperEntities(mesh_p, 2, 2)};
  for (const lf::mesh::Entity &node : mesh_p->Entities(2)) {
    std::cout << cells_at_nodes(node) << " cells @ " << node << ' '
              << mesh_p->Index(node) << std::endl;
  }
}

TEST(test_mesh_utils, bd_flag_test) {
  std::cout << "Test for flagging of boundary entities" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  auto bd_flags{flagEntitiesOnBoundary(mesh_p)};
  for (const lf::mesh::Entity &edge : mesh_p->Entities(1)) {
    std::cout << edge << ' ' << mesh_p->Index(edge)
              << (bd_flags(edge) ? (" ON BOUNDARY") : (" INTERIOR"))
              << std::endl;
  }
  for (const lf::mesh::Entity &node : mesh_p->Entities(2)) {
    std::cout << node << ' ' << mesh_p->Index(node)
              << (bd_flags(node) ? (" ON BOUNDARY") : (" INTERIOR"))
              << std::endl;
  }
}

}  // namespace lf::mesh::utils::test
