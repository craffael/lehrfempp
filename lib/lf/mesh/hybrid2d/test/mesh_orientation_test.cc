/**
 * @file
 * @brief Tests for RelativeOrientations() method
 * @author Raffael Casagrande
 * @date   2018-07-01 01:15:55
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <iostream>

#include "lf/mesh/test_utils/test_meshes.h"
#include "mesh_factory_test.h"

namespace lf::mesh::hybrid2d::test {
// Output orientations of edges for test mesh
TEST(lf_hybrid2d, lf_orientation) {
  std::cout << "### TEST: mesh orientation" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  EXPECT_EQ(mesh_p->NumEntities(0), 9) << "Test mesh: 9 cells expected!";
  EXPECT_EQ(mesh_p->NumEntities(1), 18) << "Test mesh: 18 edges expected!";
  EXPECT_EQ(mesh_p->NumEntities(2), 10) << "Test mesh: 10 nodes expected!";

  std::vector<std::vector<lf::mesh::Orientation>> o{
      {lf::mesh::Orientation::positive, lf::mesh::Orientation::positive,
       lf::mesh::Orientation::positive},
      {lf::mesh::Orientation::positive, lf::mesh::Orientation::positive,
       lf::mesh::Orientation::positive},
      {lf::mesh::Orientation::positive, lf::mesh::Orientation::positive,
       lf::mesh::Orientation::negative},
      {lf::mesh::Orientation::positive, lf::mesh::Orientation::positive,
       lf::mesh::Orientation::negative},
      {lf::mesh::Orientation::positive, lf::mesh::Orientation::positive,
       lf::mesh::Orientation::negative},
      {lf::mesh::Orientation::negative, lf::mesh::Orientation::positive,
       lf::mesh::Orientation::negative, lf::mesh::Orientation::positive},
      {lf::mesh::Orientation::negative, lf::mesh::Orientation::positive,
       lf::mesh::Orientation::positive, lf::mesh::Orientation::positive},
      {lf::mesh::Orientation::negative, lf::mesh::Orientation::positive,
       lf::mesh::Orientation::negative},
      {lf::mesh::Orientation::negative, lf::mesh::Orientation::negative,
       lf::mesh::Orientation::negative}};

  // Run through cells and output edge orientations
  for (const mesh::Entity* cell : mesh_p->Entities(0)) {
    const lf::base::glb_idx_t cell_idx(mesh_p->Index(*cell));
    std::cout << *cell << ' ' << cell_idx << ": ";
    // Array of edges
    auto edges = cell->SubEntities(1);
    // Array of orientations
    auto oris = cell->RelativeOrientations();
    for (int ed_idx = 0; ed_idx < cell->RefEl().NumSubEntities(1); ed_idx++) {
      std::cout << ", edge " << mesh_p->Index(*edges[ed_idx])
                << ": or = " << to_char(oris[ed_idx]) << ' ';
      EXPECT_EQ(o[cell_idx][ed_idx], oris[ed_idx])
          << "Orientation mismatch, cell " << cell_idx << ", edge " << ed_idx;
    }
    std::cout << std::endl;
  }
}

TEST(lf_hybrid_2d, Orientation) {
  std::cout << "### TEST: opposite orientations" << std::endl;
  // Building a triangular mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  // For distinguishing boundary edges
  auto ed_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  // For accumulating orientations
  lf::mesh::utils::CodimMeshDataSet<int> sum_ori(mesh_p, 1, 0);
  // Visit all cells
  for (const mesh::Entity* cell : mesh_p->Entities(0)) {
    // Array of edges
    auto edges = cell->SubEntities(1);
    // Array of orientations
    auto oris = cell->RelativeOrientations();
    for (int ed_idx = 0; ed_idx < cell->RefEl().NumSubEntities(1); ed_idx++) {
      sum_ori(*edges[ed_idx]) += (int)oris[ed_idx];
    }
  }
  for (const mesh::Entity* edge : mesh_p->Entities(1)) {
    if (ed_bd_flags(*edge)) {
      EXPECT_EQ(std::abs(sum_ori(*edge)), 1)
          << "Wrong orientation of boundary edge";
    } else {
      EXPECT_EQ(sum_ori(*edge), 0)
          << "Orientation problem: edge " << mesh_p->Index(*edge);
    }
  }
}

}  // namespace lf::mesh::hybrid2d::test
