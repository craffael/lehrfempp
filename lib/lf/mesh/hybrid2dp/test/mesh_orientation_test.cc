/**
 * @file
 * @brief Tests for the lf::mesh::hybrid2dp::MeshFactory
 * @author Raffael Casagrande
 * @date   2018-07-01 01:15:55
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <Eigen/Eigen>
#include <iostream>

#include <lf/mesh/hybrid2dp/hybrid2dp.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"
#include "mesh_factory_test.h"

namespace lf::mesh::hybrid2dp::test {
// Output orientations of edges for test mesh
TEST(lf_hybrid2dp, lf_orientation) {
  std::cout << "### TEST: mesh orientation" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  EXPECT_EQ(mesh_p->Size(0), 9) << "Test mesh: 9 cells expected!";
  EXPECT_EQ(mesh_p->Size(1), 18) << "Test mesh: 18 edges expected!";
  EXPECT_EQ(mesh_p->Size(2), 10) << "Test mesh: 10 nodes expected!";

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
  for (const mesh::Entity &cell : mesh_p->Entities(0)) {
    const lf::base::glb_idx_t cell_idx(mesh_p->Index(cell));
    std::cout << cell << ' ' << cell_idx << ": ";
    // Array of edges
    lf::base::RandomAccessRange<const lf::mesh::Entity> edges(
        cell.SubEntities(1));
    // Array of orientations
    lf::base::RandomAccessRange<const lf::mesh::Orientation> oris(
        cell.RelativeOrientations());
    for (int ed_idx = 0; ed_idx < cell.RefEl().NumSubEntities(1); ed_idx++) {
      std::cout << ", edge " << mesh_p->Index(edges[ed_idx])
                << ": or = " << to_char(oris[ed_idx]) << ' ';
      EXPECT_EQ(o[cell_idx][ed_idx], oris[ed_idx])
          << "Orientation mismatch, cell " << cell_idx << ", edge " << ed_idx;
    }
    std::cout << std::endl;
  }
}

}  // namespace lf::mesh::hybrid2dp::test
