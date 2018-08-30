/**
 * @file
 * @brief Tests for the lf::mesh::hybrid2dp::MeshFactory
 * @author Raffael Casagrande
 * @date   2018-07-01 01:15:55
 * @copyright MIT License
 */

#include <Eigen/Eigen>
#include <iostream>
#include <gtest/gtest.h>

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

  // Run through cells and output edge orientations
  for (const mesh::Entity &cell : mesh_p->Entities(0)) {
    std::cout << "Cell " << mesh_p->Index(cell) << ": ";
    // Array of edges
    lf::base::RandomAccessRange<const lf::mesh::Entity> edges(
        cell.SubEntities(1));
    // Array of orientations
    lf::base::RandomAccessRange<const lf::mesh::Orientation> oris(
        cell.RelativeOrientations());
    for (int ed_idx = 0; ed_idx < cell.RefEl().NumSubEntities(1); ed_idx++) {
      std::cout << ", edge " << mesh_p->Index(edges[ed_idx])
		<< ": orientation = " << to_char(oris[ed_idx]) << ' ';
    }
    std::cout << std::endl;
  }
}

}  // namespace lf::mesh::hybrid2dp::test
