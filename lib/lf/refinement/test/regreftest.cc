/** @file regrefttest.cc
 * test for local refinement of linear geometries
 */

#include <gtest/gtest.h>
#include <lf/refinement/refinement_hierarchy.h>
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"
#include <iostream>

namespace lf::refinement::test {
  TEST(RegRefTest,RegRef) {
    // Generate test mesh
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
    // Output mesh information
    lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
    // Build mesh hierarchy
    std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
    lf::refinement::MeshHierarchy multi_mesh(mesh_p,*mesh_factory_ptr);
    multi_mesh.RefineRegular();
    lf::mesh::utils::PrintInfo(multi_mesh.getMesh(1),std::cout);
  }
} 
