/** @file regrefttest.cc
 * test for local refinement of linear geometries
 */

#include <gtest/gtest.h>
#include <lf/refinement/refinement_hierarchy.h>
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
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

    auto &fine_mesh = multi_mesh.getMesh(1);
    std::cout << "Checking mesh completeness" << std::endl;
    lf::mesh::test_utils::checkMeshCompleteness(fine_mesh);
    
    std::cout << "Writing MATLAB file" << std::endl;
    lf::mesh::utils::writeMatlab(fine_mesh, "fine_mesh.m");

    // Printing mesh information
    lf::mesh::utils::PrintInfo(fine_mesh,std::cout);
    
  }
} 
