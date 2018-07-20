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
    lf::refinement::MeshHierarchy::output_ctrl_ = 1000;
    std::cout << "TEST: Uniform regular refinement" << std::endl;
    // Generate test mesh
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
    // Output mesh information
    lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
    // Build mesh hierarchy
    std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
    lf::refinement::MeshHierarchy multi_mesh(mesh_p,*mesh_factory_ptr);

    std::cout << "RegRefTEST: Regular refinement" << std::endl;
    multi_mesh.RefineRegular();

    auto &fine_mesh = multi_mesh.getMesh(1);
    std::cout << "Checking mesh completeness" << std::endl;
    lf::mesh::test_utils::checkMeshCompleteness(fine_mesh);
    
    std::cout << "Writing MATLAB file" << std::endl;
    lf::mesh::utils::writeMatlab(fine_mesh, "fine_mesh.m");

    // Printing mesh information
    lf::mesh::utils::PrintInfo(fine_mesh,std::cout);
    
  }

  TEST(RegRefTest,BarycentricRef) {
    std::cout << "TEST: Uniform barycentric refinement" << std::endl;
    // Generate test mesh
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
    // Output mesh information
    lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
    // Build mesh hierarchy
    std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
    lf::refinement::MeshHierarchy multi_mesh(mesh_p,*mesh_factory_ptr);

    multi_mesh.RefineRegular(lf::refinement::RefPat::rp_barycentric);

    auto &fine_mesh = multi_mesh.getMesh(1);
    std::cout << "Checking mesh completeness" << std::endl;
    lf::mesh::test_utils::checkMeshCompleteness(fine_mesh);
    
    std::cout << "Writing MATLAB file" << std::endl;
    lf::mesh::utils::writeMatlab(fine_mesh, "barycentric_ref.m");

    // Printing mesh information
    lf::mesh::utils::PrintInfo(fine_mesh,std::cout);
  }

  TEST(RegRefTest,AllMarkedRefinement) {
    std::cout << "TEST: All edges marked" << std::endl;
    
    // Generate test mesh
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
    // Output mesh information
    lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
    // Build mesh hierarchy
    std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
    lf::refinement::MeshHierarchy multi_mesh(mesh_p,*mesh_factory_ptr);

    // Mark all edges
    std::function<bool(const lf::mesh::Mesh &,const lf::mesh::Entity &)>  marker
      = [](const lf::mesh::Mesh &mesh,const lf::mesh::Entity &edge) -> bool { return true; };

    std::cout << "#### Marking edges" << std::endl;
    multi_mesh.MarkEdges(marker);

    // Refine uniformly
    std::cout << "#### Refining mesh" << std::endl;
    multi_mesh.RefineMarked();
    
    auto &fine_mesh = multi_mesh.getMesh(1);
    std::cout << "Checking mesh completeness" << std::endl;
    lf::mesh::test_utils::checkMeshCompleteness(fine_mesh);
    
    std::cout << "Writing MATLAB file" << std::endl;
    lf::mesh::utils::writeMatlab(fine_mesh, "allref.m");

    // Printing mesh information
    lf::mesh::utils::PrintInfo(fine_mesh,std::cout);
  }
} 
