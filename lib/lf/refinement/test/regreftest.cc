/** @file regrefttest.cc
 * test for local refinement of linear geometries
 */

#include <gtest/gtest.h>
#include <lf/refinement/refinement_hierarchy.h>
#include <lf/refinement/refutils.h>
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
#include "lf/mesh/utils/utils.h"
#include <iostream>

namespace lf::refinement::test {
  TEST(RegRefTest,RegRef) {
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
    // lf::refinement::MeshHierarchy::output_ctrl_ = 1000;
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

    std::cout << "Writing parent information" << std::endl;
    WriteMatlabLevel(multi_mesh,1,"allref_pi.m");
    
    // Printing mesh information
    lf::mesh::utils::PrintInfo(fine_mesh,std::cout);
  }

  TEST(RegRefTest,LocalRefinement) {
    std::cout << "TEST: Marked edges in the unit square" << std::endl;
    
    // Generate test mesh
    auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
    // Output mesh information
    lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
    // Build mesh hierarchy
    std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
    lf::refinement::MeshHierarchy multi_mesh(mesh_p,*mesh_factory_ptr);

    // Mark all edges
    std::function<bool(const lf::mesh::Mesh &,const lf::mesh::Entity &edge)>  marker
      = [](const lf::mesh::Mesh &mesh,const lf::mesh::Entity &edge) -> bool {
      Eigen::MatrixXd ref_c(1,1); ref_c(0,0) = 0.5;
      Eigen::VectorXd c(edge.Geometry()->Global(ref_c));
      return ((c[0] > 0.0) && (c[0] < 1.0) && (c[1] > 0.0) && (c[1] < 1.0));
    };

    std::cout << "#### Marking edges" << std::endl;
    multi_mesh.MarkEdges(marker);

    // Refine uniformly
    std::cout << "#### Refining mesh locally" << std::endl;
    multi_mesh.RefineMarked();
    
    auto &fine_mesh = multi_mesh.getMesh(1);
    std::cout << "Checking mesh completeness" << std::endl;
    lf::mesh::test_utils::checkMeshCompleteness(fine_mesh);
    
    std::cout << "Writing MATLAB file" << std::endl;
    lf::mesh::utils::writeMatlab(fine_mesh, "locref.m");

    std::cout << "Writing parent information" << std::endl;
    WriteMatlabLevel(multi_mesh,1,"locref_pi.m");
    
    // Printing mesh information
    lf::mesh::utils::PrintInfo(fine_mesh,std::cout);
  }
} 
