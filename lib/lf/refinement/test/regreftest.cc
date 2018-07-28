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
    
    std::cout << "Checking geometry compatibulity: " << std::flush;
    lf::mesh::test_utils::watertight_mesh_ctrl = 100;
    auto fails = lf::mesh::test_utils::isWatertightMesh(fine_mesh,false);
    EXPECT_EQ(fails.size(),0) << "Inconsistent geometry!";
    if (fails.size() == 0) {
      std::cout << "consistent!" << std::endl;
    }
    else {
      std::cout << "INCONSISTENT!" << std::endl;
      for (auto & geo_errs : fails) {
	std::cout  << geo_errs.first.ToString() << "(" << geo_errs.second << ")" << std::endl;
      }
    }
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

    std::cout << "Checking geometry compatibulity: " << std::flush;
    lf::mesh::test_utils::watertight_mesh_ctrl = 100;
    auto fails = lf::mesh::test_utils::isWatertightMesh(fine_mesh,false);
    EXPECT_EQ(fails.size(),0) << "Inconsistent geometry!";
    if (fails.size() == 0) {
      std::cout << "consistent!" << std::endl;
    }
    else {
      std::cout << "INCONSISTENT!" << std::endl;
      for (auto & geo_errs : fails) {
	std::cout  << geo_errs.first.ToString() << "(" << geo_errs.second << ")" << std::endl;
      }
    }
    
    std::cout << "Writing MATLAB file" << std::endl;
    lf::mesh::utils::writeMatlab(fine_mesh, "barycentric_ref.m");

    // Printing mesh information
    lf::mesh::utils::PrintInfo(fine_mesh,std::cout);
  } 


  TEST(RegRefTest,AllMarkedRefinement) {
    lf::refinement::MeshHierarchy::output_ctrl_ = 0;
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
    
    std::cout << "Checking geometry compatibulity: " << std::flush;
    lf::mesh::test_utils::watertight_mesh_ctrl = 100;
    auto fails = lf::mesh::test_utils::isWatertightMesh(fine_mesh,false);
    EXPECT_EQ(fails.size(),0) << "Inconsistent geometry!";
    if (fails.size() == 0) {
      std::cout << "consistent!" << std::endl;
    }
    else {
      std::cout << "INCONSISTENT!" << std::endl;
      for (auto & geo_errs : fails) {
	std::cout  << geo_errs.first.ToString() << "(" << geo_errs.second << ")" << std::endl;
      }
    }

    std::cout << "Writing MATLAB file" << std::endl;
    lf::mesh::utils::writeMatlab(fine_mesh, "allref.m");

    std::cout << "Writing parent information" << std::endl;
    WriteMatlabLevel(multi_mesh,1,"allref_pi.m");
    
    // Printing mesh information
    lf::mesh::utils::PrintInfo(fine_mesh,std::cout);
  } 

  TEST(LocRefTest,LocalRefinement) {
    lf::refinement::MeshHierarchy::output_ctrl_ = 0;
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
      return ((c[0] > 1.0) && (c[0] < 2.0) && (c[1] > 1.0) && (c[1] < 2.0));
    };

    std::cout << "#### Marking edges" << std::endl;
    multi_mesh.MarkEdges(marker);

    // Refine uniformly
    std::cout << "#### Refining mesh locally" << std::endl;
    multi_mesh.RefineMarked();
    
    auto &fine_mesh = multi_mesh.getMesh(1);
    std::cout << "Checking mesh completeness" << std::endl;
    lf::mesh::test_utils::checkMeshCompleteness(fine_mesh);

    std::cout << "Checking geometry compatibulity: " << std::flush;
    lf::mesh::test_utils::watertight_mesh_ctrl = 100;
    auto fails = lf::mesh::test_utils::isWatertightMesh(fine_mesh,false);
    EXPECT_EQ(fails.size(),0) << "Inconsistent geometry!";
    if (fails.size() == 0) {
      std::cout << "consistent!" << std::endl;
    }
    else {
      std::cout << "INCONSISTENT!" << std::endl;
      for (auto & geo_errs : fails) {
	std::cout  << geo_errs.first.ToString() << "(" << geo_errs.second << ")" << std::endl;
      }
    }

    std::cout << "Writing MATLAB file" << std::endl;
    lf::mesh::utils::writeMatlab(fine_mesh, "locref.m");

    std::cout << "Writing parent information" << std::endl;
    WriteMatlabLevel(multi_mesh,1,"locref_pi.m");
    
    // Printing mesh information
    lf::mesh::utils::PrintInfo(fine_mesh,std::cout);
  } 


  /* MATLAB script for visualizing the output of the next test
 
     % MATLAB script for ploting testing refinement hierarchy

     display('Level 0');
     opts = []; opts.numbers = true; plot_lf_mesh(@multiref_L0,opts);
     title(sprintf('Multirev level %i',0));

     display('Level 1');
     opts = []; opts.numbers = true; plot_lf_mesh(@multiref_L1,opts);
     title(sprintf('Multirev level %i',1));
     opts = []; opts.parents = @multiref_L1_pi; plot_lf_mesh(@multiref_L1,opts);
     title(sprintf('Multirev parent into %i',1));

     display('Level 2');
     opts = []; opts.numbers = true; plot_lf_mesh(@multiref_L2,opts);
     title(sprintf('Multirev level %i',2));
     opts = []; opts.parents = @multiref_L2_pi; plot_lf_mesh(@multiref_L2,opts);
     title(sprintf('Multirev parent into %i',2));

     display('Level 3');
     opts = []; opts.numbers = true; plot_lf_mesh(@multiref_L3,opts);
     title(sprintf('Multirev level %i',3));
     opts = []; opts.parents = @multiref_L3_pi; plot_lf_mesh(@multiref_L3,opts);
     title(sprintf('Multirev parent into %i',3));

     display('Level 4');
     opts = []; opts.numbers = true; plot_lf_mesh(@multiref_L4,opts);
     title(sprintf('Multirev level %i',4));
     opts = []; opts.parents = @multiref_L4_pi; plot_lf_mesh(@multiref_L4,opts);
     title(sprintf('Multirev parent into %i',4));

  */
  
  TEST(LocRefTest,MultipleRefinement) {
    lf::refinement::MeshHierarchy::output_ctrl_ = 0;
    const size_type Nrefs = 4;
    std::cout << "TEST: Multiple refinement with Marked edges in a square" << std::endl;
    
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
      return ((c[0] > 1.0) && (c[0] < 2.0) && (c[1] > 1.0) && (c[1] < 2.0));
    };

    // Main refinement loop
    std::cout << "#### Entering main refinement loop" << std::endl;
    for (int refstep = 0; refstep < Nrefs; refstep++) {
      std::cout << "#### Marking edges" << std::endl;
      multi_mesh.MarkEdges(marker);
      // Refine locally
      std::cout << "#### Refining mesh locally" << std::endl;
      multi_mesh.RefineMarked();

      const size_type n_levels = multi_mesh.numLevels();
      const lf::mesh::Mesh &mesh(multi_mesh.getMesh(n_levels-1));
      std::cout << "#### Mesh on level " << n_levels-1 << ": "
		<< mesh.Size(2) << " nodes, "
		<< mesh.Size(1) << " nodes, "
		<< mesh.Size(0) << " cells," << std::endl;

      std::cout << "Checking mesh completeness" << std::endl;
      lf::mesh::test_utils::checkMeshCompleteness(mesh);
      
      // Printing mesh information
      lf::mesh::utils::PrintInfo(mesh,std::cout);
    }
    WriteMatlab(multi_mesh,"multiref");
  } 

}
