/** @file regrefttest.cc
 * test for local refinement of linear geometries
 */

#include <iostream>
#include "lf/io/io.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"
#include "refinement_test_utils.h"

namespace lf::refinement::test {

TEST(RegRefTest, RegRef) {
  std::cout << "TEST: Uniform regular refinement" << std::endl;
  // Generate test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Output mesh information
  lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
  // Build mesh hierarchy
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(mesh_p, mesh_factory_ptr);

  std::cout << "RegRefTEST: Regular refinement" << std::endl;
  multi_mesh.RefineRegular();

  std::shared_ptr<const mesh::Mesh> fine_mesh = multi_mesh.getMesh(1);

  std::cout << "Checking mesh completeness" << std::endl;
  lf::mesh::test_utils::checkMeshCompleteness(*fine_mesh);

  std::cout << "Checking geometry compatibulity: " << std::flush;
  lf::mesh::test_utils::watertight_mesh_ctrl = 100;
  auto fails = lf::mesh::test_utils::isWatertightMesh(*fine_mesh, false);
  EXPECT_EQ(fails.size(), 0) << "Inconsistent geometry!";
  if (fails.size() == 0) {
    std::cout << "consistent!" << std::endl;
  } else {
    std::cout << "INCONSISTENT!" << std::endl;
    for (auto &geo_errs : fails) {
      std::cout << geo_errs.first.ToString() << "(" << geo_errs.second << ")"
                << std::endl;
    }
  }
  // Verify that child mesh is properly connected to father mesh
  checkFatherChildRelations(multi_mesh, 0);

  std::cout << "Writing MATLAB file" << std::endl;
  lf::mesh::utils::writeMatlab(*fine_mesh, "fine_mesh.m");

  std::cout << "Writing parent information" << std::endl;
  WriteMatlabLevel(multi_mesh, 1, "fine_mesh_pi.m");

  // Printing mesh information
  lf::mesh::utils::PrintInfo(*fine_mesh, std::cout);

  // Output mesh geometry in TikZ format
  // Enable, once function is available in master branch
  // lf::mesh::utils::writeTikZ(*fine_mesh, "fine_mesh.tikz",
  // 			     lf::mesh::utils::TikzOutputCtrl::RenderCells |
  // 			     lf::mesh::utils::TikzOutputCtrl::CellNumbering |
  // 			     lf::mesh::utils::TikzOutputCtrl::VerticeNumbering |
  // 			     lf::mesh::utils::TikzOutputCtrl::NodeNumbering |
  // 			     lf::mesh::utils::TikzOutputCtrl::EdgeNumbering
  // 			     );
}  // End RegRefTest::RegRef

TEST(RegRefTest, BarycentricRef) {
  std::cout << "TEST: Uniform barycentric refinement" << std::endl;
  // Generate test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Output mesh information
  lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
  // Build mesh hierarchy
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(mesh_p, mesh_factory_ptr);

  multi_mesh.RefineRegular(lf::refinement::RefPat::rp_barycentric);

  std::shared_ptr<const mesh::Mesh> fine_mesh = multi_mesh.getMesh(1);

  std::cout << "Checking mesh completeness" << std::endl;
  lf::mesh::test_utils::checkMeshCompleteness(*fine_mesh);

  std::cout << "Checking geometry compatibulity: " << std::flush;
  lf::mesh::test_utils::watertight_mesh_ctrl = 100;
  auto fails = lf::mesh::test_utils::isWatertightMesh(*fine_mesh, false);
  EXPECT_EQ(fails.size(), 0) << "Inconsistent geometry!";
  if (fails.size() == 0) {
    std::cout << "consistent!" << std::endl;
  } else {
    std::cout << "INCONSISTENT!" << std::endl;
    for (auto &geo_errs : fails) {
      std::cout << geo_errs.first.ToString() << "(" << geo_errs.second << ")"
                << std::endl;
    }
  }
  // Verify that child mesh is properly connected to father mesh
  checkFatherChildRelations(multi_mesh, 0);

  std::cout << "Writing MATLAB file" << std::endl;
  lf::mesh::utils::writeMatlab(*fine_mesh, "barycentric_ref.m");

  std::cout << "Writing parent information" << std::endl;
  WriteMatlabLevel(multi_mesh, 1, "barycentric_ref_pi.m");

  // Printing mesh information
  lf::mesh::utils::PrintInfo(*fine_mesh, std::cout);
}

TEST(RegRefTest, AllMarkedRefinement) {
  lf::refinement::MeshHierarchy::output_ctrl_ = 0;
  std::cout << "TEST: All edges marked" << std::endl;

  // Generate test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Output mesh information
  lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
  // Build mesh hierarchy
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(mesh_p, mesh_factory_ptr);

  // Mark all edges
  std::function<bool(const lf::mesh::Mesh &, const lf::mesh::Entity &)> marker =
      [](const lf::mesh::Mesh &mesh, const lf::mesh::Entity &edge) -> bool {
    return true;
  };

  std::cout << "#### Marking edges" << std::endl;
  multi_mesh.MarkEdges(marker);

  // Refine uniformly
  std::cout << "#### Refining mesh" << std::endl;
  multi_mesh.RefineMarked();

  std::shared_ptr<const mesh::Mesh> fine_mesh = multi_mesh.getMesh(1);

  std::cout << "Checking mesh completeness" << std::endl;
  lf::mesh::test_utils::checkMeshCompleteness(*fine_mesh);

  std::cout << "Checking geometry compatibulity: " << std::flush;
  lf::mesh::test_utils::watertight_mesh_ctrl = 100;
  auto fails = lf::mesh::test_utils::isWatertightMesh(*fine_mesh, false);
  EXPECT_EQ(fails.size(), 0) << "Inconsistent geometry!";
  if (fails.size() == 0) {
    std::cout << "consistent!" << std::endl;
  } else {
    std::cout << "INCONSISTENT!" << std::endl;
    for (auto &geo_errs : fails) {
      std::cout << geo_errs.first.ToString() << "(" << geo_errs.second << ")"
                << std::endl;
    }
  }
  // Verify that child mesh is properly connected to father mesh
  checkFatherChildRelations(multi_mesh, 0);

  std::cout << "Writing MATLAB file" << std::endl;
  lf::mesh::utils::writeMatlab(*fine_mesh, "allref.m");

  std::cout << "Writing parent information" << std::endl;
  WriteMatlabLevel(multi_mesh, 1, "allref_pi.m");

  // Printing mesh information
  lf::mesh::utils::PrintInfo(*fine_mesh, std::cout);
}

TEST(LocRefTest, LocalRefinement) {
  lf::refinement::MeshHierarchy::output_ctrl_ = 0;
  std::cout << "TEST: Marked edges in the unit square" << std::endl;

  // Generate test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Output mesh information
  lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
  // Build mesh hierarchy
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(mesh_p, mesh_factory_ptr);

  // Mark edges whose center lies inside a square
  std::function<bool(const lf::mesh::Mesh &, const lf::mesh::Entity &edge)>
      marker =
          [](const lf::mesh::Mesh &mesh, const lf::mesh::Entity &edge) -> bool {
    Eigen::MatrixXd ref_c(1, 1);
    ref_c(0, 0) = 0.5;
    Eigen::VectorXd c(edge.Geometry()->Global(ref_c));
    return ((c[0] > 1.0) && (c[0] < 2.0) && (c[1] > 1.0) && (c[1] < 2.0));
  };

  std::cout << "#### Marking edges" << std::endl;
  multi_mesh.MarkEdges(marker);

  // Refine uniformly
  std::cout << "#### Refining mesh locally" << std::endl;
  multi_mesh.RefineMarked();

  std::shared_ptr<const mesh::Mesh> fine_mesh = multi_mesh.getMesh(1);

  std::cout << "Checking mesh completeness" << std::endl;
  lf::mesh::test_utils::checkMeshCompleteness(*fine_mesh);

  std::cout << "Checking geometry compatibulity: " << std::flush;
  lf::mesh::test_utils::watertight_mesh_ctrl = 100;
  auto fails = lf::mesh::test_utils::isWatertightMesh(*fine_mesh, false);
  EXPECT_EQ(fails.size(), 0) << "Inconsistent geometry!";
  if (fails.size() == 0) {
    std::cout << "consistent!" << std::endl;
  } else {
    std::cout << "INCONSISTENT!" << std::endl;
    for (auto &geo_errs : fails) {
      std::cout << geo_errs.first.ToString() << "(" << geo_errs.second << ")"
                << std::endl;
    }
  }
  // Verify that child mesh is properly connected to father mesh
  checkFatherChildRelations(multi_mesh, 0);

  std::cout << "Writing MATLAB file" << std::endl;
  lf::mesh::utils::writeMatlab(*fine_mesh, "locref.m");

  std::cout << "Writing parent information" << std::endl;
  WriteMatlabLevel(multi_mesh, 1, "locref_pi.m");

  // Printing mesh information
  lf::mesh::utils::PrintInfo(*fine_mesh, std::cout);
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

TEST(LocRefTest, MultipleRefinement) {
  lf::refinement::MeshHierarchy::output_ctrl_ = 0;
  lf::mesh::test_utils::watertight_mesh_ctrl = 100;

  const size_type Nrefs = 4;
  std::cout << "TEST: Multiple refinement with Marked edges in a square"
            << std::endl;

  // Generate the standard hybrid test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Output mesh information
  lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
  // Build mesh hierarchy
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(mesh_p, mesh_factory_ptr);

  // Mark edges whose midpoints are located in a certain region
  std::function<bool(const lf::mesh::Mesh &, const lf::mesh::Entity &edge)>
      marker =
          [](const lf::mesh::Mesh &mesh, const lf::mesh::Entity &edge) -> bool {
    Eigen::MatrixXd ref_c(1, 1);
    ref_c(0, 0) = 0.5;
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

    const size_type n_levels = multi_mesh.NumLevels();

    std::shared_ptr<const mesh::Mesh> mesh = multi_mesh.getMesh(n_levels - 1);

    std::cout << "#### Mesh on level " << n_levels - 1 << ": " << mesh->Size(2)
              << " nodes, " << mesh->Size(1) << " nodes, " << mesh->Size(0)
              << " cells," << std::endl;

    std::cout << "Checking mesh completeness" << std::endl;
    lf::mesh::test_utils::checkMeshCompleteness(*mesh);

    std::cout << "RefStep " << refstep + 1
              << ": Checking geometry compatibulity: " << std::flush;
    auto fails = lf::mesh::test_utils::isWatertightMesh(*mesh, false);
    EXPECT_EQ(fails.size(), 0) << "Inconsistent geometry!";
    if (fails.size() == 0) {
      std::cout << "consistent!" << std::endl;
    } else {
      std::cout << "INCONSISTENT!" << std::endl;
      for (auto &geo_errs : fails) {
        std::cout << geo_errs.first.ToString() << "(" << geo_errs.second << ")"
                  << std::endl;
      }
    }

    // Printing mesh information
    lf::mesh::utils::PrintInfo(*mesh, std::cout);
  }
  WriteMatlab(multi_mesh, "multiref");
}

TEST(LocRefTest, MixedRefinement) {
  lf::mesh::test_utils::watertight_mesh_ctrl = 100;
  lf::refinement::MeshHierarchy::output_ctrl_ = 0;
  std::cout << "TEST: Mixed local/global refinement" << std::endl;

  // Generate the standard hybrid test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Output mesh information
  lf::mesh::utils::PrintInfo(*mesh_p, std::cout);
  // Build mesh hierarchy
  std::shared_ptr<lf::mesh::hybrid2dp::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2dp::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(mesh_p, mesh_factory_ptr);

  // Mark edges whose midpoints are located in a certain region
  std::function<bool(const lf::mesh::Mesh &, const lf::mesh::Entity &edge)>
      marker =
          [](const lf::mesh::Mesh &mesh, const lf::mesh::Entity &edge) -> bool {
    Eigen::MatrixXd ref_c(1, 1);
    ref_c(0, 0) = 0.5;
    Eigen::VectorXd c(edge.Geometry()->Global(ref_c));
    return ((c[0] > 1.0) && (c[0] < 2.0) && (c[1] > 1.0) && (c[1] < 2.0));
  };

  // First step regular refinement
  multi_mesh.RefineRegular();
  // Second step local refinement
  multi_mesh.MarkEdges(marker);
  multi_mesh.RefineMarked();
  // Third step: barycentric refinement
  multi_mesh.RefineRegular(lf::refinement::RefPat::rp_barycentric);
  // Fourth step: local refinement again
  multi_mesh.MarkEdges(marker);
  multi_mesh.RefineMarked();

  // Check mesh integrity
  const size_type n_levels = multi_mesh.NumLevels();
  EXPECT_EQ(n_levels, 5) << "After four steps of refinement " << n_levels
                         << " level?";
  for (int l = 0; l < n_levels; ++l) {
    std::shared_ptr<mesh::Mesh> mesh = multi_mesh.getMesh(l);

    std::cout << "### LEVEL " << l << " ####" << std::endl;
    std::cout << "#### Mesh on level " << n_levels - 1 << ": " << mesh->Size(2)
              << " nodes, " << mesh->Size(1) << " nodes, " << mesh->Size(0)
              << " cells," << std::endl;

    std::cout << "Checking mesh completeness" << std::endl;
    lf::mesh::test_utils::checkMeshCompleteness(*mesh);

    std::cout << ": Checking geometry compatibulity: " << std::flush;
    auto fails = lf::mesh::test_utils::isWatertightMesh(*mesh, false);
    EXPECT_EQ(fails.size(), 0) << "Inconsistent geometry!";
    if (fails.size() == 0) {
      std::cout << "consistent!" << std::endl;
    } else {
      std::cout << "INCONSISTENT!" << std::endl;
      for (auto &geo_errs : fails) {
        std::cout << geo_errs.first.ToString() << "(" << geo_errs.second << ")"
                  << std::endl;
      }
    }
    // Write mesh to VTK file
    std::stringstream level_asc;
    level_asc << l;
    std::string filename = std::string("mixedref") + level_asc.str() + ".vtk";
    lf::io::VtkWriter(mesh, filename);
  }  // end loop over levels
  WriteMatlab(multi_mesh, "mixedref");
}  // end mixed refinement test

}  // namespace lf::refinement::test
