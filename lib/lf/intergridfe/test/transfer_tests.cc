#include <gtest/gtest.h>

#include <lf/intergridfe/intergridfe.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/refinement/mesh_hierarchy.h>

#include <lf/io/vtk_writer.h>

using lf::uscalfe::operator-;

TEST(lf_intergridfe, NormLinear) {
  // Generate a test mesh, refine it and compare the norm of the coarse and the
  // fine version
  for (int selector = 0; selector < 9; ++selector) {
    // Construct a mesh hierarchy from a test mesh
    const auto mesh_coarse =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
    auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    auto hierarchy =
        lf::refinement::MeshHierarchy(mesh_coarse, std::move(factory));
    hierarchy.RefineRegular();
    const auto mesh_fine = hierarchy.getMesh(1);
    // Construct the fe spaces on the coarse and on the fine mesh
    const auto fes_coarse =
        std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_coarse);
    const auto fes_fine =
        std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_fine);
    // Iterate over all basis functions on the coarse mesh and transfer them to
    // the fine mesh
    const auto& dofh_coarse = fes_coarse->LocGlobMap();
    const auto& dofh_fine = fes_fine->LocGlobMap();
    const lf::base::size_type dofcount_coarse = dofh_coarse.NumDofs();
    for (lf::base::size_type dofidx = 0; dofidx < dofcount_coarse; ++dofidx) {
      // Construct the mesh function on the coarse mesh
      Eigen::VectorXd dofvector_coarse = Eigen::VectorXd::Zero(dofcount_coarse);
      dofvector_coarse[dofidx] = 1;
      // Transfer the mesh function one level down in the hierarchy
      const auto dofvector_fine = lf::intergridfe::prolongate(hierarchy, fes_coarse, fes_fine, dofvector_coarse, 0);
      lf::uscalfe::MeshFunctionFE mf_fine(fes_fine, dofvector_fine);
      // Compare the norms of the two mesh functions
      lf::uscalfe::MeshFunctionFE mf_coarse(fes_coarse, dofvector_coarse);
      const auto qr_provider = [](const lf::mesh::Entity& e) {
        return lf::quad::make_QuadRule(e.RefEl(), 2);
      };
      const double norm_coarse = lf::uscalfe::IntegrateMeshFunction(
          *mesh_coarse, lf::uscalfe::squaredNorm(mf_coarse), qr_provider);
      const double norm_fine = lf::uscalfe::IntegrateMeshFunction(
          *mesh_fine, lf::uscalfe::squaredNorm(mf_fine), qr_provider);
      ASSERT_DOUBLE_EQ(norm_coarse, norm_fine)
          << "dofidx = " << dofidx << "\ndofvector_coarse = ["
          << dofvector_coarse.transpose() << "]\ndofvector_fine = ["
          << dofvector_fine.transpose() << "]\n";
    }
  }
}

TEST(lf_intergridfe, NormQuadratic) {
  // Generate a test mesh, refine it and compare the norm of the coarse and the
  // fine version
  for (int selector = 0; selector < 9; ++selector) {
    // Construct a mesh hierarchy from a test mesh
    const auto mesh_coarse =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
    auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    auto hierarchy =
        lf::refinement::MeshHierarchy(mesh_coarse, std::move(factory));
    hierarchy.RefineRegular();
    const auto mesh_fine = hierarchy.getMesh(1);
    // Construct the fe spaces on the coarse and on the fine mesh
    const auto fes_coarse =
        std::make_shared<const lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_coarse);
    const auto fes_fine =
        std::make_shared<const lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_fine);
    // Iterate over all basis functions on the coarse mesh and transfer them to
    // the fine mesh
    const auto& dofh_coarse = fes_coarse->LocGlobMap();
    const auto& dofh_fine = fes_fine->LocGlobMap();
    const lf::base::size_type dofcount_coarse = dofh_coarse.NumDofs();
    for (lf::base::size_type dofidx = 0; dofidx < dofcount_coarse; ++dofidx) {
      // Construct the mesh function on the coarse mesh
      Eigen::VectorXd dofvector_coarse = Eigen::VectorXd::Zero(dofcount_coarse);
      dofvector_coarse[dofidx] = 1;
      // Transfer the mesh function one level down in the hierarchy
      const auto dofvector_fine = lf::intergridfe::prolongate(hierarchy, fes_coarse, fes_fine, dofvector_coarse, 0);
      lf::uscalfe::MeshFunctionFE mf_fine(fes_fine, dofvector_fine);
      // Compare the norms of the two mesh functions
      lf::uscalfe::MeshFunctionFE mf_coarse(fes_coarse, dofvector_coarse);
      const auto qr_provider = [](const lf::mesh::Entity& e) {
        return lf::quad::make_QuadRule(e.RefEl(), 4);
      };
      const double norm_coarse = lf::uscalfe::IntegrateMeshFunction(
          *mesh_coarse, lf::uscalfe::squaredNorm(mf_coarse), qr_provider);
      const double norm_fine = lf::uscalfe::IntegrateMeshFunction(
          *mesh_fine, lf::uscalfe::squaredNorm(mf_fine), qr_provider);
      ASSERT_TRUE(std::fabs(norm_coarse - norm_fine) < 1e-10)
          << "dofidx = " << dofidx << "\ndofvector_coarse = ["
          << dofvector_coarse.transpose() << "]\ndofvector_fine = ["
          << dofvector_fine.transpose() << "]\n";
    }
  }
}
