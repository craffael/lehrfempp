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
        std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(
            mesh_coarse);
    const auto fes_fine =
        std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(
            mesh_fine);
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
      const auto dofvector_fine = lf::intergridfe::prolongate(
          hierarchy, fes_coarse, fes_fine, dofvector_coarse, 0);
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
        std::make_shared<const lf::uscalfe::FeSpaceLagrangeO2<double>>(
            mesh_coarse);
    const auto fes_fine =
        std::make_shared<const lf::uscalfe::FeSpaceLagrangeO2<double>>(
            mesh_fine);
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
      const auto dofvector_fine = lf::intergridfe::prolongate(
          hierarchy, fes_coarse, fes_fine, dofvector_coarse, 0);
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

template <typename FES_COARSE, typename FES_FINE>
void check_lagr_interp_nodes(std::shared_ptr<lf::mesh::Mesh> mesh_coarse,
                             std::shared_ptr<const FES_COARSE> fes_coarse,
                             const Eigen::VectorXd& dofs_coarse,
                             lf::refinement::RefPat ref_pat) {
  // Refine the mesh using the given refinement pattern
  auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  auto mh = lf::refinement::MeshHierarchy(mesh_coarse, std::move(factory));
  mh.RefineRegular(ref_pat);
  const auto mesh_fine = mh.getMesh(1);
  // Initialize the fine fe space
  const auto fes_fine = std::make_shared<const FES_FINE>(mesh_fine);
  // Prolongate the given dof vector
  const auto dofs_fine =
      lf::intergridfe::prolongate(mh, fes_coarse, fes_fine, dofs_coarse, 0);
  // Compare the resulting mesh functions at the interpolation nodes
  lf::uscalfe::MeshFunctionFE mf_coarse(fes_coarse, dofs_coarse);
  lf::uscalfe::MeshFunctionFE mf_fine(fes_fine, dofs_fine);
  for (const auto cell : mesh_fine->Entities(0)) {
    const auto idx = mesh_fine->Index(*cell);
    const auto layout = fes_fine->ShapeFunctionLayout(cell->RefEl());
    const auto rel_geom = mh.GeometryInParent(1, *cell);
    const Eigen::MatrixXd eval_nodes = layout->EvaluationNodes();
    const Eigen::MatrixXd parent_eval_nodes = rel_geom->Global(eval_nodes);
    const auto parent = mh.ParentInfos(1, 0)[idx].parent_ptr;
    const auto val_fine = mf_fine(*cell, eval_nodes);
    const auto val_coarse = mf_coarse(*parent, parent_eval_nodes);
    for (int i = 0; i < val_fine.size(); ++i) {
      ASSERT_TRUE(std::fabs(val_fine[i] - val_coarse[i]) < 1e-10);
    }
  }
}

TEST(lf_intergridfe, LagrInterpNodes) {
  for (int selector = 0; selector < 9; ++selector) {
    // Initialize a coarse test mesh
    const auto mesh_coarse =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(selector);
    // Initialize the coarse fe spaces
    const auto fes_coarse_o1 =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_coarse);
    const auto fes_coarse_o2 =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh_coarse);
    // Test for all possible refinement patterns
    for (const auto ref_pat :
         {lf::refinement::rp_regular, lf::refinement::rp_barycentric}) {
      const auto ndofs_coarse_o1 = fes_coarse_o1->LocGlobMap().NumDofs();
      const auto ndofs_coarse_o2 = fes_coarse_o2->LocGlobMap().NumDofs();
      // Check interpolation of first order lagrangian fe
      for (unsigned long dof = 0; dof < ndofs_coarse_o1; ++dof) {
        Eigen::VectorXd dofs = Eigen::VectorXd::Zero(ndofs_coarse_o1);
        dofs[dof] = 1;
        check_lagr_interp_nodes<lf::uscalfe::FeSpaceLagrangeO1<double>,
                                lf::uscalfe::FeSpaceLagrangeO1<double>>(
            mesh_coarse, fes_coarse_o1, dofs, ref_pat);
        check_lagr_interp_nodes<lf::uscalfe::FeSpaceLagrangeO1<double>,
                                lf::uscalfe::FeSpaceLagrangeO2<double>>(
            mesh_coarse, fes_coarse_o1, dofs, ref_pat);
      }
      // Check interpolation of second order lagrangian fe
      for (unsigned long dof = 0; dof < ndofs_coarse_o2; ++dof) {
        Eigen::VectorXd dofs = Eigen::VectorXd::Zero(ndofs_coarse_o2);
        dofs[dof] = 1;
        check_lagr_interp_nodes<lf::uscalfe::FeSpaceLagrangeO2<double>,
                                lf::uscalfe::FeSpaceLagrangeO1<double>>(
            mesh_coarse, fes_coarse_o2, dofs, ref_pat);
        check_lagr_interp_nodes<lf::uscalfe::FeSpaceLagrangeO2<double>,
                                lf::uscalfe::FeSpaceLagrangeO2<double>>(
            mesh_coarse, fes_coarse_o2, dofs, ref_pat);
      }
    }
  }
}
