/**
 * @file
 * @brief Check that the classes in fe_tools.h work as expected.
 * @author Raffael Casagrande
 * @date   2019-01-19 09:18:35
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

namespace lf::fe::test {

TEST(lf_fe_feTools, IntegrateMeshFunction) {
  // Obtain a hybrid mesh of the square [0,3]^2
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  // Write mesh data to file for visualization
  io::VtkWriter vtk_writer(mesh, "mesh.vtk");

  // scalar valued mesh function sin(\pi*x)*y
  auto mfScalar = mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
    return std::sin(base::kPi * x[0]) * x[1];
  });
  // Compute the integral with a very high-order quadrature rule
  auto intScalar = lf::fe::IntegrateMeshFunction(*mesh, mfScalar, 20);
  EXPECT_FLOAT_EQ(intScalar, 9 / base::kPi);

  // vector valued mesh function: return type is automatically deduced from that
  // of the lambda function
  auto mfVec = mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d& x) {
    return Eigen::Vector3d(x[0], x[0] * x[1], std::cos(base::kPi * x[0]));
  });
  // Integrate with high-order quadrature rule and get back a vector!
  auto intVec = lf::fe::IntegrateMeshFunction(*mesh, mfVec, 20);
  EXPECT_FLOAT_EQ(intVec[0], 27. / 2.);
  EXPECT_FLOAT_EQ(intVec[1], 81. / 4.);
  EXPECT_LT(intVec[2], 1e-10);

  // dynamic-matrix-valued mesh function: Even in this case return type
  // deduction works!
  auto mfMatrixDyn = mesh::utils::MeshFunctionGlobal(
      [](const Eigen::Vector2d& x) -> Eigen::MatrixXd {
        return (Eigen::MatrixXd(2, 3) << 1, x[0], x[1], x[0] * x[1],
                x[0] * x[0], x[1] * x[1])
            .finished();
      });
  // Numerical integration returns a matrix
  auto intMatrixDyn = lf::fe::IntegrateMeshFunction(*mesh, mfMatrixDyn, 10);
  EXPECT_FLOAT_EQ(intMatrixDyn(0, 0), 9.);
  EXPECT_FLOAT_EQ(intMatrixDyn(0, 1), 27. / 2.);
  EXPECT_FLOAT_EQ(intMatrixDyn(0, 2), 27. / 2.);
  EXPECT_FLOAT_EQ(intMatrixDyn(1, 0), 81. / 4.);
  EXPECT_FLOAT_EQ(intMatrixDyn(1, 1), 27.);
  EXPECT_FLOAT_EQ(intMatrixDyn(1, 2), 27.);
}

TEST(lf_fe_feTools, NodalProjection) {
  // First stage: obtain a test mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3.0);
  // Build the finite element space of piecewise linear finite element functions
  std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_p =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Second stage: prepare function to be interpolated
  // scalar valued function
  auto expfn = [](Eigen::Vector2d x) -> double {
    return std::exp(x[0] * x[1]);
  };
  // Build a mesh function representing the lambda function
  auto expmf{mesh::utils::MeshFunctionGlobal(expfn)};

  // Third stage: perform nodal interpolation
  Eigen::VectorXd mu = lf::fe::NodalProjection(*fe_space_p, expmf);

  // Fourth stage: validation
  const lf::mesh::Mesh& mesh{*(fe_space_p->Mesh())};
  const lf::assemble::DofHandler& dofh{fe_space_p->LocGlobMap()};
  LF_VERIFY_MSG(
      dofh.NumDofs() == mesh.NumEntities(2),
      "mismatch " << dofh.NumDofs() << " <-> " << mesh.NumEntities(2));
  // Traverse nodes of the mesh
  for (const lf::mesh::Entity* node : mesh.Entities(2)) {
    // Obtain location of node
    const lf::geometry::Geometry& node_geo{*(node->Geometry())};
    const Eigen::MatrixXd p{lf::geometry::Corners(node_geo)};
    LF_VERIFY_MSG((p.cols() == 1) && (p.rows() == 2),
                  "Wrong vertex matrix size");
    // Obtain number of dof associated with current node
    EXPECT_EQ(dofh.NumLocalDofs(*node), 1)
        << "Too many dofs for " << *node << std::endl;
    const auto dof_idx_array{dofh.GlobalDofIndices(*node)};
    const lf::base::glb_idx_t dof_idx = dof_idx_array[0];
    const double dof_val = mu[dof_idx];
    EXPECT_FLOAT_EQ(dof_val, expfn(p.col(0)))
        << "dof_val = " << dof_val << " <-> mu_val = " << mu[dof_idx];
  }
}

TEST(lf_fe_feTools, InitEssentialConditionFromFunction) {
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3.0);

  // create a random high order space with varying polynomial degrees
  fe::HierarchicScalarFESpace<double> fes(
      mesh_p, [&](const mesh::Entity& e) { return mesh_p->Index(e) + 1; });

  auto on_boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1);

  // select edges where we want to impose the boundary condition
  auto edge_sel = [&](const mesh::Entity& e) -> bool {
    EXPECT_TRUE(e.RefEl() == base::RefEl::kSegment());
    return on_boundary(e) && mesh_p->Index(e) < 9;
  };

  // mesh function that contains prescribed values:
  auto mf = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) {
    return std::sin(4 * x.x()) * std::cos(3 * x.y());
  });

  // test InitEssentialConditionFromFunction
  auto result = InitEssentialConditionFromFunction(fes, edge_sel, mf);

  // for reference make a nodal projection:
  auto nodal_p = NodalProjection(fes, mf);

  const auto& dofh = fes.LocGlobMap();
  ASSERT_EQ(result.size(), dofh.NumDofs());

  // afterwards constrained_dofs[i] should be the same as result[i].first
  std::vector<bool> constrained_dofs(result.size(), false);
  for (const auto* ep : mesh_p->Entities(1)) {
    auto dofs = dofh.GlobalDofIndices(*ep);
    if (edge_sel(*ep) == false) {
      continue;
    }
    for (auto dof : dofs) {
      constrained_dofs[dof] = true;
      ASSERT_NEAR(result[dof].second, nodal_p(dof), 1e-7);
    }
  }

  for (base::size_type i = 0; i < result.size(); ++i) {
    ASSERT_EQ(constrained_dofs[i], result[i].first);
  }
}

}  // namespace lf::fe::test
