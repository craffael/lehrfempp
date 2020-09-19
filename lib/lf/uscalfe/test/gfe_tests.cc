/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Unit tests for parametric finite element facilities
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/uscalfe/uscalfe.h>
#include <iostream>

#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include "fe_testutils.h"
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::uscalfe::test {

TEST(lf_gfe, lf_gfe_l2norm) {
  // LocCompLagrFEPreprocessor::ctrl_ = 255;
  std::cout << "### TEST Computation of L2 norm" << std::endl;
  // This test computes an approximation of the L2 norm of a function
  // by local quadrature on a finite element mesh, using the facilities
  // for computation of norms of differences of functions.

  // Building the test mesh: a purely affine mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // The function to integrate
  auto mf = mesh::utils::MeshFunctionGlobal(
      [](auto x) { return (x[0] * (1.0 - x[1])); });

  // compute norm by integrating the mesh function mf over the mesh
  double norm = std::sqrt(IntegrateMeshFunction(*mesh_p, squaredNorm(mf), 4));

  EXPECT_NEAR(norm, 5.19615, 1E-4);
}

TEST(lf_gfe, lf_gfe_L2assnorm) {
  // LocCompLagrFEPreprocessor::ctrl_ = 255;
  std::cout << "### TEST Compute L2 norm in two different ways" << std::endl;
  // This test computes the L2 norm of a FE function in two ways

  // Building the test mesh: a purely affine mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // Set up global FE space
  auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
      mesh_p, std::make_shared<FeLagrangeO1Tria<double>>(),
      std::make_shared<FeLagrangeO1Quad<double>>());

  // Local-to-global index map
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  // Dimension of finite element space`
  const lf::assemble::size_type N_dofs(dofh.NumDofs());
  // Matrix in triplet format holding Galerkin matrix
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Assemble finite element Galerkin matrix (mass matrix)
  uscalfe::test::SecOrdBVPLagrFEFullInteriorGalMat(
      fe_space, mesh::utils::MeshFunctionConstant(0.0),
      mesh::utils::MeshFunctionConstant(1.0), A);

  // First optin for setting the coefficient vector
  // Model linear function
  // auto u = [](auto x) { return (x[0] - 2 * x[1]); };
  // Interpolation
  // auto coeffvec{NodalProjection(fe_space, u)};

  // Alternative:
  // Random coefficient vector of FE function
  Eigen::VectorXd coeffvec{Eigen::VectorXd::Random(dofh.NumDofs())};

  double normsq_from_A = coeffvec.dot(A.MatVecMult(1.0, coeffvec));

  // Directly compute L2 norm of FE function by integrating the
  // associated mesh function:
  MeshFunctionFE<double, double> mf_fe(fe_space, coeffvec);

  double normsq_direct = IntegrateMeshFunction(*mesh_p, squaredNorm(mf_fe), 2);

  std::cout << "Norm^2 through A = " << normsq_from_A << std::endl;
  std::cout << "Norm^2 directly = " << normsq_direct << std::endl;
  EXPECT_NEAR(normsq_from_A, normsq_direct, 1E-4);
}

TEST(lf_gfe, lf_gfe_H1assnorm) {
  // LocCompLagrFEPreprocessor::ctrl_ = 255;
  std::cout << "### TEST Compute H1 seminorm in two different ways"
            << std::endl;
  // This test computes the H1 seminorm of a FE function in two ways

  // Building the test mesh: a purely affine mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // Set up global FE space
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Local-to-global index map
  auto &dofh = fe_space->LocGlobMap();

  // Dimension of finite element space`
  const lf::assemble::size_type N_dofs(dofh.NumDofs());
  // Matrix in triplet format holding Galerkin matrix
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Assemble finite element Galerkin matrix
  std::shared_ptr<UniformScalarFESpace<double>> temp = fe_space;
  SecOrdBVPLagrFEFullInteriorGalMat(temp,
                                    mesh::utils::MeshFunctionConstant(1.0),
                                    mesh::utils::MeshFunctionConstant(0.0), A);

  // First optin for setting the coefficient vector
  // Model linear function
  // auto u = [](auto x) { return (x[0] - 2 * x[1]); };
  // Interpolation
  // auto coeffvec{NodalProjection(fe_space, u)};

  // Alternative:
  // Random coefficient vector of FE function
  Eigen::VectorXd coeffvec{Eigen::VectorXd::Random(dofh.NumDofs())};

  const double normsq_from_A = coeffvec.dot(A.MatVecMult(1.0, coeffvec));

  // Directly compute H1 seminorm by integrating the mesh function
  MeshFunctionGradFE mf_grad_fe(fe_space, coeffvec);

  double normsq_direct =
      IntegrateMeshFunction(*mesh_p, squaredNorm(mf_grad_fe), 2);

  std::cout << "H1 norm^2 through A = " << normsq_from_A << std::endl;
  std::cout << "H1 norm^2 directly = " << normsq_direct << std::endl;
  EXPECT_NEAR(normsq_from_A, normsq_direct, 1E-4);
}

TEST(lf_gfe, lf_gfe_l2norm_vf) {
  // LocCompLagrFEPreprocessor::ctrl_ = 255;
  std::cout << "### TEST Computation of L2 norm of vectorfield" << std::endl;
  // This test computes an approximation of the L2 norm of a function
  // by local quadrature on a finite element mesh, using the facilities
  // for computation of norms of differences of functions.

  // Building the test mesh: a purely affine mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // mesh function to integrate:
  auto mf = mesh::utils::MeshFunctionGlobal(
      [](auto x) { return (Eigen::Vector2d() << x[1], -x[0]).finished(); });

  // integrate mesh function:
  double norm = std::sqrt(IntegrateMeshFunction(*mesh_p, squaredNorm(mf), 4));

  std::cout << "Norm of vectorfield = " << norm << std::endl;
  EXPECT_NEAR(norm, 7.34847, 1E-4);
}

TEST(lf_gfe, lf_gfe_lintp) {
  std::cout << "### TEST: Linear Interpolation" << std::endl;
  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Set up global FE space
  auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
      mesh_p, std::make_shared<FeLagrangeO1Tria<double>>(),
      std::make_shared<FeLagrangeO1Quad<double>>());

  auto u = [](auto x) { return std::exp(x[0] * (1.0 - x[1])); };
  auto coeffvec =
      NodalProjection(*fe_space, mesh::utils::MeshFunctionGlobal(u));

  // Local-to-global index mapped
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  EXPECT_EQ(coeffvec.size(), dofh.NumDofs())
      << "Coefficient vector length mismatch";

  // Check agreement of nodal values
  for (lf::assemble::gdof_idx_t j = 0; j < coeffvec.size(); ++j) {
    const lf::mesh::Entity &node{dofh.Entity(j)};
    EXPECT_TRUE(node.RefEl() == lf::base::RefEl::kPoint());
    auto tmp_coord(
        node.Geometry()->Global(lf::base::RefEl::kPoint().NodeCoords()));
    auto pt_coord(tmp_coord.col(0));
    // std::cout << "@ [" << pt_coord.transpose() << "]: u = " << u(pt_coord)
    //           << " <-> u_h = " << coeffvec[j] << std::endl;
    EXPECT_DOUBLE_EQ(coeffvec[j], u(pt_coord))
        << " @ [" << pt_coord.transpose() << "]:";
  }
}

TEST(lf_gfe, set_dirbdc) {
  std::cout << "****** TEST: Setting Dirichlet boundary values *****"
            << std::endl;
  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Set up global FE space
  auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
      mesh_p, std::make_shared<FeLagrangeO1Tria<double>>(),
      std::make_shared<FeLagrangeO1Quad<double>>(),
      std::make_shared<FeLagrangeO1Segment<double>>());
  // Specification of local shape functions for a edge
  auto fe_spec_edge_p{
      fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment())};
  LF_ASSERT_MSG(fe_spec_edge_p,
                "Reference Finite Element for segment is missing.");
  // Local to global index mapping
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  // Retrieve edges on the boundary
  lf::mesh::utils::CodimMeshDataSet<bool> bd_flags(
      lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1));
  // Selector for boundary edges
  auto bd_edge_sel = [&bd_flags](const lf::mesh::Entity &edge) -> bool {
    return bd_flags(edge);
  };
  // Function defining the Dirichlet boundary values
  auto g = [](Eigen::Vector2d x) -> double {
    return (x[0] * x[0] + x[1] * x[1]);
  };

  auto flag_val_vec(InitEssentialConditionFromFunction(
      dofh, *fe_spec_edge_p, bd_edge_sel, mesh::utils::MeshFunctionGlobal(g)));
  // Checking agreement of values
  for (lf::assemble::gdof_idx_t j = 0; j < flag_val_vec.size(); ++j) {
    const lf::mesh::Entity &node{dofh.Entity(j)};
    EXPECT_TRUE(node.RefEl() == lf::base::RefEl::kPoint());
    auto tmp_coord(
        node.Geometry()->Global(lf::base::RefEl::kPoint().NodeCoords()));
    auto pt_coord(tmp_coord.col(0));
    std::cout << "@ [" << pt_coord.transpose()
              << "]: bd = " << ((flag_val_vec[j].first) ? "TRUE" : "FALSE")
              << ", val = " << flag_val_vec[j].second
              << " <-> g = " << g(pt_coord) << std::endl;
    if (flag_val_vec[j].first) {
      EXPECT_DOUBLE_EQ(flag_val_vec[j].second, g(pt_coord))
          << " @ [" << pt_coord.transpose() << "]:";
    } else {
      EXPECT_DOUBLE_EQ(flag_val_vec[j].second, 0.0)
          << " @ [" << pt_coord.transpose() << "]:";
    }
  }
}

// check whether linear function is interpolated exactly
bool checkInterpolationLinear(
    const std::shared_ptr<UniformScalarFESpace<double>> &fe_space) {
  // Model linear function
  auto u =
      mesh::utils::MeshFunctionGlobal([](auto x) { return (x[0] - 2 * x[1]); });
  // Interpolation
  Eigen::VectorXd coeffvec = NodalProjection(*fe_space, u);
  auto mf_fe = MeshFunctionFE<double, double>(fe_space, coeffvec);
  auto norm =
      IntegrateMeshFunction(*fe_space->Mesh(), squaredNorm(mf_fe - u), 4);

  EXPECT_NEAR(norm, 0.0, 1E-6);
  return (std::fabs(norm) < 1.0E-6);
}

TEST(lf_gfe, lf_gfe_lintp_exact) {
  std::cout << "### TEST: Reproduction of linear functions" << std::endl;
  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Set up global FE space
  auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
      mesh_p, std::make_shared<FeLagrangeO1Tria<double>>(),
      std::make_shared<FeLagrangeO1Quad<double>>());
  EXPECT_TRUE(checkInterpolationLinear(fe_space));
}

TEST(lf_gfe, lf_gfe_intperrcvg) {
  // Four levels of refinement
  const int reflevels = 6;
  std::cout << "### TEST: Convergence of interpolation error" << std::endl;
  // Building the test mesh: a general hybrid mesh
  // This serves as the coarsest mesh of the hierarchy
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p,
                                                              reflevels);
  lf::refinement::MeshHierarchy &multi_mesh{*multi_mesh_p};

  // output of mesh hierarchy
  multi_mesh.PrintInfo(std::cout);

  // Function
  auto f = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); });
  auto grad_f =
      mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) -> Eigen::Vector2d {
        return (std::exp(x[0] * x[1]) *
                (Eigen::Vector2d() << x[1], x[0]).finished())
            .eval();
      });

  auto errs{InterpolationErrors(multi_mesh, f, grad_f,
                                std::make_shared<FeLagrangeO1Tria<double>>(),
                                std::make_shared<FeLagrangeO1Quad<double>>())};

  size_type L = errs.size();
  for (int l = 0; l < L; ++l) {
    std::cout << "Level" << l << ": L2 error = " << errs[l].first
              << ", H1 error = " << errs[l].second << std::endl;
  }
}

}  // namespace lf::uscalfe::test
