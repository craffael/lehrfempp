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
#include <iostream>
#include "lf/fe/fe_testutils.h"
#include "lf/fe/lin_fe.h"
#include "lf/fe/sec_ord_ell_bvp.h"

#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::fe::test {

TEST(lf_gfe, lf_gfe_l2norm) {
  // LocCompLagrFEPreprocessor::ctrl_ = 255;
  std::cout << "### TEST Computation of L2 norm" << std::endl;
  // This test computes an approximation of the L2 norm of a function
  // by local quadrature on a finite element mesh, using the facilities
  // for computation of norms of differences of functions.

  // Building the test mesh: a purely affine mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // Set up global FE space
  UniformScalarFiniteElementSpace fe_space(
      mesh_p, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>());

  // Dummy vector for coefficients of FE function
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
  std::vector<double> zerovec(dofh.NoDofs(), 0.0);

  // Helper object for computation of norm
  // Function passed as a generic lambda expression
  LocalL2NormDifference loc_comp(
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()),
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()),
      [](auto x) { return (x[0] * (1.0 - x[1])); }, 4);
  // loc_comp.ctrl_ = 255;

  /*
  // Alternative implementation
  // Function whose L2 should be computed; a generic lambda here
  auto u = [](auto x) { return (x[0] * (1.0 - x[1])); };
   LocalL2NormDifference<decltype(u)> loc_comp(
        fe_space.TriaShapeFunctionLayout(), fe_space.QuadShapeFunctionLayout(),
        [](auto x) { return (x[0] * (1.0 - x[1])); }, 4);
  */

  // Actual compuation of norm
  double norm = NormOfDifference(dofh, loc_comp, zerovec);
  std::cout << "Norm = " << norm << std::endl;
  EXPECT_NEAR(norm, 5.19615, 1E-4);
}

TEST(lf_gfe, lf_gfe_L2assnorm) {
  // LocCompLagrFEPreprocessor::ctrl_ = 255;
  std::cout << "### TEST Compute L2 norm in two different ways" << std::endl;
  // This test computes the L2 norm of a FE function in two ways

  // Building the test mesh: a purely affine mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // Set up global FE space
  UniformScalarFiniteElementSpace fe_space(
      mesh_p, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>());

  // Local-to-global index map
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

  // Dimension of finite element space`
  const lf::assemble::size_type N_dofs(dofh.NoDofs());
  // Matrix in triplet format holding Galerkin matrix
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Assemble finite element Galerkin matrix
  lf::fe::SecOrdBVPLagrFEFullInteriorGalMat(
      fe_space, [](auto x) -> double { return 0.0; },
      [](auto x) -> double { return 1.0; }, A);

  // First optin for setting the coefficient vector
  // Model linear function
  // auto u = [](auto x) { return (x[0] - 2 * x[1]); };
  // Interpolation
  // auto coeffvec{NodalProjection(fe_space, u)};

  // Alternative:
  // Random coefficient vector of FE function
  Eigen::VectorXd coeffvec{Eigen::VectorXd::Random(dofh.NoDofs())};

  const double normsq_from_A = coeffvec.dot(A.MatVecMult(1.0, coeffvec));

  // Directly compute L2 norm of FE function
  // Helper object for computation of norm
  // Function passed as a generic lambda expression
  LocalL2NormDifference loc_comp(
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()),
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()),
      [](auto x) -> double { return 0.0; });

  const double normsq_direct =
      std::pow(NormOfDifference(dofh, loc_comp, coeffvec), 2);

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
  LinearLagrangianFESpace fe_space(mesh_p);

  // Local-to-global index map
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

  // Dimension of finite element space`
  const lf::assemble::size_type N_dofs(dofh.NoDofs());
  // Matrix in triplet format holding Galerkin matrix
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Assemble finite element Galerkin matrix
  lf::fe::SecOrdBVPLagrFEFullInteriorGalMat(
      fe_space, [](auto x) -> double { return 1.0; },
      [](auto x) -> double { return 0.0; }, A);

  // First optin for setting the coefficient vector
  // Model linear function
  // auto u = [](auto x) { return (x[0] - 2 * x[1]); };
  // Interpolation
  // auto coeffvec{NodalProjection(fe_space, u)};

  // Alternative:
  // Random coefficient vector of FE function
  Eigen::VectorXd coeffvec{Eigen::VectorXd::Random(dofh.NoDofs())};

  const double normsq_from_A = coeffvec.dot(A.MatVecMult(1.0, coeffvec));

  // Directly compute H1 seminorm of FE function
  // Helper object for computation of norm
  LocL2GradientFEDifference loc_comp(
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()),
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()),
      [](auto x) { return Eigen::Vector2d::Zero(2); });

  const double normsq_direct =
      std::pow(NormOfDifference(dofh, loc_comp, coeffvec), 2);

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

  // Set up global FE space
  UniformScalarFiniteElementSpace fe_space(
      mesh_p, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>());

  // Dummy vector for coefficients of FE function
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
  std::vector<double> zerovec(dofh.NoDofs(), 0.0);

  // Helper object for computation of norm
  // Function passed as a generic lambda expression
  LocL2GradientFEDifference loc_comp(
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()),
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()),
      [](auto x) { return (Eigen::Vector2d() << x[1], -x[0]).finished(); }, 4);
  // loc_comp.ctrl_ = 255;

  // Actual compuation of norm
  double norm = NormOfDifference(dofh, loc_comp, zerovec);
  std::cout << "Norm of vectorfield = " << norm << std::endl;
  EXPECT_NEAR(norm, 7.34847, 1E-4);
}

TEST(lf_gfe, lf_gfe_lintp) {
  std::cout << "### TEST: Linear Interpolation" << std::endl;
  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Set up global FE space
  UniformScalarFiniteElementSpace fe_space(
      mesh_p, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>());

  auto u = [](auto x) { return std::exp(x[0] * (1.0 - x[1])); };
  auto coeffvec{NodalProjection(fe_space, u)};

  // Local-to-global index mapped
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};
  EXPECT_EQ(coeffvec.size(), dofh.NoDofs())
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
  UniformScalarFiniteElementSpace fe_space(
      mesh_p, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>(),
      std::make_shared<SegmentLinearLagrangeFE<double>>());
  // Specification of local shape functions for a edge
  auto fe_spec_edge_p{
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kSegment())};
  // Local to global index mapping
  const lf::assemble::DofHandler &dofh{fe_space.LocGlobMap()};

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

  auto flag_val_vec(InitEssentialConditionFromFunction(dofh, *fe_spec_edge_p,
                                                       bd_edge_sel, g));
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
bool checkInterpolationLinear(const UniformScalarFiniteElementSpace &fe_space) {
  // Model linear function
  auto u = [](auto x) { return (x[0] - 2 * x[1]); };
  // Interpolation
  auto coeffvec{NodalProjection(fe_space, u)};
  // Helper class for error computation
  LocalL2NormDifference loc_comp(
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kTria()),
      fe_space.ShapeFunctionLayout(lf::base::RefEl::kQuad()), u, 4);
  // Actual compuation of error norm
  double norm = NormOfDifference(fe_space.LocGlobMap(), loc_comp, coeffvec);
  std::cout << "Norm = " << norm << std::endl;
  EXPECT_NEAR(norm, 0.0, 1E-6);
  return (std::fabs(norm) < 1.0E-6);
}

TEST(lf_gfe, lf_gfe_lintp_exact) {
  std::cout << "### TEST: Reproduction of linear functions" << std::endl;
  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Set up global FE space
  UniformScalarFiniteElementSpace fe_space(
      mesh_p, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>());
  EXPECT_TRUE(checkInterpolationLinear(fe_space));
}

TEST(lf_gfe, lf_gfe_intperrcvg) {
  // Four levels of refinement
  const int reflevels = 6;
  std::cout << "### TEST: Convergence of interpolation error" << std::endl;
  // Building the test mesh: a general hybrid mesh
  // This serves as the coarsest mesh of the hierarchy
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Prepare for creating a hierarchy of meshes
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(std::move(mesh_p), mesh_factory_ptr);

  // Perform several steps of regular refinement of the given mesh
  for (int refstep = 0; refstep < reflevels; ++refstep) {
    // Barycentric refinement is the other option
    multi_mesh.RefineRegular(/*lf::refinement::RefPat::rp_barycentric*/);
  }

  // Function
  auto f = [](Eigen::Vector2d x) -> double { return std::exp(x[0] * x[1]); };
  auto grad_f = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return (std::exp(x[0] * x[1]) *
            (Eigen::Vector2d() << x[1], x[0]).finished())
        .eval();
  };

  auto errs{InterpolationErrors(
      multi_mesh, f, grad_f, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>())};

  size_type L = errs.size();
  for (int l = 0; l < L; ++l) {
    std::cout << "Level" << l << ": L2 error = " << errs[l].first
              << ", H1 error = " << errs[l].second << std::endl;
  }
}

TEST(lf_gfe, Neu_BVP_ass) {
  std::cout << "#### TEST: Assembly for Neumann BVP ###" << std::endl;
  // Right-hand side source function
  auto f = [](Eigen::Vector2d x) -> double {
    return (std::sin(2 * M_PI * x[0]) * std::sin(2 * M_PI * x[1]));
  };
  // Neumann data
  auto h = [](Eigen::Vector2d /*x*/) -> double { return 1.0; };

  // Class describing the Neumann boundary value problem
  std::shared_ptr<PureNeumannProblemLaplacian<decltype(f), decltype(h)>> bvp_p =
      std::make_shared<PureNeumannProblemLaplacian<decltype(f), decltype(h)>>(
          f, h);

  // Obtain finite Element space
  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Set up global FE space
  UniformScalarFiniteElementSpace fe_space(
      mesh_p, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>());

  auto [A, phi] = SecOrdEllBVPLagrFELinSys<double>(fe_space, bvp_p);
}

}  // namespace lf::fe::test
