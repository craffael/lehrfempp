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
#include "lf/fe/fe_space.h"
#include "lf/fe/fe_testutils.h"
#include "lf/fe/fe_tools.h"
#include "lf/fe/lin_fe.h"
#include "lf/fe/loc_comp_ellbvp.h"

#include <lf/mesh/utils/utils.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace lf::fe::test {

TEST(lf_fe, lf_fe_linfe) {
  // Three points in the reference element
  Eigen::MatrixXd refcoords{
      (Eigen::MatrixXd(2, 3) << 0.3, 0.1, 0.7, 0.2, 0.5, 0.1).finished()};
  std::cout << "Points in reference cell\n" << refcoords << std::endl;

  // Testing triangular element
  {
    TriaLinearLagrangeFE<double> tlfe{};
    EXPECT_EQ(tlfe.NumRefShapeFunctions(), 3);
    EXPECT_EQ(tlfe.NumRefShapeFunctions(0, 0), 0);
    EXPECT_EQ(tlfe.NumRefShapeFunctions(1, 0), 0);
    EXPECT_EQ(tlfe.NumRefShapeFunctions(2, 0), 1);

    auto rsf_vals = tlfe.EvalReferenceShapeFunctions(refcoords);
    for (const auto &v : rsf_vals) {
      std::cout << "Tria: RSF values: " << v << std::endl;
    }
    for (const auto &v : tlfe.GradientsReferenceShapeFunctions(refcoords)) {
      std::cout << "Tria: RSF gradients:\n " << v << std::endl;
    }
    std::cout << "Tria: Evaluation nodes\n"
              << tlfe.EvaluationNodes() << std::endl;
  }

  // Testing quadrilateral element
  {
    QuadLinearLagrangeFE<double> qlfe{};
    EXPECT_EQ(qlfe.NumRefShapeFunctions(), 4);
    EXPECT_EQ(qlfe.NumRefShapeFunctions(0, 0), 0);
    EXPECT_EQ(qlfe.NumRefShapeFunctions(1, 0), 0);
    EXPECT_EQ(qlfe.NumRefShapeFunctions(2, 0), 1);

    auto rsf_vals = qlfe.EvalReferenceShapeFunctions(refcoords);
    for (const auto &v : rsf_vals) {
      std::cout << "Quad: RSF values: " << v << std::endl;
    }
    for (const auto &v : qlfe.GradientsReferenceShapeFunctions(refcoords)) {
      std::cout << "Quad: RSF gradients:\n " << v << std::endl;
    }
    std::cout << "Quad: Evaluation nodes\n"
              << qlfe.EvaluationNodes() << std::endl;
  }
}

TEST(lf_fe, lf_fe_segment) {
  // Three points in unit interval
  Eigen::MatrixXd refcoords{
      (Eigen::MatrixXd(1, 3) << 0.3, 0.1, 0.7).finished()};
  std::cout << "Points in reference cell\n" << refcoords << std::endl;

  SegmentLinearLagrangeFE<double> slfe{};
  EXPECT_EQ(slfe.NumRefShapeFunctions(), 2);
  EXPECT_EQ(slfe.NumRefShapeFunctions(0, 0), 0);
  EXPECT_EQ(slfe.NumRefShapeFunctions(1, 0), 1);

  auto rsf_vals = slfe.EvalReferenceShapeFunctions(refcoords);
  for (const auto &v : rsf_vals) {
    std::cout << "Segment: RSF values: " << v << std::endl;
  }
  for (const auto &v : slfe.GradientsReferenceShapeFunctions(refcoords)) {
    std::cout << "Segment: RSF gradients:\n " << v << std::endl;
  }
  std::cout << "Segment: Evaluation nodes\n"
            << slfe.EvaluationNodes() << std::endl;
}

TEST(lf_fe, lf_fe_ellbvp) {
  std::cout << "### TEST: Computation of element matrices" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  std::shared_ptr<TriaLinearLagrangeFE<double>> tlfe_p{
      std::make_shared<TriaLinearLagrangeFE<double>>()};
  std::shared_ptr<QuadLinearLagrangeFE<double>> qlfe_p{
      std::make_shared<QuadLinearLagrangeFE<double>>()};

  // Set up objects taking care of local computations
  auto alpha = [](Eigen::Vector2d) -> double { return 1.0; };
  auto gamma = [](Eigen::Vector2d) -> double { return 0.0; };
  LagrangeFEEllBVPElementMatrix<double, decltype(alpha), decltype(gamma)>
      comp_elem_mat{tlfe_p, qlfe_p, alpha, gamma};
  // Set debugging flags
  // comp_elem_mat.ctrl_ = 255;
  // lf::quad::QuadRule::out_ctrl_ = 1;

  // For comparison
  LinearFELaplaceElementMatrix lfe_elem_mat{};

  // Loop over cells and compute element matrices;
  for (const lf::mesh::Entity &cell : mesh_p->Entities(0)) {
    const lf::base::size_type n(cell.RefEl().NumNodes());
    std::cout << "CELL " << cell << ":" << std::endl;
    std::cout << "Element matrix from LinearFELaplaceElementMatrix:"
              << std::endl;
    LinearFELaplaceElementMatrix::ElemMat lfe_mat{lfe_elem_mat.Eval(cell)};
    std::cout << lfe_mat << std::endl;
    std::cout << "Element matrix from LagrangeFEEllBVPElementMatrix:"
              << std::endl;
    typename decltype(comp_elem_mat)::ElemMat quad_mat{
        comp_elem_mat.Eval(cell)};
    std::cout << quad_mat << std::endl;
    EXPECT_NEAR((lfe_mat.block(0, 0, n, n) - quad_mat).norm(), 0.0, 1E-2);
  }
}

TEST(lf_fe, lf_fe_edgemass) {
  std::cout << "### TEST: Computation of local edge" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  std::shared_ptr<SegmentLinearLagrangeFE<double>> fe_p{
      std::make_shared<SegmentLinearLagrangeFE<double>>()};

  // Set up objects taking care of local computations
  auto gamma = [](Eigen::Vector2d) -> double { return 1.0; };
  LagrangeFEEdgeMassMatrix<double, decltype(gamma)> comp_elem_mat{fe_p, gamma};
  // Set debugging flags
  // comp_elem_mat.ctrl_ = 255;
  // lf::quad::QuadRule::out_ctrl_ = 1;

  // Reference mass matrix
  Eigen::Matrix2d RefM(
      (Eigen::Matrix2d() << 1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0)
          .finished());

  // Loop over edges and compute element matrices;
  for (const lf::mesh::Entity &edge : mesh_p->Entities(1)) {
    const lf::base::size_type n(edge.RefEl().NumNodes());
    std::cout << "Edge " << edge << ":" << std::endl;
    typename decltype(comp_elem_mat)::ElemMat quad_mat{
        comp_elem_mat.Eval(edge)};
    std::cout << quad_mat << std::endl;
    const double diffnorm =
        (quad_mat - lf::geometry::Volume(*edge.Geometry()) * RefM).norm();
    EXPECT_NEAR(diffnorm, 0.0, 1E-6);
  }
}

TEST(lf_fe, lf_fe_loadvec) {
  std::cout << "### TEST: Computation of element vectors" << std::endl;
  // Building the test mesh: a purely triangular mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // Set up finite elements
  // Set up finite elements
  std::shared_ptr<TriaLinearLagrangeFE<double>> tlfe_p{
      std::make_shared<TriaLinearLagrangeFE<double>>()};
  std::shared_ptr<QuadLinearLagrangeFE<double>> qlfe_p{
      std::make_shared<QuadLinearLagrangeFE<double>>()};

  // Set up objects taking care of local computations
  auto f = [](Eigen::Vector2d x) -> double { return (2 * x[0] + x[1]); };
  using loc_comp_t = ScalarFELocalLoadVector<double, decltype(f)>;

  // Set debugging flags
  loc_comp_t::ctrl_ = 0;                                       // 255;
  lf::quad::QuadRule::out_ctrl_ = 0;                           // 1;
  LinearFELocalLoadVector<double, decltype(f)>::dbg_ctrl = 0;  // 3;

  // Instantiate object for local computations
  loc_comp_t comp_elem_vec(tlfe_p, qlfe_p, f);

  // For comparison
  LinearFELocalLoadVector<double, decltype(f)> lfe_elem_vec(f);

  // Loop over cells and compute element matrices;
  for (const lf::mesh::Entity &cell : mesh_p->Entities(0)) {
    const lf::base::size_type n(cell.RefEl().NumNodes());
    std::cout << "CELL " << cell << ":" << std::endl;
    std::cout << "Element vector from LinearFELaplaceElementMatrix:"
              << std::endl;
    LinearFELocalLoadVector<double, decltype(f)>::ElemVec lfe_vec{
        lfe_elem_vec.Eval(cell)};
    std::cout << "[ " << lfe_vec.transpose() << "] " << std::endl;
    std::cout << "Element vector from ScalarFELocalLoadVector:" << std::endl;
    loc_comp_t::ElemVec quad_vec{comp_elem_vec.Eval(cell)};
    std::cout << "[ " << quad_vec.transpose() << "] " << std::endl;
    EXPECT_NEAR((lfe_vec.head(n) - quad_vec).norm(), 0.0, 1E-2);
  }
}

TEST(lf_fe, lf_fe_l2norm) {
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

TEST(lf_fe, lf_fe_L2assnorm) {
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

TEST(lf_fe, lf_fe_H1assnorm) {
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

TEST(lf_fe, lf_fe_l2norm_vf) {
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

TEST(lf_fe, lf_fe_lintp) {
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

TEST(lf_fe, lf_fe_lintp_exact) {
  std::cout << "### TEST: Reproduction of linear functions" << std::endl;
  // Building the test mesh: a general hybrid mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);

  // Set up global FE space
  UniformScalarFiniteElementSpace fe_space(
      mesh_p, std::make_shared<TriaLinearLagrangeFE<double>>(),
      std::make_shared<QuadLinearLagrangeFE<double>>());
  EXPECT_TRUE(checkInterpolationLinear(fe_space));
}

TEST(lf_fe, lf_fe_intperrcvg) {
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

}  // end namespace lf::fe::test
