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

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>

namespace lf::uscalfe::test {

template <typename SCALAR>
bool scalarFEEvalNodeTest(const ScalarReferenceFiniteElement<SCALAR> &fe_desc) {
  // Evaluates a random linear combination of reference shape functions
  // at the evaluation nodes for the finite element and then reconstructs
  // the "interpolant" which must agree with what we started from

  // Fetch evaluation nodes
  const Eigen::MatrixXd evl_nodes{fe_desc.EvaluationNodes()};
  const size_type N_evln = fe_desc.NumEvaluationNodes();
  EXPECT_EQ(evl_nodes.cols(), N_evln) << "No. evl nodes mismatch";

  // Evaluate reference shape functions in evaluation nodes
  const size_type N_rsf = fe_desc.NumRefShapeFunctions();
  auto rsf_at_evln = fe_desc.EvalReferenceShapeFunctions(evl_nodes);
  EXPECT_EQ(rsf_at_evln.rows(), N_rsf) << "No. rsf mismatch";

  // Form random linear combination and store it values
  // in the evaluation nodes
  Eigen::RowVectorXd rand_coeffs{Eigen::RowVectorXd::Random(N_rsf)};
  Eigen::RowVectorXd nodvals = rand_coeffs * rsf_at_evln;

  // Reconstruct linear combination
  Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> coeffs{
      fe_desc.NodalValuesToDofs(nodvals)};
  // Check agreement of coefficients
  EXPECT_DOUBLE_EQ((coeffs - rand_coeffs).norm(), 0.0)
      << "Coefficient mismatch" << coeffs << " <-> " << rand_coeffs;
  return true;
}

TEST(lf_fe_linear, scal_fe_coeff_node) {
  // Test of consistency of nodal interpolation
  std::cout << ">>> Linear FE: test of consistency of nodal interpolation"
            << std::endl;

  FeLagrangeO1Tria<double> tlfe{};
  std::cout << tlfe << std::endl;
  EXPECT_TRUE(scalarFEEvalNodeTest(tlfe));

  FeLagrangeO1Quad<double> qlfe{};
  std::cout << qlfe << std::endl;
  EXPECT_TRUE(scalarFEEvalNodeTest(tlfe));

  FeLagrangeO1Segment<double> slfe{};
  std::cout << slfe << std::endl;
  EXPECT_TRUE(scalarFEEvalNodeTest(tlfe));
}

TEST(lf_fe_quadratic, scal_fe_coeff_node) {
  // Test of consistency of nodal interpolation
  std::cout << ">>> Quadratic FE: test of consistency of nodal interpolation"
            << std::endl;

  FeLagrangeO2Tria<double> tfe{};
  std::cout << tfe << std::endl;
  EXPECT_TRUE(scalarFEEvalNodeTest(tfe));

  FeLagrangeO2Quad<double> qfe{};
  std::cout << qfe << std::endl;
  EXPECT_TRUE(scalarFEEvalNodeTest(qfe));

  FeLagrangeO2Segment<double> sfe{};
  std::cout << sfe << std::endl;
  EXPECT_TRUE(scalarFEEvalNodeTest(sfe));
}

template <typename SCALAR>
bool scalarFEInterpTest(const ScalarReferenceFiniteElement<SCALAR> &fe_desc) {
  // Interpolates random values at interpolation nodes
  // and checks whether the resulting linear combination of
  // basis functions reproduces those values

  // Fetch evaluation nodes
  const Eigen::MatrixXd evl_nodes{fe_desc.EvaluationNodes()};
  const size_type N_evln = fe_desc.NumEvaluationNodes();
  EXPECT_EQ(evl_nodes.cols(), N_evln) << "No. evl nodes mismatch";

  // Evaluate reference shape functions in evaluation nodes
  const size_type N_rsf = fe_desc.NumRefShapeFunctions();
  auto rsf_at_evln = fe_desc.EvalReferenceShapeFunctions(evl_nodes);
  EXPECT_EQ(rsf_at_evln.rows(), N_rsf) << "No. rsf mismatch";

  // Test makes sense only, if the number of local shape functions
  // agrees with the number of evaluation nodes
  EXPECT_EQ(N_evln, N_rsf) << "Nos of rsf and evaluation nodes must agree";

  // Vector of random nodal values
  Eigen::RowVectorXd rand_vals{Eigen::RowVectorXd::Random(N_evln)};
  // Obtain corresponding linear combination of shape functions
  Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> coeffs{
      fe_desc.NodalValuesToDofs(rand_vals)};

  // Evaluate linear combination of basis functions at evaluation nodes
  Eigen::RowVectorXd nodvals = coeffs * rsf_at_evln;

  // Check agreement of values
  EXPECT_DOUBLE_EQ((nodvals - rand_vals).norm(), 0.0)
      << "Value mismatch" << nodvals << " <-> " << rand_vals;
  return true;
}

TEST(lf_fe_linear, scal_fe_val_node) {
  // Test of exactness of nodal reconstruction
  std::cout << ">>> Linear FE: test of consistency of nodal interpolation"
            << std::endl;
  FeLagrangeO1Tria<double> tlfe{};
  EXPECT_TRUE(scalarFEInterpTest(tlfe));
  std::cout << tlfe << std::endl;

  FeLagrangeO1Quad<double> qlfe{};
  std::cout << qlfe << std::endl;
  EXPECT_TRUE(scalarFEInterpTest(tlfe));

  FeLagrangeO1Segment<double> slfe{};
  std::cout << slfe << std::endl;
  EXPECT_TRUE(scalarFEInterpTest(tlfe));
}

TEST(lf_fe_quadratic, scal_fe_val_node) {
  std::cout << ">>> Quadratic FE: test of consistency of nodal interpolation"
            << std::endl;

  FeLagrangeO2Tria<double> tfe{};
  std::cout << tfe << std::endl;
  EXPECT_TRUE(scalarFEInterpTest(tfe));

  FeLagrangeO2Quad<double> qfe{};
  std::cout << qfe << std::endl;
  EXPECT_TRUE(scalarFEInterpTest(qfe));

  FeLagrangeO2Segment<double> sfe{};
  std::cout << sfe << std::endl;
  EXPECT_TRUE(scalarFEInterpTest(sfe));
}

TEST(lf_fe_linear, lf_fe_linfe) {
  // Three points in the reference element
  Eigen::MatrixXd refcoords{
      (Eigen::MatrixXd(2, 3) << 0.3, 0.1, 0.7, 0.2, 0.5, 0.1).finished()};
  std::cout << "Points in reference cell\n" << refcoords << std::endl;

  // Testing triangular element
  {
    FeLagrangeO1Tria<double> tlfe{};
    EXPECT_EQ(tlfe.NumRefShapeFunctions(), 3);
    EXPECT_EQ(tlfe.NumRefShapeFunctions(0, 0), 0);
    EXPECT_EQ(tlfe.NumRefShapeFunctions(1, 0), 0);
    EXPECT_EQ(tlfe.NumRefShapeFunctions(2, 0), 1);

    auto rsf_vals = tlfe.EvalReferenceShapeFunctions(refcoords);
    EXPECT_EQ(rsf_vals.rows(), 3);
    EXPECT_EQ(rsf_vals.cols(), 3);
    EXPECT_DOUBLE_EQ(rsf_vals(0, 0), 0.5);
    EXPECT_DOUBLE_EQ(rsf_vals(1, 0), 0.3);
    EXPECT_DOUBLE_EQ(rsf_vals(2, 0), 0.2);
    EXPECT_DOUBLE_EQ(rsf_vals(0, 1), 0.4);
    EXPECT_DOUBLE_EQ(rsf_vals(1, 1), 0.1);
    EXPECT_DOUBLE_EQ(rsf_vals(2, 1), 0.5);
    EXPECT_DOUBLE_EQ(rsf_vals(0, 2), 0.2);
    EXPECT_DOUBLE_EQ(rsf_vals(1, 2), 0.7);
    EXPECT_DOUBLE_EQ(rsf_vals(2, 2), 0.1);

    EXPECT_TRUE(
        tlfe.EvalReferenceShapeFunctions(base::RefEl::kTria().NodeCoords())
            .isApprox(Eigen::MatrixXd::Identity(3, 3)));

    EXPECT_TRUE(
        tlfe.EvaluationNodes().isApprox(base::RefEl::kTria().NodeCoords()));
  }

  // Testing quadrilateral element
  {
    FeLagrangeO1Quad<double> qlfe{};
    EXPECT_EQ(qlfe.NumRefShapeFunctions(), 4);
    EXPECT_EQ(qlfe.NumRefShapeFunctions(0, 0), 0);
    EXPECT_EQ(qlfe.NumRefShapeFunctions(1, 0), 0);
    EXPECT_EQ(qlfe.NumRefShapeFunctions(2, 0), 1);

    auto rsf_vals = qlfe.EvalReferenceShapeFunctions(refcoords);
    EXPECT_EQ(rsf_vals.rows(), 4);
    EXPECT_EQ(rsf_vals.cols(), 3);
    EXPECT_DOUBLE_EQ(rsf_vals(0, 0), 0.56);
    EXPECT_DOUBLE_EQ(rsf_vals(1, 0), 0.24);
    EXPECT_DOUBLE_EQ(rsf_vals(2, 0), 0.06);
    EXPECT_DOUBLE_EQ(rsf_vals(3, 0), 0.14);
    EXPECT_DOUBLE_EQ(rsf_vals(0, 1), 0.45);
    EXPECT_DOUBLE_EQ(rsf_vals(1, 1), 0.05);
    EXPECT_DOUBLE_EQ(rsf_vals(2, 1), 0.05);
    EXPECT_DOUBLE_EQ(rsf_vals(3, 1), 0.45);
    EXPECT_DOUBLE_EQ(rsf_vals(0, 2), 0.27);
    EXPECT_DOUBLE_EQ(rsf_vals(1, 2), 0.63);
    EXPECT_DOUBLE_EQ(rsf_vals(2, 2), 0.07);
    EXPECT_DOUBLE_EQ(rsf_vals(3, 2), 0.03);

    EXPECT_TRUE(
        qlfe.EvalReferenceShapeFunctions(base::RefEl::kQuad().NodeCoords())
            .isApprox(Eigen::MatrixXd::Identity(4, 4)));

    EXPECT_TRUE(
        qlfe.EvaluationNodes().isApprox(base::RefEl::kQuad().NodeCoords()));

    std::cout << "Quad: RSF gradients:\n "
              << qlfe.GradientsReferenceShapeFunctions(refcoords) << std::endl;
  }
}

TEST(lf_fe_quadratic, lf_fe_quadrfe) {
  // triangular finite element:
  {
    FeLagrangeO2Tria<double> tfe{};
    EXPECT_EQ(tfe.NumRefShapeFunctions(), 6);
    EXPECT_EQ(tfe.NumRefShapeFunctions(0, 0), 0);
    EXPECT_EQ(tfe.NumRefShapeFunctions(1, 0), 1);
    EXPECT_EQ(tfe.NumRefShapeFunctions(2, 0), 1);

    // cardinal basis property
    EXPECT_TRUE(tfe.EvalReferenceShapeFunctions(tfe.EvaluationNodes())
                    .isApprox(Eigen::MatrixXd::Identity(6, 6)));
  }

  // quadrilateral finite element:
  {
    FeLagrangeO2Quad<double> qfe{};
    EXPECT_EQ(qfe.NumRefShapeFunctions(), 9);
    EXPECT_EQ(qfe.NumRefShapeFunctions(0, 0), 1);
    EXPECT_EQ(qfe.NumRefShapeFunctions(1, 0), 1);
    EXPECT_EQ(qfe.NumRefShapeFunctions(2, 0), 1);

    EXPECT_TRUE(qfe.EvalReferenceShapeFunctions(qfe.EvaluationNodes())
                    .isApprox(Eigen::MatrixXd::Identity(9, 9)));
  }
}

TEST(lf_fe_linear, lf_fe_segment) {
  // Three points in unit interval
  Eigen::MatrixXd refcoords{
      (Eigen::MatrixXd(1, 3) << 0.3, 0.1, 0.7).finished()};
  std::cout << "Points in reference cell\n" << refcoords << std::endl;

  FeLagrangeO1Segment<double> slfe{};
  EXPECT_EQ(slfe.NumRefShapeFunctions(), 2);
  EXPECT_EQ(slfe.NumRefShapeFunctions(0, 0), 0);
  EXPECT_EQ(slfe.NumRefShapeFunctions(1, 0), 1);

  auto rsf_vals = slfe.EvalReferenceShapeFunctions(refcoords);
  EXPECT_EQ(rsf_vals.rows(), 2);
  EXPECT_EQ(rsf_vals.cols(), 3);
  EXPECT_DOUBLE_EQ(rsf_vals(0, 0), 0.7);
  EXPECT_DOUBLE_EQ(rsf_vals(1, 0), 0.3);
  EXPECT_DOUBLE_EQ(rsf_vals(0, 1), 0.9);
  EXPECT_DOUBLE_EQ(rsf_vals(1, 1), 0.1);
  EXPECT_DOUBLE_EQ(rsf_vals(0, 2), 0.3);
  EXPECT_DOUBLE_EQ(rsf_vals(1, 2), 0.7);

  auto rsf_grads = slfe.GradientsReferenceShapeFunctions(refcoords);
  EXPECT_TRUE(rsf_grads.row(0).isConstant(-1.0));
  EXPECT_TRUE(rsf_grads.row(1).isConstant(1.0));

  EXPECT_TRUE(
      slfe.EvaluationNodes().isApprox(base::RefEl::kSegment().NodeCoords()));
}

TEST(lf_fe_quadratic, lf_fe_segment) {
  FeLagrangeO2Segment<double> sfe{};
  EXPECT_EQ(sfe.NumRefShapeFunctions(), 3);
  EXPECT_EQ(sfe.NumRefShapeFunctions(0, 0), 1);
  EXPECT_EQ(sfe.NumRefShapeFunctions(1, 0), 1);

  EXPECT_TRUE(sfe.EvalReferenceShapeFunctions(sfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(3, 3)));
}

TEST(lf_fe_linear, lf_fe_ellbvp) {
  std::cout << "### TEST: Computation of element matrices" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Set up objects taking care of local computations
  auto alpha = MeshFunctionGlobal([](Eigen::Vector2d) { return 1.0; });
  auto gamma = MeshFunctionGlobal([](Eigen::Vector2d) { return 0.0; });
  ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                         decltype(gamma)>
      comp_elem_mat{fe_space, alpha, gamma};
  // Set debugging flags
  // comp_elem_mat.ctrl_ = 255;
  // lf::quad::QuadRule::out_ctrl_ = 1;

  // For comparison
  LinearFELaplaceElementMatrix lfe_elem_mat{};

  // Loop over cells and compute element matrices;
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    const lf::base::size_type n(cell->RefEl().NumNodes());
    std::cout << "CELL " << *cell << ":" << std::endl;
    std::cout << "Element matrix from LinearFELaplaceElementMatrix:"
              << std::endl;
    LinearFELaplaceElementMatrix::ElemMat lfe_mat{lfe_elem_mat.Eval(*cell)};
    std::cout << lfe_mat << std::endl;
    std::cout << "Element matrix from ReactionDiffusionElementMatrixProvider:"
              << std::endl;
    typename decltype(comp_elem_mat)::ElemMat quad_mat{
        comp_elem_mat.Eval(*cell)};
    std::cout << quad_mat << std::endl;
    EXPECT_NEAR((lfe_mat.block(0, 0, n, n) - quad_mat).norm(), 0.0, 1E-2);
  }
}

// Test that the ReactionDiffusionElementMatrixProvider works as expected
// for tensor valued coefficients.
TEST(lf_fe_linear, ReactionDiffusionEMPTensor) {
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Set up objects taking care of local computations
  auto alpha = MeshFunctionGlobal([](Eigen::Vector2d x) {
    return (Eigen::Matrix2d() << 1, x[0], x[1], x[0] * x[1]).finished();
  });
  auto gamma =
      MeshFunctionGlobal([](Eigen::Vector2d x) { return x[0] * x[1]; });

  // specify quad rule ( default of degree 2 is not enough for a quadratic
  // gamma)
  quad_rule_collection_t quad_rules{
      {lf::base::RefEl::kTria(),
       lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 4)},
      {lf::base::RefEl::kQuad(),
       lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 4)}};

  ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                         decltype(gamma)>
      emp{fe_space, alpha, gamma, quad_rules};

  assemble::COOMatrix<double> matrix(fe_space->LocGlobMap().NumDofs(),
                                     fe_space->LocGlobMap().NumDofs());
  AssembleMatrixLocally(0, fe_space->LocGlobMap(), fe_space->LocGlobMap(), emp,
                        matrix);

  // project two linear functions onto the fespace:
  MeshFunctionGlobal a([](Eigen::Vector2d x) { return 1 + x[0] + 2 * x[1]; });
  MeshFunctionGlobal b([](Eigen::Vector2d x) { return 3 * x[0]; });
  auto a_vec = NodalProjection<double>(*fe_space, a);
  auto b_vec = NodalProjection<double>(*fe_space, b);

  auto product = (a_vec.transpose() * matrix.makeSparse() * b_vec).eval();
  EXPECT_NEAR(product(0, 0), 7911. / 8., 1.0E-2);
}

TEST(lf_fe_linear, lf_fe_edgemass) {
  std::cout << "### TEST: Computation of local edge" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);
  std::shared_ptr<FeLagrangeO1Segment<double>> fe_p{
      std::make_shared<FeLagrangeO1Segment<double>>()};

  // Set up objects taking care of local computations
  auto gamma = MeshFunctionConstant(1.0);
  MassEdgeMatrixProvider comp_elem_mat(fe_space, gamma);
  // Set debugging flags
  // comp_elem_mat.ctrl_ = 255;
  // lf::quad::QuadRule::out_ctrl_ = 1;

  // Reference mass matrix
  Eigen::Matrix2d RefM(
      (Eigen::Matrix2d() << 1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 3.0)
          .finished());

  // Loop over edges and compute element matrices;
  for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
    const lf::base::size_type n(edge->RefEl().NumNodes());
    std::cout << "Edge " << *edge << ":" << std::endl;
    typename decltype(comp_elem_mat)::ElemMat quad_mat{
        comp_elem_mat.Eval(*edge)};
    std::cout << quad_mat << std::endl;
    const double diffnorm =
        (quad_mat - lf::geometry::Volume(*edge->Geometry()) * RefM).norm();
    EXPECT_NEAR(diffnorm, 0.0, 1E-6);
  }
}

TEST(lf_fe_linear, lf_fe_loadvec) {
  std::cout << "### TEST: Computation of element vectors" << std::endl;
  // Building the test mesh: a purely triangular mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Set up objects taking care of local computations
  auto f = MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2 * x[0] + x[1]); });
  using loc_comp_t = ScalarLoadElementVectorProvider<double, decltype(f)>;

  // Set debugging flags
  loc_comp_t::ctrl_ = 0;                                       // 255;
  lf::quad::QuadRule::out_ctrl_ = 0;                           // 1;
  LinearFELocalLoadVector<double, decltype(f)>::dbg_ctrl = 0;  // 3;

  // Instantiate object for local computations
  loc_comp_t comp_elem_vec(fe_space, f);

  // For comparison
  LinearFELocalLoadVector<double, decltype(f)> lfe_elem_vec(f);

  // Loop over cells and compute element matrices;
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    const lf::base::size_type n(cell->RefEl().NumNodes());
    std::cout << "CELL " << *cell << ":" << std::endl;
    std::cout << "Element vector from LinearFELaplaceElementMatrix:"
              << std::endl;
    LinearFELocalLoadVector<double, decltype(f)>::ElemVec lfe_vec{
        lfe_elem_vec.Eval(*cell)};
    std::cout << "[ " << lfe_vec.transpose() << "] " << std::endl;
    std::cout << "Element vector from ScalarLoadElementVectorProvider:"
              << std::endl;
    loc_comp_t::ElemVec quad_vec{comp_elem_vec.Eval(*cell)};
    std::cout << "[ " << quad_vec.transpose() << "] " << std::endl;
    EXPECT_NEAR((lfe_vec.head(n) - quad_vec).norm(), 0.0, 1E-2);
  }
}

TEST(lf_fe_linear, lf_fe_edgeload) {
  std::cout << "### TEST: Computation of local edge load vector" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Set up objects taking care of local computations
  auto g = MeshFunctionConstant(1.0);
  ScalarLoadEdgeVectorProvider comp_elem_vec{fe_space, g};
  // Set debugging flags
  // comp_elem_mat.ctrl_ = 255;
  // lf::quad::QuadRule::out_ctrl_ = 1;

  Eigen::Vector2d Ref_vec(2);
  Ref_vec[0] = Ref_vec[1] = 0.5;

  // Loop over edges and compute element vectors
  for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
    const lf::base::size_type n(edge->RefEl().NumNodes());
    std::cout << "Edge " << *edge << ": [";
    typename decltype(comp_elem_vec)::ElemVec elem_vec{
        comp_elem_vec.Eval(*edge)};
    std::cout << elem_vec.transpose() << " ]" << std::endl;
    const double diffnorm =
        (elem_vec - lf::geometry::Volume(*edge->Geometry()) * Ref_vec).norm();
    EXPECT_NEAR(diffnorm, 0.0, 1E-6);
  }
}

TEST(lf_fe_quadratic, mass_mat_test) {
  std::cout << " >>> Quadratic FE: Test of computation of element matrices "
            << std::endl;

  // build test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // set up finite elements
  auto rfs_tria = std::make_shared<FeLagrangeO2Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO2Quad<double>>();
  auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
      mesh_p, rfs_tria, rfs_quad);

  // coefficients:
  auto alpha =
      MeshFunctionGlobal([](Eigen::Vector2d) -> double { return 1.0; });
  auto gamma =
      MeshFunctionGlobal([](Eigen::Vector2d) -> double { return 1.0; });

  // specify quad rule
  quad_rule_collection_t quad_rules{
      {lf::base::RefEl::kTria(),
       lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 4)},
      {lf::base::RefEl::kQuad(),
       lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 4)}};

  // set up object for local computations
  ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                         decltype(gamma)>
      provider(fe_space, alpha, gamma, quad_rules);

  // loop over cells and compute element matrices:
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    auto M = provider.Eval(*cell);

    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);
    const double vol = one_vec_r.dot(M * one_vec_c);

    EXPECT_NEAR(vol, lf::geometry::Volume(*cell->Geometry()), 1.0E-3)
        << "missmatch for cell " << cell << std::endl;
  }
}

TEST(lf_fe_quadratic, ReactionDiffusionEMPTensor) {
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // set up finite elements
  auto rfs_tria = std::make_shared<FeLagrangeO2Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO2Quad<double>>();
  auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
      mesh_p, rfs_tria, rfs_quad);

  // Set up objects taking care of local computations
  auto alpha = MeshFunctionGlobal([](Eigen::Vector2d x) {
    return (Eigen::Matrix2d() << 1, x[0], x[1], x[0] * x[1]).finished();
  });
  auto gamma =
      MeshFunctionGlobal([](Eigen::Vector2d x) { return x[0] * x[1]; });

  // specify quad rule ( default of degree 2 is not enough for a quadratic
  // gamma)
  quad_rule_collection_t quad_rules{
      {lf::base::RefEl::kTria(),
       lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 4)},
      {lf::base::RefEl::kQuad(),
       lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 4)}};

  ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                         decltype(gamma)>
      emp{fe_space, alpha, gamma, quad_rules};

  assemble::COOMatrix<double> matrix(fe_space->LocGlobMap().NumDofs(),
                                     fe_space->LocGlobMap().NumDofs());
  AssembleMatrixLocally(0, fe_space->LocGlobMap(), fe_space->LocGlobMap(), emp,
                        matrix);

  // project two linear functions onto the fespace:
  MeshFunctionGlobal a([](Eigen::Vector2d x) { return 1 + x[0] + 2 * x[1]; });
  MeshFunctionGlobal b([](Eigen::Vector2d x) { return 3 * x[0]; });
  auto a_vec = NodalProjection<double>(*fe_space, a);
  auto b_vec = NodalProjection<double>(*fe_space, b);

  auto product = (a_vec.transpose() * matrix.makeSparse() * b_vec).eval();
  EXPECT_NEAR(product(0, 0), 7911. / 8., 1E-2);
}

TEST(lf_fe_quadratic, lf_fe_edgemass) {
  std::cout << " >>> Quadratic FE: Test of computation of local edge matrices "
            << std::endl;

  // build test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // set up finite elements
  auto rfs_tria = std::make_shared<FeLagrangeO2Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO2Quad<double>>();
  auto rfs_segment = std::make_shared<FeLagrangeO2Segment<double>>();
  auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
      mesh_p, rfs_tria, rfs_quad, rfs_segment);

  // coefficientt
  auto gamma =
      MeshFunctionGlobal([](Eigen::Vector2d) -> double { return 1.0; });

  // Set up objects taking care of local computations
  auto g = MeshFunctionConstant(1.0);
  MassEdgeMatrixProvider provider{fe_space, g};

  // Loop over edges and compute element vectors
  for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
    auto M = provider.Eval(*edge);

    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);

    const double length = one_vec_r.dot(M * one_vec_c);
    EXPECT_NEAR(length, lf::geometry::Volume(*edge->Geometry()), 1.0E-3)
        << "missmatch for edge " << edge << std::endl;
  }
}

TEST(lf_fe_quadratic, lf_fe_loadvec) {
  std::cout << " >>> Quadratic FE: Test of computation of element load vectors "
            << std::endl;

  // build test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // set up finite elements
  auto rfs_tria = std::make_shared<FeLagrangeO2Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO2Quad<double>>();
  auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
      mesh_p, rfs_tria, rfs_quad);

  // coefficient:
  auto f =
      MeshFunctionGlobal([](Eigen::Vector2d x) -> double { return (1.0); });

  // set up object for local computations
  ScalarLoadElementVectorProvider<double, decltype(f)> provider(fe_space, f);

  // loop over cells and compute element vectors:
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    auto v = provider.Eval(*cell);

    Eigen::VectorXd one_vec = Eigen::VectorXd::Constant(v.size(), 1.0);

    const double vol = one_vec.dot(v);
    EXPECT_NEAR(vol, lf::geometry::Volume(*cell->Geometry()), 1.0E-3)
        << "missmatch for cell " << cell << std::endl;
  }
}

TEST(lf_fe_quadratic, lf_fe_edgeload) {
  std::cout << " >>> Quadratic FE: Test of computation of local edge vectors "
            << std::endl;

  // build test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // set up finite elements
  auto rfs_tria = std::make_shared<FeLagrangeO2Tria<double>>();
  auto rfs_quad = std::make_shared<FeLagrangeO2Quad<double>>();
  auto rfs_segment = std::make_shared<FeLagrangeO2Segment<double>>();
  auto fe_space = std::make_shared<UniformScalarFESpace<double>>(
      mesh_p, rfs_tria, rfs_quad, rfs_segment);

  // Set up objects taking care of local computations
  auto g = MeshFunctionConstant(1.0);
  ScalarLoadEdgeVectorProvider provider{fe_space, g};

  // Loop over edges and compute element vectors
  for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
    auto v = provider.Eval(*edge);

    Eigen::VectorXd one_vec = Eigen::VectorXd::Constant(v.size(), 1.0);

    const double length = one_vec.dot(v);
    EXPECT_NEAR(length, lf::geometry::Volume(*edge->Geometry()), 1.0E-3)
        << "missmatch for edge " << edge << std::endl;
  }
}

}  // end namespace lf::uscalfe::test
