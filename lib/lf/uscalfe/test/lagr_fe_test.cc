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
#include <lf/uscalfe/lagrfe.h>
#include <iostream>

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/lagrfe.h>

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

TEST(lf_fe, scal_fe_coeff_node) {
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

TEST(lf_fe, scal_fe_val_node) {
  // Test of exactness of nodal reconstruction
  std::cout << ">>> Linear FE: test of consistency of nodal interpolation"
            << std::endl;
  FeLagrangeO1Tria<double> tlfe{};
  EXPECT_TRUE(scalarFEInterpTest(tlfe));
  FeLagrangeO1Quad<double> qlfe{};
  EXPECT_TRUE(scalarFEInterpTest(tlfe));
  FeLagrangeO1Segment<double> slfe{};
  EXPECT_TRUE(scalarFEInterpTest(tlfe));
}

TEST(lf_fe, lf_fe_linfe) {
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

TEST(lf_fe, lf_fe_segment) {
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

TEST(lf_fe, lf_fe_ellbvp) {
  std::cout << "### TEST: Computation of element matrices" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Set up objects taking care of local computations
  auto alpha = MeshFunctionGlobal([](Eigen::Vector2d) { return 1.0; });
  auto gamma = MeshFunctionGlobal([](Eigen::Vector2d) { return 0.0; });
  LagrangeFEEllBVPElementMatrix<double, decltype(alpha), decltype(gamma)>
      comp_elem_mat{fe_space, alpha, gamma};
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
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);
  std::shared_ptr<FeLagrangeO1Segment<double>> fe_p{
      std::make_shared<FeLagrangeO1Segment<double>>()};

  // Set up objects taking care of local computations
  auto gamma = MeshFunctionConstant(1.0);
  LagrangeFEEdgeMassMatrix comp_elem_mat(fe_space, gamma);
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
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Set up objects taking care of local computations
  auto f = MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2 * x[0] + x[1]); });
  using loc_comp_t = ScalarFELocalLoadVector<double, decltype(f)>;

  // Set debugging flags
  loc_comp_t::ctrl_ = 0;                                       // 255;
  lf::quad::QuadRule::out_ctrl_ = 0;                           // 1;
  LinearFELocalLoadVector<double, decltype(f)>::dbg_ctrl = 0;  // 3;

  // Instantiate object for local computations
  loc_comp_t comp_elem_vec(fe_space, f);

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

TEST(lf_fe, lf_fe_edgeload) {
  std::cout << "### TEST: Computation of local edge load vector" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Set up objects taking care of local computations
  auto g = MeshFunctionConstant(1.0);
  ScalarFEEdgeLocalLoadVector comp_elem_vec{fe_space, g};
  // Set debugging flags
  // comp_elem_mat.ctrl_ = 255;
  // lf::quad::QuadRule::out_ctrl_ = 1;

  Eigen::Vector2d Ref_vec(2);
  Ref_vec[0] = Ref_vec[1] = 0.5;

  // Loop over edges and compute element vectors
  for (const lf::mesh::Entity &edge : mesh_p->Entities(1)) {
    const lf::base::size_type n(edge.RefEl().NumNodes());
    std::cout << "Edge " << edge << ": [";
    typename decltype(comp_elem_vec)::ElemVec elem_vec{
        comp_elem_vec.Eval(edge)};
    std::cout << elem_vec.transpose() << " ]" << std::endl;
    const double diffnorm =
        (elem_vec - lf::geometry::Volume(*edge.Geometry()) * Ref_vec).norm();
    EXPECT_NEAR(diffnorm, 0.0, 1E-6);
  }
}

}  // end namespace lf::uscalfe::test
