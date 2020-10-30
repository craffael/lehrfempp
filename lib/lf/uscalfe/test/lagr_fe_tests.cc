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
#include <iostream>

#include <gtest/gtest.h>

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

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
  EXPECT_NEAR((coeffs - rand_coeffs).norm(), 0.0, 1e-13)
      << "Coefficient mismatch" << coeffs << " <-> " << rand_coeffs;
  return true;
}

TEST(lf_fe_linear, scal_fe_coeff_node) {
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

TEST(lf_fe_cubic, scalf_fe_coeff_node) {
  std::cout << ">>> Cubuic FE: test of consistency of nodal interpolation"
            << std::endl;
  FeLagrangeO3Tria<double> tfe{};
  std::cout << tfe << std::endl;
  EXPECT_TRUE(scalarFEEvalNodeTest(tfe));

  FeLagrangeO3Quad<double> qfe{};
  std::cout << qfe << std::endl;
  EXPECT_TRUE(scalarFEEvalNodeTest(qfe));

  FeLagrangeO3Segment<double> sfe{};
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
  EXPECT_NEAR((nodvals - rand_vals).norm(), 0.0, 1e-13)
      << "Value mismatch" << nodvals << " <-> " << rand_vals;
  return true;
}

TEST(lf_fe_linear, scal_fe_val_node) {
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

TEST(lf_fe_cubic, scal_fe_val_node) {
  std::cout << ">>> Cubic FE: test of consistency of nodal interpolation"
            << std::endl;

  FeLagrangeO3Tria<double> tfe{};
  std::cout << tfe << std::endl;
  EXPECT_TRUE(scalarFEInterpTest(tfe));

  FeLagrangeO3Quad<double> qfe{};
  std::cout << qfe << std::endl;
  EXPECT_TRUE(scalarFEInterpTest(qfe));

  FeLagrangeO3Segment<double> sfe{};
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

    // cardinal basis property
    EXPECT_TRUE(qfe.EvalReferenceShapeFunctions(qfe.EvaluationNodes())
                    .isApprox(Eigen::MatrixXd::Identity(9, 9)));
  }
}

TEST(lf_fe_cubic, lf_fe_cubfe) {
  // triangular finite element:
  {
    FeLagrangeO3Tria<double> tfe{};
    EXPECT_EQ(tfe.NumRefShapeFunctions(), 10);
    EXPECT_EQ(tfe.NumRefShapeFunctions(0, 0), 1);
    EXPECT_EQ(tfe.NumRefShapeFunctions(1, 0), 2);
    EXPECT_EQ(tfe.NumRefShapeFunctions(2, 0), 1);

    // cardinal basis property
    EXPECT_TRUE(tfe.EvalReferenceShapeFunctions(tfe.EvaluationNodes())
                    .isApprox(Eigen::MatrixXd::Identity(10, 10)));
  }

  // quadrilateral finite element:
  {
    FeLagrangeO3Quad<double> qfe{};
    EXPECT_EQ(qfe.NumRefShapeFunctions(), 16);
    EXPECT_EQ(qfe.NumRefShapeFunctions(0, 0), 4);
    EXPECT_EQ(qfe.NumRefShapeFunctions(1, 0), 2);
    EXPECT_EQ(qfe.NumRefShapeFunctions(2, 0), 1);

    // cardinal basis property
    EXPECT_TRUE(qfe.EvalReferenceShapeFunctions(qfe.EvaluationNodes())
                    .isApprox(Eigen::MatrixXd::Identity(16, 16)));
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

  // cardinal basis property
  EXPECT_TRUE(sfe.EvalReferenceShapeFunctions(sfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(3, 3)));
}

TEST(lf_fe_cubic, lf_fe_segment) {
  FeLagrangeO3Segment<double> sfe{};
  EXPECT_EQ(sfe.NumRefShapeFunctions(), 4);
  EXPECT_EQ(sfe.NumRefShapeFunctions(0, 0), 2);
  EXPECT_EQ(sfe.NumRefShapeFunctions(1, 0), 1);

  // cardinal basis property
  EXPECT_TRUE(sfe.EvalReferenceShapeFunctions(sfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(4, 4)));
}

/**
 * @brief computes a projection of the function g into the fe space and
 * retruns the squared H1 norm of the difference.
 *
 * @tparam SCALAR field like double
 * @tparam FUNCTION \ref mesh_function "MeshFunction" returning a scalar
 * @tparam FUNCTION_GRAD \ref mesh_function "MeshFunction" returning a
 * Eigen::Vector of scalars
 *
 *
 * @param fe_space  FE space onto which the function is projected
 * @param g  function which is projected onto the fe space
 * @param grad_g gradient of the function g
 * @param quad_degre degree of the quadrature rule used to approximate the
 * norm of the error
 */
template <typename SCALAR, typename FUNCTION, typename FUNCTION_GRAD>
SCALAR nodalProjectionTest(
    std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space, FUNCTION g,
    FUNCTION_GRAD grad_g, int quad_degree) {
  auto dof_vector = NodalProjection(*fe_space, g);
  auto mf_fe = MeshFunctionFE<double, double>(fe_space, dof_vector);
  auto grad_mf_fe = MeshFunctionGradFE<double, double>(fe_space, dof_vector);

  return IntegrateMeshFunction(
      *(fe_space->Mesh()),
      squaredNorm(g - mf_fe) + squaredNorm(grad_mf_fe - grad_g), quad_degree);
}

TEST(lf_fe_linear, projection_test) {
  std::cout << ">>> Linear FE: Projection into linear fe space" << std::endl;

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  std::shared_ptr<const UniformScalarFESpace<double>> fe_space =
      std::make_shared<const FeSpaceLagrangeO1<double>>(mesh_p);

  auto f_1 = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return x[0] + x[1]; });
  auto grad_1 = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d /*x*/) {
    return (Eigen::VectorXd(2) << 1.0, 1.0).finished();
  });

  auto f_2 = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return 5 * x[0] - 2 * x[1]; });
  auto grad_2 = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d /*x*/) {
    return (Eigen::VectorXd(2) << 5.0, -2.0).finished();
  });

  EXPECT_NEAR(nodalProjectionTest(fe_space, f_1, grad_1, 20), 0.0, 1e-10)
      << "projection error for f(x,y) = x + y";
  EXPECT_NEAR(nodalProjectionTest(fe_space, f_2, grad_2, 20), 0.0, 1e-10)
      << "projection error for f(x,y) = 5x - 2y";
}

TEST(lf_fe_quadratic, projection_test) {
  std::cout << ">>> Quadratic FE: Projection into quadratic fe space"
            << std::endl;

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  std::shared_ptr<const UniformScalarFESpace<double>> fe_space =
      std::make_shared<const FeSpaceLagrangeO2<double>>(mesh_p);

  auto f_1 = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return x[0] + x[1]; });
  auto grad_1 = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d /*x*/) {
    return (Eigen::VectorXd(2) << 1.0, 1.0).finished();
  });

  auto f_2 = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return 5 * x[0] * x[0] + 2 * x[1] * x[1]; });
  auto grad_2 = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) {
    return (Eigen::VectorXd(2) << 10 * x[0], 4 * x[1]).finished();
  });

  auto f_3 = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return x[0] * x[1]; });
  auto grad_3 = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) {
    return (Eigen::VectorXd(2) << x[1], x[0]).finished();
  });

  EXPECT_NEAR(nodalProjectionTest(fe_space, f_1, grad_1, 20), 0.0, 1e-10)
      << "projection error for f(x,y) = x + y";
  EXPECT_NEAR(nodalProjectionTest(fe_space, f_2, grad_2, 20), 0.0, 1e-10)
      << "projection error for f(x,y) = 5x^2 + 2y^2";
  EXPECT_NEAR(nodalProjectionTest(fe_space, f_3, grad_3, 20), 0.0, 1e-10)
      << "projection error for f(x,y) = x*y";
}

TEST(lf_fe_cubic, projection_test) {
  std::cout << ">>> Cubic FE: Projection into cubic fe space" << std::endl;

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  std::shared_ptr<const UniformScalarFESpace<double>> fe_space =
      std::make_shared<const FeSpaceLagrangeO3<double>>(mesh_p);

  auto f_1 = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return x[0] + x[1]; });
  auto grad_1 = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d /*x*/) {
    return (Eigen::VectorXd(2) << 1.0, 1.0).finished();
  });

  auto f_2 = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) {
    return 5 * x[0] * x[0] * x[0] + 2 * x[1] * x[1] * x[1];
  });
  auto grad_2 = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) {
    return (Eigen::VectorXd(2) << 15 * x[0] * x[0], 6 * x[1] * x[1]).finished();
  });

  auto f_3 = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) {
    return x[0] * x[0] * x[1] + 2 * x[0] * x[1] * x[1];
  });
  auto grad_3 = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) {
    return (Eigen::VectorXd(2) << 2 * x[0] * x[1] + 2 * x[1] * x[1],
            x[0] * x[0] + 4 * x[0] * x[1])
        .finished();
  });
  EXPECT_NEAR(nodalProjectionTest(fe_space, f_1, grad_1, 20), 0.0, 1e-10)
      << "projection error for f(x,y) = x + y";
  EXPECT_NEAR(nodalProjectionTest(fe_space, f_2, grad_2, 20), 0.0, 1e-10)
      << "projection error for f(x,y) = 5x^3 + 2y^3";
  EXPECT_NEAR(nodalProjectionTest(fe_space, f_3, grad_3, 20), 0.0, 1e-10)
      << "projection error for f(x,y) = x^2*y + 2*x*y^2";
}

TEST(lf_fe_linear, lf_fe_ellbvp) {
  std::cout << "### TEST: Computation of element matrices" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Set up objects taking care of local computations
  auto alpha =
      lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d) { return 1.0; });
  auto gamma =
      lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d) { return 0.0; });
  ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                         decltype(gamma)>
      comp_elem_mat{fe_space, alpha, gamma};

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

TEST(lf_fe_linear, lf_fe_edgemass) {
  std::cout << "### TEST: Computation of local edge" << std::endl;
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);
  std::shared_ptr<FeLagrangeO1Segment<double>> fe_p{
      std::make_shared<FeLagrangeO1Segment<double>>()};

  // Set up objects taking care of local computations
  auto gamma = lf::mesh::utils::MeshFunctionConstant(1.0);
  MassEdgeMatrixProvider comp_elem_mat(fe_space, gamma);

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
  auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return (2 * x[0] + x[1]); });
  using loc_comp_t = ScalarLoadElementVectorProvider<double, decltype(f)>;

  // Set debugging flags
  scalar_load_element_vector_provider_logger->set_level(spdlog::level::info);
  linear_fe_local_load_vector_logger->set_level(spdlog::level::info);  // 3;

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
  auto g = lf::mesh::utils::MeshFunctionConstant(1.0);
  ScalarLoadEdgeVectorProvider comp_elem_vec{fe_space, g};

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

/**
 * @brief approximates the value of  \f$ b(a_{FE},b_{FE}) \f$, where
 * \f$b(\dot,\dot)\f$ is the bilinear form described in
 * ReactionDiffusionElementMatrixProvider, and \f$ a_{FE}, b_{FE} \f$ are the
 * projections of a and b onto a finite element space.
 *
 * @tparam SCALAR field type such as double
 * @tparam MF_ALPHA a \ref mesh_function "MeshFunction" that defines the
 * diffusion coefficient in the bilinear form
 * @tparam MF_ALPHA a \ref mesh_function "MeshFunction" that defines the
 * reaction coefficient in the bilinear form
 * @tparam MF_A a SCALAR-valued \ref mesh_function "MeshFunction"
 * @tparam MF_A a SCALAR-valued \ref mesh_function "MeshFunction"
 *
 * @param fe_space finite element space onto which a and b are projected
 * @param alpha diffusion coefficient of the bilinear form
 * @param gamma reaction coefficient of the bilinear form
 * @param a function, whose projection onto the fe space is the first argument
 * of the evaluated bilinear form
 * @param b function, whose projection onto the fe sapce is the second argument
 * of the evalauted bilinear form
 * @param quad_degree degree of the quadrature rules used in local computations.
 */
template <typename SCALAR, typename MF_ALPHA, typename MF_GAMMA, typename MF_A,
          typename MF_B>
SCALAR reactionDiffusionTest(
    std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space,
    MF_ALPHA alpha, MF_GAMMA gamma, MF_A a, MF_B b, int quad_degree) {
  // construct quad rules:
  quad_rule_collection_t quad_rules{
      {lf::base::RefEl::kTria(),
       lf::quad::make_QuadRule(lf::base::RefEl::kTria(), quad_degree)},
      {lf::base::RefEl::kQuad(),
       lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), quad_degree)}};

  ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                         decltype(gamma)>
      provider{fe_space, alpha, gamma, quad_rules};

  // assemble system matrix
  assemble::COOMatrix<double> matrix(fe_space->LocGlobMap().NumDofs(),
                                     fe_space->LocGlobMap().NumDofs());
  AssembleMatrixLocally(0, fe_space->LocGlobMap(), fe_space->LocGlobMap(),
                        provider, matrix);

  // project functions onto the fe space
  auto a_vec = NodalProjection<double>(*fe_space, a);
  auto b_vec = NodalProjection<double>(*fe_space, b);

  // evaluate bilinear form on projected functions:
  auto product = (a_vec.transpose() * matrix.makeSparse() * b_vec).eval();
  return product(0, 0);
}

// Test that the ReactionDiffusionElementMatrixProvider works as expected
// for tensor valued coefficients.
TEST(lf_fe_linear, ReactionDiffusion) {
  std::cout << "Linear FE >>> Computation of bilinear forms" << std::endl;

  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  std::shared_ptr<const UniformScalarFESpace<double>> fe_space =
      std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Set parameter functions
  auto alpha = lf::mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) {
    return (Eigen::Matrix2d() << 1, x[0], x[1], x[0] * x[1]).finished();
  });
  auto gamma = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return x[0] * x[1]; });

  // functions projected into the fe space
  lf::mesh::utils::MeshFunctionGlobal a(
      [](Eigen::Vector2d x) { return 1 + x[0] + 2 * x[1]; });
  lf::mesh::utils::MeshFunctionGlobal b(
      [](Eigen::Vector2d x) { return 3 * x[0]; });

  auto product = reactionDiffusionTest(fe_space, alpha, gamma, a, b, 4);
  EXPECT_NEAR(product, 7911. / 8., 1.0E-2);
}

TEST(lf_fe_quadratic, ReactionDiffusion) {
  std::cout << "Quadratic FE >>> Computation of bilinear forms" << std::endl;

  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Set up finite elements
  std::shared_ptr<const UniformScalarFESpace<double>> fe_space =
      std::make_shared<FeSpaceLagrangeO2<double>>(mesh_p);

  // Set parameter functions
  auto alpha = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return x[0]; });
  auto gamma = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return x[0] * x[1]; });

  // functions projected into the fe space
  lf::mesh::utils::MeshFunctionGlobal a(
      [](Eigen::Vector2d x) { return x[0] * x[0] + x[1] * x[1]; });
  lf::mesh::utils::MeshFunctionGlobal b(
      [](Eigen::Vector2d x) { return x[0] * x[0] - x[1] * x[1]; });

  auto product = reactionDiffusionTest(fe_space, alpha, gamma, a, b, 6);
  EXPECT_NEAR(product, 81., 1.0E-2);
}

TEST(lf_fe_cubic, ReactionDiffusion) {
  std::cout << "Cubic FE >>> Computation of bilinear forms" << std::endl;

  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  // Set up finite elements
  std::shared_ptr<const UniformScalarFESpace<double>> fe_space =
      std::make_shared<FeSpaceLagrangeO3<double>>(mesh_p);

  // Set parameter functions
  auto alpha = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return x[1]; });
  auto gamma = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) { return x[0] * x[1]; });

  // functions projected into the fe space
  lf::mesh::utils::MeshFunctionGlobal a([](Eigen::Vector2d x) {
    return x[0] * x[0] * x[0] + x[1] * x[1] * x[1];
  });
  lf::mesh::utils::MeshFunctionGlobal b(
      [](Eigen::Vector2d x) { return x[0] * x[1] * x[1]; });

  auto product = reactionDiffusionTest(fe_space, alpha, gamma, a, b, 8);
  EXPECT_NEAR(product, 1996731 / 280., 1.0E-2);
}

/**
 * @brief checks products of the form \f$ x^T *M * y \f$ for certain
 * element/edge matrices/vectors \f$ M \f$ and vectors \f$ x,y \f$.
 * @param fe_space finite element space used for the computation of the element
 * matrices.
 */
template <typename SCALAR>
void locCompProductsTest(
    std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space) {
  auto mesh_p = fe_space->Mesh();

  // verification based on local element matrices
  {
    auto alpha = lf::mesh::utils::MeshFunctionConstant(1.0);
    auto gamma = lf::mesh::utils::MeshFunctionConstant(1.0);
    ReactionDiffusionElementMatrixProvider provider(fe_space, alpha, gamma);

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

  // verification based on local edge matrices
  {
    auto g = lf::mesh::utils::MeshFunctionConstant(1.0);
    MassEdgeMatrixProvider provider(fe_space, g);

    // loop over edges and compute element matrices:
    for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
      auto M = provider.Eval(*edge);

      Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
      Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);

      const double length = one_vec_r.dot(M * one_vec_c);
      EXPECT_NEAR(length, lf::geometry::Volume(*edge->Geometry()), 1.0E-3)
          << "missmatch for edge " << edge << std::endl;
    }
  }

  // verification based on local element load vectors:
  {
    auto f = lf::mesh::utils::MeshFunctionConstant(1.0);
    ScalarLoadElementVectorProvider provider(fe_space, f);

    // loop over cells and compute element vectors:
    for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
      auto v = provider.Eval(*cell);

      Eigen::VectorXd one_vec = Eigen::VectorXd::Constant(v.size(), 1.0);

      const double vol = one_vec.dot(v);
      EXPECT_NEAR(vol, lf::geometry::Volume(*cell->Geometry()), 1.0E-3)
          << "missmatch for cell " << cell << std::endl;
    }
  }

  // verification based on local edge vectors:
  {
    auto g = lf::mesh::utils::MeshFunctionConstant(1.0);
    ScalarLoadEdgeVectorProvider provider(fe_space, g);

    // Loop over edges and compute element vectors
    for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
      auto v = provider.Eval(*edge);

      Eigen::VectorXd one_vec = Eigen::VectorXd::Constant(v.size(), 1.0);

      const double length = one_vec.dot(v);
      EXPECT_NEAR(length, lf::geometry::Volume(*edge->Geometry()), 1.0E-3)
          << "missmatch for edge " << edge << std::endl;
    }
  }
}

TEST(lf_fe_linear, loc_comp_products_test) {
  std::cout << "Linear FE >>> Computation of local quantities" << std::endl;

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  std::shared_ptr<const UniformScalarFESpace<double>> fe_space =
      std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);
  locCompProductsTest(fe_space);
}

TEST(lf_fe_quadratic, loc_comp_products_test) {
  std::cout << "Quadratic FE >>> Computation of local quantities " << std::endl;

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  std::shared_ptr<const UniformScalarFESpace<double>> fe_space =
      std::make_shared<FeSpaceLagrangeO2<double>>(mesh_p);
  locCompProductsTest(fe_space);
}

TEST(lf_fe_cubic, loc_comp_products_test) {
  std::cout << "Cubic FE >>> Computation of local quantities " << std::endl;

  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  std::shared_ptr<const UniformScalarFESpace<double>> fe_space =
      std::make_shared<FeSpaceLagrangeO3<double>>(mesh_p);
  locCompProductsTest(fe_space);
}
}  // end namespace lf::uscalfe::test
