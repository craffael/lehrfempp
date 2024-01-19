/* **************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Unit tests for local element matrix builders
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <iostream>

namespace lf::uscalfe {

// Make sure the ElementMatrix/VectorProviders in the lf::fe namespace only
// depend on the archetype functionality by explicitly instantiating them:
template class ReactionDiffusionElementMatrixProvider<
    double, mesh::utils::MeshFunctionAT<double>,
    mesh::utils::MeshFunctionAT<double>>;
template class ReactionDiffusionElementMatrixProvider<
    double, mesh::utils::MeshFunctionAT<Eigen::Vector2d>,
    mesh::utils::MeshFunctionAT<Eigen::Vector2d>>;

template class MassEdgeMatrixProvider<double, mesh::utils::MeshFunctionAT<double>, base::PredicateTrue>;

template class ScalarLoadElementVectorProvider<double, mesh::utils::MeshFunctionAT<double>>;

template class ScalarLoadEdgeVectorProvider<double, mesh::utils::MeshFunctionAT<double>, base::PredicateTrue>;

}  // namespace lf::uscalfe

namespace lf::uscalfe::test {

TEST(lf_uscalfe, mass_mat_test) {
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Quadrature rules
  quad_rule_collection_t quad_rules{
      {lf::base::RefEl::kTria(), lf::quad::make_TriaQR_P3O3()},
      {lf::base::RefEl::kQuad(), lf::quad::make_QuadQR_P1O2()}};

  // Coefficients
  auto alpha = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> double { return 1.0; });
  auto gamma = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d) -> double { return 1.0; });

  ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                         decltype(gamma)>
      elmat_builder(fe_space, alpha, gamma, quad_rules);

  // Traverse the cells of the mesh and compute element matrices
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    const typename ReactionDiffusionElementMatrixProvider<
        double, decltype(alpha), decltype(gamma)>::ElemMat M{
        elmat_builder.Eval(*cell)};
    // Vectors with all 1 components
    Eigen::VectorXd one_vec_c = Eigen::VectorXd::Constant(M.cols(), 1.0);
    Eigen::VectorXd one_vec_r = Eigen::VectorXd::Constant(M.rows(), 1.0);
    // Multiply element matrix with these vectors from left and right.
    const double volint = one_vec_r.dot(M * one_vec_c);
    // This operations supresses the contribution of the second-order term
    // the value amounts to the integral of gamma over the cell.
    // In this case this is just the volume.
    EXPECT_NEAR(volint, lf::geometry::Volume(*(cell->Geometry())), 1.0E-10)
        << " mismatch for cell " << *cell;
  }
}

TEST(lf_uscalfe, cross_val) {
  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // In this test we compare the results of the implementation with
  // user-supplied quadrature rules and with the internal allocation of
  // quadrature rules. We choose the quadrature rules in the same way as the
  // internal selection so that the results should be exactly the same

  // Quadrature rules of degree of exactness = 2
  quad_rule_collection_t quad_rules{
      {lf::base::RefEl::kTria(),
       lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2)},
      {lf::base::RefEl::kQuad(),
       lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 2)}};

  // Coefficients, reasonably complicated
  auto alpha = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return x[0] * x[0] + x[1] * x[1]; });

  auto gamma = mesh::utils::MeshFunctionGlobal([](Eigen::Vector2d x) -> double {
    return 1.0 / (1 + x[0] * x[0] + x[1] * x[1]);
  });

  auto f = mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return x[0] * x[0] - x[1] * x[1]; });

  ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                         decltype(gamma)>
      elmat_builder_qr(fe_space, alpha, gamma, quad_rules);

  ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                         decltype(gamma)>
      elmat_builder_noqr(fe_space, alpha, gamma);

  ScalarLoadElementVectorProvider<double, decltype(f)> elvec_builder_qr(
      fe_space, f, quad_rules);

  ScalarLoadElementVectorProvider<double, decltype(f)> elvec_builder_noqr(
      fe_space, f);

  // Traverse the cells of the mesh and compute element matrices
  for (const lf::mesh::Entity *cell : mesh_p->Entities(0)) {
    const typename ReactionDiffusionElementMatrixProvider<
        double, decltype(alpha), decltype(gamma)>::ElemMat M_qr{
        elmat_builder_qr.Eval(*cell)};

    const typename ReactionDiffusionElementMatrixProvider<
        double, decltype(alpha), decltype(gamma)>::ElemMat M_no{
        elmat_builder_noqr.Eval(*cell)};

    EXPECT_NEAR((M_qr - M_no).norm(), 0.0, 1.0E-10)
        << "M_qr = " << M_qr << " M_no = " << M_no;

    const typename ScalarLoadElementVectorProvider<double, decltype(f)>::ElemVec
        phi_qr{elvec_builder_qr.Eval(*cell)};

    const typename ScalarLoadElementVectorProvider<double, decltype(f)>::ElemVec
        phi_no{elvec_builder_noqr.Eval(*cell)};

    EXPECT_NEAR((phi_qr - phi_no).norm(), 0.0, 1.0E-10)
        << "phi_qr = " << phi_qr << " phi_no = " << phi_no;
  }
}

TEST(lf_uscalfe, missingQuadRule) {
  // We test that an exception is thrown when a quadrature rule is missing
  // for a hybrid mesh and we call Eval() for a quad cell

  // Building the test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  // Set up finite elements
  auto fe_space = std::make_shared<FeSpaceLagrangeO1<double>>(mesh_p);

  // Quadrature rules of degree of exactness = 2
  quad_rule_collection_t quad_rules{
      {lf::base::RefEl::kTria(),
       lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2)}};

  // Coefficients, reasonably complicated
  auto coeff = mesh::utils::MeshFunctionConstant(2.0);

  ReactionDiffusionElementMatrixProvider rdemp(fe_space, coeff, coeff,
                                               quad_rules);

  for (auto e : mesh_p->Entities(0)) {
    if (e->RefEl() == base::RefEl::kTria()) {
      // Check NoThrow:
      EXPECT_NO_THROW(rdemp.Eval(*e));
    } else {
      EXPECT_THROW(rdemp.Eval(*e), base::LfException);
    }
  }
}

}  // namespace lf::uscalfe::test
