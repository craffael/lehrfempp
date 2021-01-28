/**
 * @file
 * @brief adapted test utilities from the lehrfem framework
 * to test lagrangian shape functions both locally and globally.
 * @author
 * @date June 2019
 * @copyright MIT License
 */

#ifndef PROJECTS_DPG_TEST_LAGR_TEST_UTILS
#define PROJECTS_DPG_TEST_LAGR_TEST_UTILS

#include <lf/fe/fe.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/uscalfe/test/fe_testutils.h>
#include <lf/uscalfe/uscalfe.h>

#include "../dpg.h"

namespace projects::dpg::test {

// adapted test from lf/uscalfe/lagr_fe_test.cc
// since now the evaluation nodes are somtimes placed at 1/3 or 2/3
// the stritct double equality  in the original test could no longer be
// satisfied and the constraint was weakened.
template <typename SCALAR>
bool scalarFEEvalNodeTest(
    const lf::fe::ScalarReferenceFiniteElement<SCALAR> &fe_desc) {
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

// adapted test from lf/uscalfe/lagr_fe_test.cc
// since now the evaluation nodes are somtimes placed at 1/3 or 2/3
// the double equality in the original test could no longer be satisfied
// and the constraint was weakened.
template <typename SCALAR>
bool scalarFEInterpTest(
    const lf::fe::ScalarReferenceFiniteElement<SCALAR> &fe_desc) {
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

// checks if the gradients sum up to 0
template <typename SCALAR>
bool scalarFEEvalGRadTest(
    const lf::fe::ScalarReferenceFiniteElement<SCALAR> &fe_desc) {
  lf::base::RefEl refel = fe_desc.RefEl();
  unsigned dim = refel.Dimension();

  const Eigen::MatrixXd evl_nodes{fe_desc.EvaluationNodes()};
  const size_type N_evln = fe_desc.NumEvaluationNodes();
  EXPECT_EQ(evl_nodes.cols(), N_evln);

  const size_type N_grads = fe_desc.NumRefShapeFunctions();
  auto grads_at_evln = fe_desc.GradientsReferenceShapeFunctions(evl_nodes);
  EXPECT_EQ(grads_at_evln.rows(), N_grads);
  EXPECT_EQ(grads_at_evln.cols(), dim * N_evln);

  // check, that the sum of gradients is equal to zero, in all the evaluaiton
  // nodes.
  for (int i = 0; i < dim * N_evln; ++i) {
    double sum = grads_at_evln.col(i).sum();
    EXPECT_NEAR(sum, 0.0, 1e-13);
  }
  return true;
}

using fe_ptr =
    std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<double>>;
// adaption of the lf::uscalfe::test in full_gal_tests.cc to
// support fulll galerkin test , with second and third order shape functions.
template <typename FFUNC, typename DIFF_COEFF, typename REAC_COEFF>
std::vector<double> FEEnergyTest(int reflevels, FFUNC v, DIFF_COEFF alpha,
                                 REAC_COEFF gamma, fe_ptr rfs_tria,
                                 fe_ptr rfs_quad) {
  // generate hierarchy of meshes of unit square
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(
          lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3),
          reflevels);
  std::vector<double> energies{
      lf::uscalfe::test::EnergiesOfInterpolants<double>(
          *multi_mesh_p, v, alpha, gamma, rfs_tria, rfs_quad)};

  size_type L = energies.size();
  for (int l = 0; l < L; ++l) {
    std::cout << "Level" << l << ": energy = " << energies[l] << std::endl;
  }
  return energies;
}

// adaption of the lf::uscalfe::test in full_gal_tests.cc
// to support a full galerkin test, with second and third order shape functions.
template <typename FFUNC, typename IMP_COEFF>
std::vector<double> FEInterfaceEnergyTest(int reflevels, FFUNC v, IMP_COEFF eta,
                                          fe_ptr rfs_tria, fe_ptr rfs_quad,
                                          fe_ptr rfs_segment) {
  // Generate hierarchy of meshes of the unit square
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(
          lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3),
          reflevels);
  // Select edges on top and right side of unit square
  auto pos_pred = [](Eigen::Vector2d x) -> bool {
    return ((x[0] + x[1]) > 1.0);
  };

  lf::refinement::EntityCenterPositionSelector edge_sel{pos_pred};

  // Compute energies of the  interpolants of v on the different
  // refinement levels
  std::vector<double> energies{
      lf::uscalfe::test::BoundaryEnergiesOfInterpolants<double>(
          *multi_mesh_p, v, eta, rfs_tria, rfs_quad, rfs_segment, edge_sel)};
  // Output of energies
  size_type L = energies.size();
  for (int l = 0; l < L; ++l) {
    std::cout << "Level" << l << ": energy = " << energies[l] << std::endl;
  }
  return energies;
}

}  // namespace projects::dpg::test

#endif  // PROJECTS_DPG_TEST_LAGR_TEST_UTILS
