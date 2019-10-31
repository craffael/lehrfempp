/**
 * @file
 * @brief local  tests for the implementation
 * of the discontinuous wrapper
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <iostream>

#include <lf/uscalfe/uscalfe.h>

#include "../discontinuous_scalar_reference_finite_element.h"
#include "../lagr_fe_cubic.h"
#include "lagr_test_utils.h"

namespace projects::dpg::test {

TEST(FeDisc, scal_fe_coeff_node) {
  std::cout
      << " >>> Discontinuous FE: Test of consistency of nodal interpolation \n";

  DiscontinuousScalarReferenceFiniteElement<double> tfe(
      std::make_shared<const FeLagrangeO3Tria<double>>());
  EXPECT_TRUE(scalarFEEvalNodeTest(tfe));

  DiscontinuousScalarReferenceFiniteElement<double> qfe(
      std::make_shared<const FeLagrangeO3Quad<double>>());
  EXPECT_TRUE(scalarFEEvalNodeTest(qfe));

  DiscontinuousScalarReferenceFiniteElement<double> sfe(
      std::make_shared<const FeLagrangeO3Segment<double>>());
  EXPECT_TRUE(scalarFEEvalNodeTest(sfe));
}

TEST(FeDisc, scal_fe_val_node) {
  std::cout
      << ">>> Discontinuous FE: Test of consistency of nodal interpolation \n";

  DiscontinuousScalarReferenceFiniteElement<double> tfe(
      std::make_shared<const FeLagrangeO3Tria<double>>());
  EXPECT_TRUE(scalarFEInterpTest(tfe));

  DiscontinuousScalarReferenceFiniteElement<double> qfe(
      std::make_shared<const FeLagrangeO3Quad<double>>());
  EXPECT_TRUE(scalarFEInterpTest(qfe));

  DiscontinuousScalarReferenceFiniteElement<double> sfe(
      std::make_shared<const FeLagrangeO3Segment<double>>());
  EXPECT_TRUE(scalarFEInterpTest(sfe));
}

TEST(FeDisc, shape_function_information) {
  std::cout
      << ">>> Discontinuous FE: Test of number of reference shape functions \n";

  DiscontinuousScalarReferenceFiniteElement<double> tfe(
      std::make_shared<const FeLagrangeO3Tria<double>>());
  EXPECT_EQ(tfe.NumRefShapeFunctions(), 10);
  EXPECT_EQ(tfe.NumRefShapeFunctions(0, 0), 10);
  EXPECT_EQ(tfe.NumRefShapeFunctions(1, 0), 0);
  EXPECT_EQ(tfe.NumRefShapeFunctions(2, 0), 0);

  DiscontinuousScalarReferenceFiniteElement<double> qfe(
      std::make_shared<const FeLagrangeO3Quad<double>>());
  EXPECT_EQ(qfe.NumRefShapeFunctions(), 16);
  EXPECT_EQ(qfe.NumRefShapeFunctions(0, 0), 16);
  EXPECT_EQ(qfe.NumRefShapeFunctions(1, 0), 0);
  EXPECT_EQ(qfe.NumRefShapeFunctions(2, 0), 0);

  DiscontinuousScalarReferenceFiniteElement<double> sfe(
      std::make_shared<const FeLagrangeO3Segment<double>>());
  EXPECT_EQ(sfe.NumRefShapeFunctions(), 4);
  EXPECT_EQ(sfe.NumRefShapeFunctions(0, 0), 4);
  EXPECT_EQ(sfe.NumRefShapeFunctions(1, 0), 0);
}

TEST(FeDisc, shape_function_cardinality) {
  std::cout << ">>> Discontinous Fe: Test of cardinal basis property of "
               "reference shape functions \n";

  DiscontinuousScalarReferenceFiniteElement<double> tfe(
      std::make_shared<const FeLagrangeO3Tria<double>>());
  EXPECT_TRUE(tfe.EvalReferenceShapeFunctions(tfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(10, 10)));

  DiscontinuousScalarReferenceFiniteElement<double> qfe(
      std::make_shared<const FeLagrangeO3Quad<double>>());
  EXPECT_TRUE(qfe.EvalReferenceShapeFunctions(qfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(16, 16)));

  DiscontinuousScalarReferenceFiniteElement<double> sfe(
      std::make_shared<const FeLagrangeO3Segment<double>>());
  EXPECT_TRUE(sfe.EvalReferenceShapeFunctions(sfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(4, 4)));
}

TEST(FeDisc, shape_function_gradient_sum) {
  std::cout << " >>> Discontinuous FE: Test of gadients summing up to 0 \n";

  DiscontinuousScalarReferenceFiniteElement<double> tfe(
      std::make_shared<const FeLagrangeO3Tria<double>>());
  EXPECT_TRUE(scalarFEEvalGRadTest(tfe));

  DiscontinuousScalarReferenceFiniteElement<double> qfe(
      std::make_shared<const FeLagrangeO3Quad<double>>());
  EXPECT_TRUE(scalarFEEvalGRadTest(qfe));

  DiscontinuousScalarReferenceFiniteElement<double> sfe(
      std::make_shared<const FeLagrangeO3Segment<double>>());
  EXPECT_TRUE(scalarFEEvalGRadTest(sfe));
}

}  // namespace projects::dpg::test
