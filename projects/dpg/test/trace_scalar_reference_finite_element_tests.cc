/**
 * @file
 * @brief local  tests for the implementation
 * of the trace wrapper  around cubic
 * local shape functions
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <iostream>

#include <lf/uscalfe/uscalfe.h>

#include "../discontinuous_fe_constant.h"
#include "../lagr_fe_cubic.h"
#include "../lagr_fe_quadratic.h"
#include "../trace_scalar_reference_finite_element.h"
#include "lagr_test_utils.h"

namespace projects::dpg::test {

TEST(FeTrace, scal_fe_coeff_node) {
  std::cout << " >>> Trace FE: Test of consistency of nodal interpolation \n";

  TraceScalarReferenceFiniteElement<double> tfe(
      std::make_shared<const FeLagrangeO3Tria<double>>());
  EXPECT_TRUE(scalarFEEvalNodeTest(tfe));

  TraceScalarReferenceFiniteElement<double> qfe(
      std::make_shared<const FeLagrangeO3Quad<double>>());
  EXPECT_TRUE(scalarFEEvalNodeTest(qfe));

  TraceScalarReferenceFiniteElement<double> sfe(
      std::make_shared<const FeLagrangeO3Segment<double>>());
  EXPECT_TRUE(scalarFEEvalNodeTest(sfe));
}

TEST(FeTrace, scal_fe_val_node) {
  std::cout << ">>> Trace FE: Test of consistency of nodal interpolation \n";

  TraceScalarReferenceFiniteElement<double> tfe(
      std::make_shared<const FeLagrangeO3Tria<double>>());
  EXPECT_TRUE(scalarFEInterpTest(tfe));

  TraceScalarReferenceFiniteElement<double> qfe(
      std::make_shared<const FeLagrangeO3Quad<double>>());
  EXPECT_TRUE(scalarFEInterpTest(qfe));

  TraceScalarReferenceFiniteElement<double> sfe(
      std::make_shared<const FeLagrangeO3Segment<double>>());
  EXPECT_TRUE(scalarFEInterpTest(sfe));
}

TEST(FeTrace, shape_function_information) {
  std::cout << ">>> Trace FE: Test of number of reference shape functions \n";

  TraceScalarReferenceFiniteElement<double> tfe(
      std::make_shared<const FeLagrangeO3Tria<double>>());
  EXPECT_EQ(tfe.NumRefShapeFunctions(), 9);
  EXPECT_EQ(tfe.NumRefShapeFunctions(0, 0), 0);
  EXPECT_EQ(tfe.NumRefShapeFunctions(1, 0), 2);
  EXPECT_EQ(tfe.NumRefShapeFunctions(2, 0), 1);

  TraceScalarReferenceFiniteElement<double> qfe(
      std::make_shared<const FeLagrangeO3Quad<double>>());
  EXPECT_EQ(qfe.NumRefShapeFunctions(), 12);
  EXPECT_EQ(qfe.NumRefShapeFunctions(0, 0), 0);
  EXPECT_EQ(qfe.NumRefShapeFunctions(1, 0), 2);
  EXPECT_EQ(qfe.NumRefShapeFunctions(2, 0), 1);

  TraceScalarReferenceFiniteElement<double> sfe(
      std::make_shared<const FeLagrangeO3Segment<double>>());
  EXPECT_EQ(sfe.NumRefShapeFunctions(), 2);
  EXPECT_EQ(sfe.NumRefShapeFunctions(0, 0), 0);
  EXPECT_EQ(sfe.NumRefShapeFunctions(1, 0), 1);
}

TEST(FeTrace, shape_function_cardinality) {
  std::cout << ">>> Trace Fe: Test of cardinal basis property of "
               "reference shape functions \n";

  TraceScalarReferenceFiniteElement<double> tfe(
      std::make_shared<const FeLagrangeO3Tria<double>>());
  EXPECT_TRUE(tfe.EvalReferenceShapeFunctions(tfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(9, 9)));

  TraceScalarReferenceFiniteElement<double> qfe(
      std::make_shared<const FeLagrangeO3Quad<double>>());
  EXPECT_TRUE(qfe.EvalReferenceShapeFunctions(qfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(12, 12)));

  TraceScalarReferenceFiniteElement<double> sfe(
      std::make_shared<const FeLagrangeO3Segment<double>>());
  EXPECT_TRUE(sfe.EvalReferenceShapeFunctions(sfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(2, 2)));
}

// note: gradients no longer sum up to 0.

}  // namespace projects::dpg::test
