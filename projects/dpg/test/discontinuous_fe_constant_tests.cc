
/**
 * @file
 * @brief local tests for the implementation of constant
 * local shape functions
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */
#include <gtest/gtest.h>
#include <iostream>

#include <lf/uscalfe/uscalfe.h>

#include "../discontinuous_fe_constant.h"
#include "lagr_test_utils.h"

namespace projects::dpg::test {

TEST(FeLagrO0, scal_fe_coeff_node) {
  std::cout << ">>> O0 FE: Test of consistency of nodal interpolation \n";

  FeDiscontinuousO0Tria<double> tfe{};
  EXPECT_TRUE(scalarFEEvalNodeTest(tfe));

  FeDiscontinuousO0Quad<double> qfe{};
  EXPECT_TRUE(scalarFEEvalNodeTest(qfe));

  FeDiscontinuousO0Segment<double> sfe{};
  EXPECT_TRUE(scalarFEEvalNodeTest(sfe));
}

TEST(FeLagrO0, scal_fe_val_node) {
  std::cout << ">>> O0 FE: Test of consistency of nodal interpolation \n";

  FeDiscontinuousO0Tria<double> tfe{};
  EXPECT_TRUE(scalarFEInterpTest(tfe));

  FeDiscontinuousO0Quad<double> qfe{};
  EXPECT_TRUE(scalarFEInterpTest(qfe));

  FeDiscontinuousO0Segment<double> sfe{};
  EXPECT_TRUE(scalarFEInterpTest(sfe));
}

TEST(FeLagrO0, shape_function_information) {
  std::cout << ">>> O0 Fe: Test of number of reference shape functions \n";

  FeDiscontinuousO0Tria<double> tfe{};
  EXPECT_EQ(tfe.NumRefShapeFunctions(), 1);
  EXPECT_EQ(tfe.NumRefShapeFunctions(0, 0), 1);
  EXPECT_EQ(tfe.NumRefShapeFunctions(1, 0), 0);
  EXPECT_EQ(tfe.NumRefShapeFunctions(2, 0), 0);

  FeDiscontinuousO0Quad<double> qfe{};
  // test number of basis functions:
  EXPECT_EQ(qfe.NumRefShapeFunctions(), 1);
  EXPECT_EQ(qfe.NumRefShapeFunctions(0, 0), 1);
  EXPECT_EQ(qfe.NumRefShapeFunctions(1, 0), 0);
  EXPECT_EQ(qfe.NumRefShapeFunctions(2, 0), 0);

  FeDiscontinuousO0Segment<double> sfe{};
  EXPECT_EQ(sfe.NumRefShapeFunctions(), 1);
  EXPECT_EQ(sfe.NumRefShapeFunctions(0, 0), 1);
  EXPECT_EQ(sfe.NumRefShapeFunctions(1, 0), 0);
}

TEST(FeLagrO0, shape_function_cardinality) {
  std::cout << ">>> O0 Fe: Test of cardinal basis property of reference shape "
               "functions \n";

  FeDiscontinuousO0Tria<double> tfe{};
  EXPECT_TRUE(tfe.EvalReferenceShapeFunctions(tfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(1, 1)));

  FeDiscontinuousO0Quad<double> qfe{};
  EXPECT_TRUE(qfe.EvalReferenceShapeFunctions(qfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(1, 1)));

  FeDiscontinuousO0Segment<double> sfe{};
  EXPECT_TRUE(sfe.EvalReferenceShapeFunctions(sfe.EvaluationNodes())
                  .isApprox(Eigen::MatrixXd::Identity(1, 1)));
}

TEST(FeLagrO0, shape_function_gradient_sum) {
  std::cout << ">>> O0 Fe: Test of gadients summing up to 0 \n";

  FeDiscontinuousO0Tria<double> tfe{};
  EXPECT_TRUE(scalarFEEvalGRadTest(tfe));

  FeDiscontinuousO0Quad<double> qfe{};
  EXPECT_TRUE(scalarFEEvalGRadTest(qfe));

  FeDiscontinuousO0Segment<double> sfe{};
  EXPECT_TRUE(scalarFEEvalGRadTest(sfe));
}

}  // namespace projects::dpg::test
