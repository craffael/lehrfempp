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
#include "lf/fe/lagr_fe.h"

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
  std::cout << "Tria: Evaluation nodes\n" << tlfe.EvaluationNodes() << std::endl;
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
  std::cout << "Quad: Evaluation nodes\n" << qlfe.EvaluationNodes() << std::endl;
  }
}

}  // end namespace lf::fe::test
