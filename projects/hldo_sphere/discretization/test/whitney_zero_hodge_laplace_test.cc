#include <gtest/gtest.h>
#include <whitney_zero_hodge_laplacian.h>

#include <array>
#include <cmath>

TEST(projects_hldo_sphere_discretization,
     whitney_zero_hodge_laplace_basic_test) {
  // Build LSE
  projects::hldo_sphere::discretization::WhitneyZeroHodgeLaplace lse_builder;
  lse_builder.Compute();

  Eigen::SparseMatrix<double> Ae = lse_builder.GetGalerkinMatrix().makeSparse();

  Eigen::MatrixXd Ae_anal(6, 6);
  // clang-format off
  Ae_anal <<     8, -2, -2, -2, -2,  0,
                -2,  8, -2,  0, -2, -2,  
                -2, -2,  8, -2,  0, -2,  
                -2,  0, -2,  8, -2, -2,  
                -2, -2,  0, -2,  8, -2,
                 0, -2, -2, -2, -2,  8;
  // clang-format on
  Ae_anal *= 1. / 2. / std::sqrt(3);
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_DOUBLE_EQ(Ae.coeff(i, j), Ae_anal(i, j))
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }

  Eigen::VectorXd Vec = lse_builder.GetLoadVector();
  Eigen::VectorXd Vec_anal(6);
  // clang-format off
  Vec_anal <<     0, 0, 0, 0, 0,  0;
  // clang-format on

  for (int j = 0; j < Ae_anal.cols(); ++j) {
    EXPECT_DOUBLE_EQ(Vec(j), Vec_anal(j)) << "mismatch in entry (" << j << ")";
  }
}
