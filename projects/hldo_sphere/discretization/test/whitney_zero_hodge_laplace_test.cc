#include <gtest/gtest.h>
#include <whitney_one_hodge_laplacian.h>
#include <whitney_two_hodge_laplacian.h>
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

TEST(projects_hldo_sphere_discretization,
     whitney_one_hodge_laplace_basic_test) {
  // Build LSE
  projects::hldo_sphere::discretization::WhitneyOneHodgeLaplace lse_builder;
  lse_builder.Compute();

  Eigen::SparseMatrix<double> Ae = lse_builder.GetGalerkinMatrix().makeSparse();

  Eigen::MatrixXd Ae_anal(18, 18);
  // clang-format off
  Ae_anal <<     4,  2,  0,  2,  2, -2,  0,  0,  0,  0,  0,  0,    1, -1,  0,  0,  0,  0,
                 2,  4, -2,  0,  2,  0,  0, -2,  0,  0,  0,  0,   -1,  0,  1,  0,  0,  0,
                 0, -2,  4, -2,  0,  0,  0,  2,  0, -2,  0,  0,   -1,  0,  0,  1,  0,  0,
                 2,  0, -2,  4,  0, -2,  0,  0,  0,  2,  0,  0,   -1,  0,  0,  0,  1,  0,
                 2,  2,  0,  0,  4,  0, -2,  0, -2,  0,  0,  0,    0,  1, -1,  0,  0,  0,
                -2,  0,  0, -2,  0,  4,  2,  0,  0,  0,  0,  2,    0, -1,  0,  0,  1,  0,
                 0,  0,  0,  0, -2,  2,  4,  0,  2,  0,  0,  2,    0,  1,  0,  0,  0, -1, 
                 0, -2,  2,  0,  0,  0,  0,  4,  2,  0, -2,  0,    0,  0,  1, -1,  0,  0,
                 0,  0,  0,  0, -2,  0,  2,  2,  4,  0, -2,  0,    0,  0, -1,  0,  0,  1,
                 0,  0, -2,  2,  0,  0,  0,  0,  0,  4,  2, -2,    0,  0,  0,  1, -1,  0,
                 0,  0,  0,  0,  0,  0,  0, -2, -2,  2,  4, -2,    0,  0,  0, -1,  0,  1,
                 0,  0,  0,  0,  0,  2,  2,  0,  0, -2, -2,  4,    0,  0,  0,  0, -1,  1,

                -1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,    0,  0,  0,  0,  0,  0,
                 1,  0,  0,  0, -1,  1, -1,  0,  0,  0,  0,  0,    0,  0,  0,  0,  0,  0,
                 0, -1,  0,  0,  1,  0,  0, -1,  1,  0,  0,  0,    0,  0,  0,  0,  0,  0,
                 0,  0, -1,  0,  0,  0,  0,  1,  0, -1,  1,  0,    0,  0,  0,  0,  0,  0,
                 0,  0,  0, -1,  0, -1,  0,  0,  0,  1,  0,  1,    0,  0,  0,  0,  0,  0,
                 0,  0,  0,  0,  0,  0,  1,  0, -1,  0, -1, -1,    0,  0,  0,  0,  0,  0;
  // clang-format on
  Ae_anal *= 1. / std::sqrt(3);

  Eigen::MatrixXd Ae_anal_mass(6, 6);
  // clang-format off
  Ae_anal_mass <<    4, 1, 1, 1, 1, 0,
                     1, 4, 1, 0, 1, 1,
                     1, 1, 4, 1, 0, 1,
                     1, 0, 1, 4, 1, 1,
                     1, 1, 0, 1, 4, 1, 
                     0, 1, 1, 1, 1, 4;
  // clang-format on
  Ae_anal_mass *= 1. / std::sqrt(3) / 4.;
  Ae_anal.bottomRightCorner(6, 6) = Ae_anal_mass;
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae.coeff(i, j), Ae_anal(i, j), 1e-12)
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }

  Eigen::VectorXd Vec = lse_builder.GetLoadVector();
  Eigen::VectorXd Vec_anal(18);
  Vec_anal.setZero();
  // clang-format off
  // clang-format on

  for (int j = 0; j < Ae_anal.cols(); ++j) {
    EXPECT_NEAR(Vec(j), Vec_anal(j), 1e-12)
        << "mismatch in entry (" << j << ")";
  }
}

TEST(projects_hldo_sphere_discretization,
     whitney_two_hodge_laplace_basic_test) {
  // Build LSE
  projects::hldo_sphere::discretization::WhitneyTwoHodgeLaplace lse_builder;
  lse_builder.Compute();

  Eigen::SparseMatrix<double> Ae = lse_builder.GetGalerkinMatrix().makeSparse();

  Eigen::MatrixXd Ae_anal(20, 20);
  // clang-format off
  Ae_anal <<  
           10, -1,  0, -1, -1,  1,  0,  0,  0,  0,  0,  0,    -1,  0,  0,  0,  0,  0,  1,  0,
           -1, 10,  1,  0, -1,  0,  0,  1,  0,  0,  0,  0,    -1,  0,  1,  0,  0,  0,  0,  0,
            0,  1, 10,  1,  0,  0,  0, -1,  0,  1,  0,  0,     0,  0, -1,  0,  1,  0,  0,  0,
           -1,  0,  1, 10,  0,  1,  0,  0,  0, -1,  0,  0,     0,  0,  0,  0, -1,  0,  1,  0,
           -1, -1,  0,  0, 10,  0,  1,  0,  1,  0,  0,  0,    -1,  1,  0,  0,  0,  0,  0,  0, 
            1,  0,  0,  1,  0, 10, -1,  0,  0,  0,  0, -1,     0,  0,  0,  0,  0,  0, -1,  1,
            0,  0,  0,  0,  1, -1, 10,  0, -1,  0,  0, -1,     0, -1,  0,  0,  0,  0,  0,  1,
            0,  1, -1,  0,  0,  0,  0, 10, -1,  0,  1,  0,     0,  0, -1,  1,  0,  0,  0,  0,
            0,  0,  0,  0,  1,  0, -1, -1, 10,  0,  1,  0,     0, -1,  0,  1,  0,  0,  0,  0,
            0,  0,  1, -1,  0,  0,  0,  0,  0, 10, -1,  1,     0,  0,  0,  0, -1,  1,  0,  0, 
            0,  0,  0,  0,  0,  0,  0,  1,  1, -1, 10,  1,     0,  0,  0, -1,  0,  1,  0,  0,
            0,  0,  0,  0,  0, -1, -1,  0,  0,  1,  1, 10,     0,  0,  0,  0,  0, -1,  0,  1,

           -1, -1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,     0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  1,  0, -1,  0, -1,  0,  0,  0,     0,  0,  0,  0,  0,  0,  0,  0,
            0,  1, -1,  0,  0,  0,  0, -1,  0,  0,  0,  0,     0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  1,  1,  0, -1,  0,     0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  1, -1,  0,  0,  0,  0,  0, -1,  0,  0,     0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1, -1,     0,  0,  0,  0,  0,  0,  0,  0,
            1,  0,  0,  1,  0, -1,  0,  0,  0,  0,  0,  0,     0,  0,  0,  0,  0,  0,  0,  0,
            0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  1,     0,  0,  0,  0,  0,  0,  0,  0;
  // clang-format on
  Ae_anal.topLeftCorner(12, 12) *= 1. / std::sqrt(3) / 12.;
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae.coeff(i, j), Ae_anal(i, j), 1e-12)
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }

  Eigen::VectorXd Vec = lse_builder.GetLoadVector();
  Eigen::VectorXd Vec_anal(20);
  Vec_anal.setZero();
  // clang-format off
  // clang-format on

  for (int j = 0; j < Ae_anal.cols(); ++j) {
    EXPECT_NEAR(Vec(j), Vec_anal(j), 1e-12)
        << "mismatch in entry (" << j << ")";
  }
}
