#include <gtest/gtest.h>
#include <hodge_laplacians_source_problems.h>

#include <array>
#include <cmath>

TEST(projects_hldo_sphere_discretization, hodge_laplacian_zero_basic_test) {
  // Build LSE
  projects::hldo_sphere::discretization::HodgeLaplaciansSourceProblems
      lse_builder;
  lse_builder.Compute();
  lse_builder.Solve();

  Eigen::SparseMatrix<double> Ae =
      lse_builder.GetGalerkinMatrix(0).makeSparse();

  Eigen::MatrixXd Ae_anal(6, 6);
  // clang-format off
  Ae_anal <<     8, -2, -2, -2, -2,  0,
                -2,  8, -2,  0, -2, -2,  
                -2, -2,  8, -2,  0, -2,  
                -2,  0, -2,  8, -2, -2,  
                -2, -2,  0, -2,  8, -2,
                 0, -2, -2, -2, -2,  8;
  // clang-format on
  Ae_anal *= -1. / 2. / std::sqrt(3);
  Eigen::MatrixXd Ae_mass(6, 6);
  // clang-format off
  Ae_mass <<     8,  2,  2,  2,  2,  0,
                 2,  8,  2,  0,  2,  2,  
                 2,  2,  8,  2,  0,  2,  
                 2,  0,  2,  8,  2,  2,  
                 2,  2,  0,  2,  8,  2,
                 0,  2,  2,  2,  2,  8;
  // clang-format on
  Ae_mass *= 1. / 8. / sqrt(3);
  Ae_anal += Ae_mass;

  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae.coeff(i, j), Ae_anal(i, j), 1e-13)
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }

  Eigen::VectorXd Vec = lse_builder.GetLoadVector(0);
  Eigen::VectorXd Vec_anal(6);
  // clang-format off
  Vec_anal << 0, 0, 0, 0, 0, 0;
  // clang-format on

  ASSERT_EQ(Vec.rows(), Vec_anal.rows());
  for (int j = 0; j < Ae_anal.rows(); ++j) {
    EXPECT_NEAR(Vec(j), Vec_anal(j), 1e-13)
        << "mismatch in entry (" << j << ")";
  }

  // result verctor mu
  Eigen::VectorXcd mu_anal_zero = Eigen::VectorXcd::Zero(6);
  Eigen::VectorXcd mu_zero = lse_builder.GetMuZero();
  ASSERT_EQ(mu_zero.rows(), mu_anal_zero.rows());
  for (int j = 0; j < mu_anal_zero.cols(); ++j) {
    EXPECT_NEAR(mu_zero(j).imag(), mu_anal_zero(j).imag(), 1e-12)
        << "mismatch in imag entry (" << j << ")";
    EXPECT_NEAR(mu_zero(j).real(), mu_anal_zero(j).real(), 1e-12)
        << "mismatch in real entry (" << j << ")";
  }
}

TEST(projects_hldo_sphere_discretization,
     hodge_laplacians_one_form_basic_test) {
  // Build LSE
  projects::hldo_sphere::discretization::HodgeLaplaciansSourceProblems
      lse_builder;
  lse_builder.Compute();
  lse_builder.Solve();

  Eigen::SparseMatrix<double> Ae =
      lse_builder.GetGalerkinMatrix(1).makeSparse();

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
  Ae_anal *= -1. / std::sqrt(3);

  Eigen::MatrixXd Ae_anal_mass_zero(6, 6);
  // clang-format off
  Ae_anal_mass_zero <<    4, 1, 1, 1, 1, 0,
                          1, 4, 1, 0, 1, 1,
                          1, 1, 4, 1, 0, 1,
                          1, 0, 1, 4, 1, 1,
                          1, 1, 0, 1, 4, 1, 
                          0, 1, 1, 1, 1, 4;
  // clang-format on
  Ae_anal_mass_zero *= -1. / std::sqrt(3) / 4.;
  Ae_anal.bottomRightCorner(6, 6) = Ae_anal_mass_zero;

  Eigen::MatrixXd Ae_anal_mass_one(12, 12);
  // clang-format off
  Ae_anal_mass_one <<   10, -1,  0, -1, -1,  1,  0,  0,  0,  0,  0,  0, 
                        -1, 10,  1,  0, -1,  0,  0,  1,  0,  0,  0,  0, 
                         0,  1, 10,  1,  0,  0,  0, -1,  0,  1,  0,  0, 
                        -1,  0,  1, 10,  0,  1,  0,  0,  0, -1,  0,  0, 
                        -1, -1,  0,  0, 10,  0,  1,  0,  1,  0,  0,  0, 
                         1,  0,  0,  1,  0, 10, -1,  0,  0,  0,  0, -1, 
                         0,  0,  0,  0,  1, -1, 10,  0, -1,  0,  0, -1, 
                         0,  1, -1,  0,  0,  0,  0, 10, -1,  0,  1,  0, 
                         0,  0,  0,  0,  1,  0, -1, -1, 10,  0,  1,  0, 
                         0,  0,  1, -1,  0,  0,  0,  0,  0, 10, -1,  1, 
                         0,  0,  0,  0,  0,  0,  0,  1,  1, -1, 10,  1, 
                         0,  0,  0,  0,  0, -1, -1,  0,  0,  1,  1, 10;
  // clang-format on
  Ae_anal_mass_one *= 1. / std::sqrt(3) / 12.;
  Ae_anal.topLeftCorner(12, 12) += Ae_anal_mass_one;

  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae.coeff(i, j), Ae_anal(i, j), 1e-12)
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }

  Eigen::VectorXd Vec = lse_builder.GetLoadVector(1);
  Eigen::VectorXd Vec_anal(18);
  Vec_anal.setZero();

  ASSERT_EQ(Vec.rows(), Vec_anal.rows());
  for (int j = 0; j < Ae_anal.rows(); ++j) {
    EXPECT_NEAR(Vec(j), Vec_anal(j), 1e-12)
        << "mismatch in entry (" << j << ")";
  }

  // result verctor mu
  Eigen::VectorXcd mu_anal_one = Eigen::VectorXcd::Zero(12);
  Eigen::VectorXcd mu_one = std::get<0>(lse_builder.GetMuOne());
  ASSERT_EQ(mu_one.rows(), mu_anal_one.rows());
  for (int j = 0; j < mu_anal_one.rows(); ++j) {
    EXPECT_NEAR(mu_one(j).imag(), mu_anal_one(j).imag(), 1e-12)
        << "mismatch in imag entry (" << j << ")";
    EXPECT_NEAR(mu_one(j).real(), mu_anal_one(j).real(), 1e-12)
        << "mismatch in real entry (" << j << ")";
  }

  Eigen::VectorXcd rho_anal_one = Eigen::VectorXcd::Zero(6);
  Eigen::VectorXcd rho_one = std::get<1>(lse_builder.GetMuOne());
  ASSERT_EQ(rho_one.rows(), rho_anal_one.rows());
  for (int j = 0; j < rho_anal_one.rows(); ++j) {
    EXPECT_NEAR(rho_one(j).imag(), rho_anal_one(j).imag(), 1e-12)
        << "mismatch in imag entry (" << j << ")";
    EXPECT_NEAR(rho_one(j).real(), rho_anal_one(j).real(), 1e-12)
        << "mismatch in real entry (" << j << ")";
  }
}

TEST(projects_hldo_sphere_discretization, hodge_laplace_two_form_basic_test) {
  // Build LSE
  projects::hldo_sphere::discretization::HodgeLaplaciansSourceProblems
      lse_builder;
  lse_builder.Compute();
  lse_builder.Solve();

  Eigen::SparseMatrix<double> Ae =
      lse_builder.GetGalerkinMatrix(2).makeSparse();

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
  Ae_anal *= -1.;

  Eigen::MatrixXd Ae_anal_mass_two = Eigen::MatrixXd::Identity(8, 8);
  Ae_anal_mass_two *= sqrt(3) / 2.;

  Ae_anal.bottomRightCorner(8, 8) = Ae_anal_mass_two;

  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae.coeff(i, j), Ae_anal(i, j), 1e-12)
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }

  Eigen::VectorXd Vec = lse_builder.GetLoadVector(2);
  Eigen::VectorXd Vec_anal(20);
  Vec_anal.setZero();

  ASSERT_EQ(Vec.rows(), Vec_anal.rows());
  for (int j = 0; j < Ae_anal.cols(); ++j) {
    EXPECT_NEAR(Vec(j), Vec_anal(j), 1e-12)
        << "mismatch in entry (" << j << ")";
  }

  // result verctor mu
  Eigen::VectorXcd mu_anal_two = Eigen::VectorXcd::Zero(8);
  Eigen::VectorXcd mu_two = std::get<1>(lse_builder.GetMuTwo());
  ASSERT_EQ(mu_two.rows(), mu_anal_two.rows());
  for (int j = 0; j < mu_anal_two.rows(); ++j) {
    EXPECT_NEAR(mu_two(j).imag(), mu_anal_two(j).imag(), 1e-12)
        << "mismatch in imag entry (" << j << ")";
    EXPECT_NEAR(mu_two(j).real(), mu_anal_two(j).real(), 1e-12)
        << "mismatch in real entry (" << j << ")";
  }

  Eigen::VectorXcd j_anal_two = Eigen::VectorXcd::Zero(12);
  Eigen::VectorXcd j_two = std::get<0>(lse_builder.GetMuTwo());
  ASSERT_EQ(j_two.rows(), j_anal_two.rows());
  for (int j = 0; j < mu_anal_two.rows(); ++j) {
    EXPECT_NEAR(j_two(j).imag(), j_anal_two(j).imag(), 1e-12)
        << "mismatch in imag entry (" << j << ")";
    EXPECT_NEAR(j_two(j).real(), j_anal_two(j).real(), 1e-12)
        << "mismatch in real entry (" << j << ")";
  }
}
