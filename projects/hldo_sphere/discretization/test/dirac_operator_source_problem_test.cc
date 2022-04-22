#include <dirac_operator_source_problem.h>
#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <complex>

TEST(projects_hldo_sphere_discretization,
     dirac_operator_source_problem_basic_test) {
  // Build LSE
  projects::hldo_sphere::discretization::DiracOperatorSourceProblem lse_builder;
  lse_builder.Compute();
  lse_builder.Solve();

  Eigen::SparseMatrix<std::complex<double>> Ae =
      lse_builder.GetGalerkinMatrix().makeSparse();

  Eigen::MatrixXd Ae_anal(26, 26);
  // clang-format off
  Ae_anal <<     
       0, 0, 0, 0, 0, 0,   -1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    1, 0, 0, 0,-1, 1,-1, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    0,-1, 0, 0, 1, 0, 0,-1, 1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    0, 0,-1, 0, 0, 0, 0, 1, 0,-1, 1, 0,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    0, 0, 0,-1, 0,-1, 0, 0, 0, 1, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 1, 0,-1, 0,-1,-1,   0, 0, 0, 0, 0, 0, 0, 0,

       1,-1, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  -1, 0, 0, 0, 0, 0, 1, 0,
      -1, 0, 1, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  -1, 0, 1, 0, 0, 0, 0, 0,
      -1, 0, 0, 1, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0,-1, 0, 1, 0, 0, 0,
      -1, 0, 0, 0, 1, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0,-1, 0, 1, 0,
       0, 1,-1, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  -1, 1, 0, 0, 0, 0, 0, 0,
       0,-1, 0, 0, 1, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0,-1, 1,
       0, 1, 0, 0, 0,-1,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0,-1, 0, 0, 0, 0, 0, 1,
       0, 0, 1,-1, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0,-1, 1, 0, 0, 0, 0,
       0, 0,-1, 0, 0, 1,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0,-1, 0, 1, 0, 0, 0, 0,
       0, 0, 0, 1,-1, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0,-1, 1, 0, 0,
       0, 0, 0,-1, 0, 1,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0,-1, 0, 1, 0, 0,
       0, 0, 0, 0,-1, 1,    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0,-1, 0, 1,

       0, 0, 0, 0, 0, 0,   -1,-1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 1, 0,-1, 0,-1, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    0, 1,-1, 0, 0, 0, 0,-1, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 1, 1, 0,-1, 0,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    0, 0, 1,-1, 0, 0, 0, 0, 0,-1, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    1, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0,    0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1,   0, 0, 0, 0, 0, 0, 0, 0;
  // clang-format on

  Ae_anal.block(0, 6, 6, 12) *= 1. / std::sqrt(3);
  Ae_anal.block(6, 0, 12, 6) *= -1. / std::sqrt(3);
  Ae_anal.block(6, 18, 12, 8) *= -1;
  Ae_anal.block(18, 6, 8, 12) *= -1;
  Eigen::MatrixXcd Ae_anal_cd = Ae_anal.cast<std::complex<double>>();

  Eigen::MatrixXd Ae_anal_mass_zero(6, 6);
  // clang-format off
  Ae_anal_mass_zero <<     
                 8,  2,  2,  2,  2,  0,
                 2,  8,  2,  0,  2,  2,  
                 2,  2,  8,  2,  0,  2,  
                 2,  0,  2,  8,  2,  2,  
                 2,  2,  0,  2,  8,  2,
                 0,  2,  2,  2,  2,  8;
  // clang-format on
  Ae_anal_mass_zero *= 1. / 8. / sqrt(3);
  Ae_anal_cd.topLeftCorner(6, 6) +=
      Ae_anal_mass_zero.cast<std::complex<double>>() *
      std::complex<double>(0., 1.);

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
  Ae_anal_cd.block(6, 6, 12, 12) +=
      Ae_anal_mass_one.cast<std::complex<double>>() *
      std::complex<double>(0., 1.);

  Eigen::MatrixXd Ae_anal_mass_two = Eigen::MatrixXd::Identity(8, 8);
  Ae_anal_mass_two *= sqrt(3) / 2.;

  Ae_anal_cd.bottomRightCorner(8, 8) +=
      Ae_anal_mass_two.cast<std::complex<double>>() *
      std::complex<double>(0., 1.);

  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal_cd.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal_cd.cols());
  for (int i = 0; i < Ae_anal_cd.rows(); ++i) {
    for (int j = 0; j < Ae_anal_cd.cols(); ++j) {
      EXPECT_NEAR(Ae.coeff(i, j).imag(), Ae_anal_cd(i, j).imag(), 1e-12)
          << "mismatch in imag entry (" << i << ", " << j << ")";
      EXPECT_NEAR(Ae.coeff(i, j).real(), Ae_anal_cd(i, j).real(), 1e-12)
          << "mismatch in real entry (" << i << ", " << j << ")";
    }
  }

  // Check rigthhandside
  Eigen::VectorXcd vec = lse_builder.GetLoadVector();
  Eigen::VectorXd vec_anal(26);
  vec_anal.setZero();
  Eigen::VectorXcd vec_anal_cd = vec_anal.cast<std::complex<double>>();
  ASSERT_EQ(vec.rows(), vec_anal_cd.rows());
  for (int j = 0; j < Ae_anal.rows(); ++j) {
    EXPECT_NEAR(vec(j).imag(), vec_anal_cd(j).imag(), 1e-12)
        << "mismatch in imag entry (" << j << ")";
    EXPECT_NEAR(vec(j).real(), vec_anal_cd(j).real(), 1e-12)
        << "mismatch in real entry (" << j << ")";
  }

  // Check result mu
  Eigen::VectorXcd mu_anal = Eigen::VectorXcd::Zero(26);
  Eigen::VectorXcd mu = lse_builder.GetMu();
  ASSERT_EQ(mu.rows(), mu_anal.rows());
  for (int j = 0; j < mu_anal.rows(); ++j) {
    EXPECT_NEAR(mu(j).imag(), mu_anal(j).imag(), 1e-12)
        << "mismatch in imag entry (" << j << ")";
    EXPECT_NEAR(mu(j).real(), mu_anal(j).real(), 1e-12)
        << "mismatch in real entry (" << j << ")";
  }
}
