#include <dirac_operator.h>
#include <gtest/gtest.h>

#include <array>
#include <cmath>

/**
 * @brief Test the Galerkin LSE for the Dirac Operator and its load vector
 *
 * @f[
 *   \begin{pmatrix}
 *       & \int_{\partial \mathbb{S}} \bm{u} \ grad_{\Gamma} v \, dS & \li
 *       \int_{\partial \mathbb{S}} grad_{\Gamma} u \cdot \bm{v} \, dS
 *       & &
 *       \int_{\partial \mathbb{S}} \mu \ curl_{\Gamma} \bm{v} \, dS  \li
 *       & \int_{\partial \mathbb{S}} curl_{\Gamma} \bm{u} \ \nu \, dS &
 *   \end{pmatrix}
 *   &=
 *   \begin{pmatrix}
 *      \int_{\partial \mathbb{S}} f v \, dS \li
 *      \int_{\partial \mathbb{S}} \bm{f} \cdot \bm{v} \, dS \li
 *      \int_{\partial \mathbb{S}} \varphi \nu \, dS
 *   \end{pmatrix}
 *   @f]
 *
 *
 * On the octaeder with radius 1
 *
 * And load function equal to 0
 *
 */
TEST(projects_hldo_sphere_operators, dirac_operator_test) {
  // Build LSE

  projects::hldo_sphere::operators::DiracOperator lse_builder;
  lse_builder.Compute();

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

  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae.coeff(i, j).real(), Ae_anal_cd(i, j).real(), 1e-12)
          << "mismatch in real entry (" << i << ", " << j << ")";
      EXPECT_NEAR(Ae.coeff(i, j).imag(), Ae_anal_cd(i, j).imag(), 1e-12)
          << "mismatch in imag entry (" << i << ", " << j << ")";
    }
  }

  Eigen::VectorXcd Vec = lse_builder.GetLoadVector();
  Eigen::VectorXd Vec_anal(26);
  Vec_anal.setZero();

  Eigen::VectorXcd Vec_anal_cd = Vec_anal.cast<std::complex<double>>();
  for (int j = 0; j < Ae_anal.cols(); ++j) {
    EXPECT_NEAR(Vec(j).real(), Vec_anal_cd(j).real(), 1e-12)
        << "mismatch in real entry (" << j << ")";
    EXPECT_NEAR(Vec(j).imag(), Vec_anal_cd(j).imag(), 1e-12)
        << "mismatch in imag entry (" << j << ")";
  }
}
