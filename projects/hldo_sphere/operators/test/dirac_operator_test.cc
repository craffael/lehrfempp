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
 * And load function equal to
 * @f[
 *   f_0 = 5.2 \imath \\
 *   f_1 =
 *   &=
 *   \begin{pmatrix}
 *      -y // x // 0
 *   \end{pmatrix}\\
 *    f_2 = 5.2 \imath - \frac{2 z}{\sqrt{x^2 + y^2 + z^2}}
 *   @f]
 *
 * Computations of the analytical solutions are in the file
 * `DiracOperatorTest.nb`
 *
 */
TEST(projects_hldo_sphere_operators, dirac_operator_test) {
  // Build LSE

  double k = 1;

  // righthandside for the zero
  auto f_zero = [&](const Eigen::Vector3d &x_vec) -> std::complex<double> {
    return std::complex<double>(0, k * 5.2);
  };

  // righthandside for the one form
  auto f_one = [&](const Eigen::Vector3d &x_vec) -> Eigen::VectorXcd {
    Eigen::Vector3cd ret;
    ret << std::complex<double>(0, -k * x_vec(1)),
        std::complex<double>(0, k * x_vec(0)), 0;
    return ret;
  };

  auto Power = [](double a, double b) -> double { return pow(a, b); };
  auto Sqrt = [](double a) -> double { return sqrt(a); };
  auto Complex = [](double a, double b) -> std::complex<double> {
    return std::complex<double>(a, b);
  };

  // righthandside for the two form
  auto f_two = [&](const Eigen::Vector3d &x_vec) -> std::complex<double> {
    double x = x_vec(0);
    double y = x_vec(1);
    double z = x_vec(2);
    return Complex(0., 5.2) * k -
           (2 * z) / Sqrt(Power(x, 2) + Power(y, 2) + Power(z, 2));
  };

  projects::hldo_sphere::operators::DiracOperator lse_builder;
  lse_builder.SetLoadFunctions(f_zero, f_one, f_two);
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
  Ae_anal.block(6, 18, 12, 8) *= 1;
  Ae_anal.block(18, 6, 8, 12) *= 1;

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
  Eigen::VectorXcd Vec_anal(26);
  Vec_anal.setZero();
  Vec_anal.segment(0, 6) << Eigen::VectorXcd::Ones(6) *
                                std::complex<double>(0, 6.004442799572108);
  Vec_anal.segment(6, 12) << 0, 0, 0, 0,
      Complex(0, 0.8333333333333334) / Sqrt(3),
      Complex(0, 0.8333333333333334) / Sqrt(3), 0,
      Complex(0, 0.8333333333333334) / Sqrt(3), 0,
      Complex(0, 0.8333333333333334) / Sqrt(3), 0, 0;
  Vec_anal.segment(18, 8) << std::complex<double>(-0.8164965809277259,
                                                  4.50333209967908),
      std::complex<double>(0.8164965809277259, 4.50333209967908),
      std::complex<double>(-0.8164965809277259, 4.50333209967908),
      std::complex<double>(0.8164965809277259, 4.50333209967908),
      std::complex<double>(-0.8164965809277259, 4.50333209967908),
      std::complex<double>(0.8164965809277259, 4.50333209967908),
      std::complex<double>(-0.8164965809277259, 4.50333209967908),
      std::complex<double>(0.8164965809277259, 4.50333209967908);

  for (int j = 0; j < Ae_anal.cols(); ++j) {
    EXPECT_NEAR(Vec(j).real(), Vec_anal(j).real(), 1e-12)
        << "mismatch in real entry (" << j << ")";
    EXPECT_NEAR(Vec(j).imag(), Vec_anal(j).imag(), 1e-12)
        << "mismatch in imag entry (" << j << ")";
  }
}
