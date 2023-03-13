/**
 * @file
 * @brief Data structures representing HP finite elements
 * @author Tobias Rohner
 * @date May 2020
 * @copyright MIT License
 */

#include "hierarchic_fe.h"

namespace lf::fe {

/**
 * @brief computes the `n`-th degree scaled Legendre Polynomial \f$ P_n(x;t)
 *\f$
 * @param n The degree of the polynomial
 * @param x The evaluation coordinate
 * @param t The scaling parameter
 *
 * To evaluate the scaled Legendre Polynomials \f$ P_n(x;t) \f$, we use that
 * \f[ \begin{aligned}
 *	    P_0(x;t) &= 1 \\
 *	    P_1(x;t) &= 2x - t \\
 *	    nP_n(x;t) &= (2n-1)(2x-t)P_{n-1}(x;t) - (n-1)t^2P_{p-2}(x;t)
 * \end{aligned} \f]
 */
double legendre(unsigned n, double x, double t) {
  if (n == 0) {
    return 1;
  }
  if (n == 1) {
    return 2 * x - t;
  }
  double Pim1 = 1;
  double Pi = 2 * x - t;
  for (unsigned i = 2; i <= n; ++i) {
    const double Pip1 =
        (2 * i - 1) * (2 * x - t) / i * Pi - (i - 1) * t * t / i * Pim1;
    Pim1 = Pi;
    Pi = Pip1;
  }
  return Pi;
}

/**
 * @brief computes the integral of the (n-1)-th degree scaled Legendre
 *Polynomial
 * @param n The degree of the integrated polynomial
 * @param x The evaluation coordinate
 * @param t The scaling parameter
 *
 * The integral is evaluated using
 * \f[
 *	\begin{aligned}
 *	    L_1(x) &= x \\
 *	    2(2n-1)L_n(x) &= P_n(x) - t^2P_{n-2}(x)
 *	\end{aligned}
 * \f]
 */
double ilegendre(unsigned n, double x, double t) {
  if (n == 1) {
    return x;
  }
  return (legendre(n, x, t) - t * t * legendre(n - 2, x, t)) /
         (2 * (2 * n - 1));
}

/**
 * @brief Computes \f$ \frac{\partial}{\partial x} L(x;t) \f$
 * @param n The degree of the integrated scaled Legendre polynomial
 * @param x The evaluation coordinate
 * @param t The scaling parameter
 *
 * The derivative is simply given by \f$ \frac{\partial}{\partial x} L_n(x;t) =
 * P_{n-1}(x;t) \f$
 */
double ilegendre_dx(unsigned n, double x, double t) {
  return legendre(n - 1, x, t);
}

/**
 * @brief Computes \f$ \frac{\partial}{\partial t} L(x;t) \f$
 * @param n The degree of the integrated scaled Legendre polynomial
 * @param x The evaluation coordinate
 * @param t The scaling parameter
 *
 * The derivative is given by
 * \f[
 *	\begin{aligned}
 *	    \frac{\partial}{\partial t} L_1(x;t) &= 0 \\
 *	    \frac{\partial}{\partial t} L_n(x;t) &= -\frac{1}{2} \left(
 *P_{n-1}(x;t) + tP_{n-2}(x;t) \right) \end{aligned} \f]
 */
double ilegendre_dt(unsigned n, double x, double t) {
  if (n == 1) {
    return 0;
  }
  return -0.5 * (legendre(n - 1, x, t) + t * legendre(n - 2, x, t));
}

/**
 * @brief Computes the derivative of the n-th degree scaled Legendre polynomial
 * @param n The degree of the polynomial
 * @param x The evaluation coordinate
 * @param t The scaling parameter
 *
 * The derivative is given by
 * \f[
 *	\begin{aligned}
 *	    \frac{\partial}{\partial x} L_0(x;t) &= 0 \\
 *	    \frac{\partial}{\partial x} L_n(x;t) &= 2nP_{n-1}(x;t) +
 *(2x-t)\frac{\partial}{\partial x}P_{n-1}(x;t) \\ \end{aligned} \f]
 */
// NOLINTNEXTLINE(misc-no-recursion)
double legendre_dx(unsigned n, double x, double t) {
  if (n == 0) {
    return 0;
  }
  // This has quadratic complexity and could be improved to linear
  // but I don't think that's necessary
  return 2 * n * legendre(n - 1, x, t) + (2 * x - t) * legendre_dx(n - 1, x, t);
}

/**
 * @brief Computes the n-th degree shifted Jacobi polynomial
 * @param n The degree of the polynomial
 * @param alpha The \f$ \alpha \f$ parameter of the Jacobi polynomial
 * @param beta The \f$ \beta \f$ parameter of the Jacobi polynomial
 * @param x The evaluation coordinate
 *
 * We use the recurrence relation for non-shifted Jacobi polynomials
 * \f[
 *	\begin{aligned}
 *	    P_0^{(\alpha,\beta)}(x) &= 1 \\
 *	    P_1^{(\alpha,\beta)}(x) &= \frac{1}{2} \left( \alpha - \beta +
 *(\alpha + \beta + 2)x \right) \\ P_{n+1}^{(\alpha,\beta)}(x) &= \frac{1}{a_n}
 *\left( (b_n+c_nx)P_n^{(\alpha,\beta)}(x) - d_nP_{n-1}^{(\alpha,\beta)}(x)
 *\right) \end{aligned} \f] where \f[ \begin{aligned}
 *	    a_n &= 2(n+1)(n+\alpha+\beta+1)(2n+\alpha+\beta) \\
 *	    b_n &= (2n+\alpha+\beta+1)(\alpha^2-\beta^2) \\
 *	    c_n &= (2n+\alpha+\beta)(2n+\alpha+\beta+1)(2n+\alpha+\beta+2) \\
 *	    d_n &= 2(n+\alpha)(n+\beta)(2n+\alpha+\beta+2)
 *	\end{aligned}
 * \f]
 */
double jacobi(unsigned n, double alpha, double beta, double x) {
  // The recurrence relation is for the non-shifted Jacobi Polynomials
  // We thus map [0, 1] onto [-1, 1]
  x = 2 * x - 1;
  if (n == 0) {
    return 1;
  }
  if (n == 1) {
    return 0.5 * (alpha - beta + x * (alpha + beta + 2));
  }
  double p0 = 1;
  double p1 = 0.5 * (alpha - beta + x * (alpha + beta + 2));
  double p2 = 0;
  for (unsigned j = 1; j < n; ++j) {
    const double alpha1 =
        2 * (j + 1) * (j + alpha + beta + 1) * (2 * j + alpha + beta);
    const double alpha2 =
        (2 * j + alpha + beta + 1) * (alpha * alpha - beta * beta);
    const double alpha3 = (2 * j + alpha + beta) * (2 * j + alpha + beta + 1) *
                          (2 * j + alpha + beta + 2);
    const double alpha4 =
        2 * (j + alpha) * (j + beta) * (2 * j + alpha + beta + 2);
    p2 = 1. / alpha1 * ((alpha2 + alpha3 * x) * p1 - alpha4 * p0);
    p0 = p1;
    p1 = p2;
  }
  return p2;
}

/**
 * @brief Computes the n-th degree shifted Jacobi polynomial for \f$ \beta = 0
 * \f$
 * @param n The degree of the polynomial
 * @param alpha The \f$ \alpha \f$ parameter of the Jacobi polynomial
 * @param x The evaluation coordinate
 */
double jacobi(unsigned n, double alpha, double x) {
  return jacobi(n, alpha, 0, x);
}

/**
 * @brief Evaluate the integral of the (n-1)-th degree Jacobi Polynomial for \f$
 *\beta = 0 \f$
 * @param n The degree of the integrated polynomial
 * @param alpha The \f$ \alpha \f$ parameter of the Jacobi polynomial
 * @param x The evaluation coordinate
 *
 * The integral is evaluated using
 * \f[
 *	\begin{aligned}
 *	    L_1^\alpha(x) &= x \\
 *	    L_p^\alpha(x) &= a_pP_p^\alpha(x) + b_pP_{p-1}^\alpha(x) -
 *c_pP_{p-2}^\alpha(x) \end{aligned} \f] where the coefficients are defined as
 * \f[
 *	\begin{aligned}
 *	    a_p &= \frac{p+\alpha}{(2p+\alpha-1)(2p+\alpha)} \\
 *	    b_p &= \frac{\alpha}{(2p+\alpha-1)(2p+\alpha)} \\
 *	    c_p &= \frac{p-1}{(2p+\alpha-2)(2p+\alpha-1)}
 *	\end{aligned}
 * \f]
 */
double ijacobi(unsigned n, double alpha, double x) {
  if (n == 0) {
    return -1;
  }
  if (n == 1) {
    return x;
  }
  // Compute the n-th, (n-1)-th and (n-2)-th Jacobi Polynomial
  // This uses the same recurrence relation as the `eval` method
  double ajp1P = 2 * 2 * (2 + alpha) * (2 * 2 + alpha - 2);
  double bjp1P = 2 * 2 + alpha - 1;
  double cjp1P = (2 * 2 + alpha) * (2 * 2 + alpha - 2);
  double djp1P = 2 * (2 + alpha - 1) * (2 - 1) * (2 * 2 + alpha);
  double Pjm2 = 1;
  double Pjm1 = (2 + alpha) * x - 1;
  double Pj =
      (bjp1P * (cjp1P * (2 * x - 1) + alpha * alpha) * Pjm1 - djp1P * Pjm2) /
      ajp1P;
  for (unsigned j = 2; j < n; ++j) {
    ajp1P = 2 * (j + 1) * ((j + 1) + alpha) * (2 * (j + 1) + alpha - 2);
    bjp1P = 2 * (j + 1) + alpha - 1;
    cjp1P = (2 * (j + 1) + alpha) * (2 * (j + 1) + alpha - 2);
    djp1P = 2 * ((j + 1) + alpha - 1) * ((j + 1) - 1) * (2 * (j + 1) + alpha);
    double Pjp1 =
        (bjp1P * (cjp1P * (2 * x - 1) + alpha * alpha) * Pj - djp1P * Pjm1) /
        ajp1P;
    Pjm2 = Pjm1;
    Pjm1 = Pj;
    Pj = Pjp1;
  }
  // Compute the integral by the formula in the docstring
  const double anL = (n + alpha) / ((2 * n + alpha - 1) * (2 * n + alpha));
  const double bnL = alpha / ((2 * n + alpha - 2) * (2 * n + alpha));
  const double cnL = (n - 1) / ((2 * n + alpha - 2) * (2 * n + alpha - 1));
  return anL * Pj + bnL * Pjm1 - cnL * Pjm2;
}

/**
 * @brief Computes the derivative of the n-th integrated scaled Jacobi
 * polynomial
 * @param n The degree of the integrated scaled Jacobi polynomial
 * @param alpha The \f$ \alpha \f$ parameter of the Jacobi polynomial
 * @param x The evaluation coordinate
 *
 * The derivative is simply given by \f$ \frac{\partial}{\partial x}
 * L_n^{(\alpha,0)}(x) = P_{n-1}^{(\alpha,0)}(x) \f$
 */
double ijacobi_dx(unsigned n, double alpha, double x) {
  return jacobi(n - 1, alpha, x);
}

/**
 * @brief Computes the derivative of the n-th degree Jacobi Polynomial for \f$
 *\beta = 0 \f$
 * @param n The degree of the differentiated polynomial
 * @param alpha The \f$ \alpha \f$ parameter of the Jacobi Polynomial
 * @param x The evaluation coordinate
 *
 * The derivative is evaluated using
 * \f[
 *	{P^{(\alpha,0)}_n}'(x) = \frac{\alpha+n+1}{2} P^{(\alpha+1,1)}_{n-1}(x)
 * \f]
 */
double jacobi_dx(unsigned n, double alpha, double x) {
  if (n == 0) {
    return 0;
  }
  return jacobi(n - 1, alpha + 1, 1, x) * (alpha + n + 1);
}

}  // end namespace lf::fe
