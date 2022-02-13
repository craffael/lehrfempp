/**
 * @file
 * @brief Data structures representing HP finite elements
 * @author Tobias Rohner
 * @date May 2020
 * @copyright MIT License
 */

#include "hierarchic_fe.h"

namespace lf::fe {

double Legendre(unsigned n, double x, double t) {
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

double ILegendre(unsigned n, double x, double t) {
  if (n == 1) {
    return x;
  }
  return (Legendre(n, x, t) - t * t * Legendre(n - 2, x, t)) /
         (2 * (2 * n - 1));
}

double ILegendreDx(unsigned n, double x, double t) {
  return Legendre(n - 1, x, t);
}

double ILegendreDt(unsigned n, double x, double t) {
  if (n == 1) {
    return 0;
  }
  return -0.5 * (Legendre(n - 1, x, t) + t * Legendre(n - 2, x, t));
}

double LegendreDx(unsigned n, double x, double t) {
  if (n == 0) {
    return 0;
  }
  // This has quadratic complexity and could be improved to linear
  // but I don't think that's necessary
  return 2 * n * Legendre(n - 1, x, t) + (2 * x - t) * LegendreDx(n - 1, x, t);
}

double Jacobi(unsigned n, double alpha, double beta, double x) {
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

double Jacobi(unsigned n, double alpha, double x) {
  return Jacobi(n, alpha, 0, x);
}

double IJacobi(unsigned n, double alpha, double x) {
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

double IJacobiDx(unsigned n, double alpha, double x) {
  return Jacobi(n - 1, alpha, x);
}

double JacobiDx(unsigned n, double alpha, double x) {
  if (n == 0) {
    return 0;
  }
  return Jacobi(n - 1, alpha + 1, 1, x) * (alpha + n + 1);
}

}  // end namespace lf::fe
