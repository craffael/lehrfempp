#ifndef LF_USCALFE_HP_FE_H_
#define LF_USCALFE_HP_FE_H_

/**
 * @file
 * @brief Data structures representing HP finite elements
 * @author Tobias Rohner
 * @date May 2020
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>
#include <lf/quad/quad.h>

#include <cmath>
#include <memory>
#include <vector>

#include "scalar_reference_finite_element.h"

namespace lf::fe {
/** Type for indices into global matrices/vectors */
using gdof_idx_t = lf::assemble::gdof_idx_t;
/** Type for indices referring to entity matrices/vectors */
using ldof_idx_t = lf::assemble::ldof_idx_t;
/** Type for vector length/matrix sizes */
using size_type = lf::assemble::size_type;
/** Type for (co-)dimensions */
using dim_t = lf::assemble::dim_t;
/** Type for global index of entities */
using glb_idx_t = lf::assemble::glb_idx_t;
/** Type for indexing sub-entities */
using sub_idx_t = lf::base::sub_idx_t;

/**
 * @brief Struct for computing Legrendre and integrated Legendre polynomials
 * on the interval [0, 1]
 */
template <typename SCALAR>
struct LegendrePoly {
  /**
   * @brief computes the `n`-th degree shifted Legendre Polynomial \f$ P_n(x)
   *\f$
   * @param n The degree of the polynomial
   * @param x The evaluation coordinate in [0, 1]
   *
   * To evaluate the non-shifted Legendre Polynomials \f$ \tilde{P}_p(x) \f$ on
   *[-1, 1] we use that \f[ \begin{aligned}
   *	    \tilde{P}_0(x) &= 1 \\
   *	    \tilde{P}_1(x) &= x \\
   *	    (p+1)\tilde{P}_{p+1}(x) &= (2p+1)\tilde{P}_p(x) -
   *p\tilde{P}_{p-1}(x) \end{aligned} \f]
   */
  static SCALAR eval(unsigned n, double x) {
    // The recurrence relation is for the non-shifted Legendre Polynomials
    // We thus map [0, 1] onto [-1, 1]
    x = 2 * x - 1;
    double Ljm1 = 1;
    double Lj = x;
    if (n == 0) {
      return Ljm1;  // * std::sqrt(2 * n + 1);
    }
    if (n == 1) {
      return Lj;  // * std::sqrt(2 * n + 1);
    }
    // If the degree is > 1, we use the recurrence relation
    for (unsigned j = 1; j < n; ++j) {
      const double Ljp1 = ((2 * j + 1) * x * Lj - j * Ljm1) / (j + 1);
      Ljm1 = Lj;
      Lj = Ljp1;
    }
    return Lj;  // * std::sqrt(2 * n + 1);
  }

  /**
   * @brief computes the integral of the (n-1)-th degree shifted Legendre
   *Polynomial
   * @param n The degree of the integrated polynomial
   * @param x The evaluation coordinate in [0, 1]
   *
   * The integral is evaluated using
   * \f[
   *	\begin{aligned}
   *	    L_1(x) &= x \\
   *	    2(2p-1)L_p(x) &= P_p(x) - P_{p-2}(x)
   *	\end{aligned}
   * \f]
   */
  static SCALAR integral(unsigned n, double x) {
    if (n == 0) {
      return -1;
    }
    if (n == 1) {
      return x;  // * std::sqrt(2 * n - 1);
    }
    // Map to the interval [-1, 1]
    x = 2 * x - 1;
    // Compute the n-th and (n-2)-th degree shifted Legendre Polynomials
    double Ljm2 = 1;
    double Ljm1 = x;
    double Lj = (3 * x * x - 1) / 2;
    for (unsigned j = 2; j < n; ++j) {
      double Ljp1 = ((2 * j + 1) * x * Lj - j * Ljm1) / (j + 1);
      Ljm2 = Ljm1;
      Ljm1 = Lj;
      Lj = Ljp1;
    }
    // Compute the integral
    return (Lj - Ljm2) / (4 * n - 2);  // * std::sqrt(2 * n - 1);
  }

  /**
   * @brief Computes the derivative of the (n+1)-th degree shifted Legendre
   * Polynomial
   * @param n The degree of the differentiated polynomial
   * @param x The evaluation coordinate in [0, 1]
   *
   * The derivative is evaluated using
   * \f[
   *	P_n'(x) = \begin{cases}
   *	    n(n+1) &\mbox{ for } x = 1 \\
   *	    (-1)^{n+1}n(n+1) &\mbox{ for } x = 0 \\
   *	    \frac{n}{2x^2-2x}\left((2x-1)P_n(x)-P_{n-1}(x)\right)
   *	\end{cases}
   * \f]
   */
  static SCALAR derivative(unsigned n, double x) {
    if (n == 0) {
      return 2.;  // * std::sqrt(2 * n + 3);
    }
    // Special cases for x == 0 and x == 1
    if (x == 1) {
      return (n + 1) * (n + 2);  // * std::sqrt(2 * n + 3);
    }
    if (x == 0) {
      // Depends on whether the polynomial is even or odd
      return (n % 2 == 0 ? (n + 1.) * (n + 2.)
                         : -(n + 1.) * (n + 2.));  // * std::sqrt(2 * n + 3);
    }
    // Map to the interval [-1, 1]
    x = 2 * x - 1;
    // Compute the n-th and (n+1)-th degree shifted Legendre Polynomials
    double Ljm1 = 1;
    double Lj = x;
    for (unsigned j = 1; j < n + 1; ++j) {
      const double Ljp1 = ((2 * j + 1) * x * Lj - j * Ljm1) / (j + 1);
      Ljm1 = Lj;
      Lj = Ljp1;
    }
    // Compute the derivative of the (n+1)-th polynomial
    return (2 * n + 2) / (x * x - 1) *
           (x * Lj - Ljm1);  // * std::sqrt(2 * n + 3);
  }
};

/**
 * @brief Struct for computing Jacobi and integrated Jacobi polynomials with \f$
 * \beta=0 \f$
 */
template <typename SCALAR>
struct JacobiPoly {
  /**
   * @brief Evaluate the n-th degree Jacobi polynomial
   * @param n The degree of the polynomial
   * @param alpha The \f$ \alpha \f$ parameter to the Jacobi polynomial
   * @param x The evaluation point in [0, 1]
   *
   * For the evaluation, we use that
   * \f[
   *	\begin{aligned}
   *	    P_0^\alpha(x) &= 1 \\
   *	    P_1^\alpha(x) &= (2+\alpha)x - 1 \\
   *	    a_pP_p^\alpha(x) &= b_p(c_p(2x-1)+\alpha^2)P_{p-1}^\alpha(x) -
   *d_pP_{p-2}^\alpha(x) \end{aligned} \f] where the coefficients are defined as
   * \f[
   *	\begin{aligned}
   *	    a_p &= 2p(p+\alpha)(2p+\alpha-2) \\
   *	    b_p &= 2p + \alpha - 1 \\
   *	    c_p &= (2p+\alpha)(2p+\alpha-2) \\
   *	    d_p &= 2(p+\alpha-1)(p-1)(2p+\alpha)
   *	\end{aligned}
   * \f]
   */
  static SCALAR eval(unsigned n, double alpha, double x) {
    double Pjm1 = 1;
    double Pj = (2 + alpha) * x - 1;
    if (n == 0) {
      return Pjm1;
    }
    if (n == 1) {
      return Pj;
    }
    // For n >= 2 we need to use the recurrence relation
    for (unsigned j = 1; j < n; ++j) {
      const double ajp1 =
          2 * (j + 1) * ((j + 1) + alpha) * (2 * (j + 1) + alpha - 2);
      const double bjp1 = 2 * (j + 1) + alpha - 1;
      const double cjp1 = (2 * (j + 1) + alpha) * (2 * (j + 1) + alpha - 2);
      const double djp1 =
          2 * ((j + 1) + alpha - 1) * ((j + 1) - 1) * (2 * (j + 1) + alpha);
      double Pjp1 =
          (bjp1 * (cjp1 * (2 * x - 1) + alpha * alpha) * Pj - djp1 * Pjm1) /
          ajp1;
      Pjm1 = Pj;
      Pj = Pjp1;
    }
    return Pj;
  }

  /**
   * @brief Evaluate the integral of the (n-1)-th degree Jacobi Polynomial
   * @param n The degree of the integrated polynomial
   * @param alpha The \f$ \alpha \f$ parameter of the Jacobi polynomial
   * @param x The evaluation coordinate in [0, 1]
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
  static SCALAR integral(unsigned n, double alpha, double x) {
    if (n == 0) {
      return -1;
    }
    if (n == 1) {
      return x;
    }
    // Compute the n-th, (n-1)-th and (n-2)-th Jacobi Polynomial
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
   * @brief Computes the derivative of the (n+1)-th degree Jacobi Polynomial
   * @param n The degree of the differentiated polynomial
   * @param alpha The \f$ \alpha \f$ parameter of the Jacobi Polynomial
   * @param x The evaluation coordinate in [0, 1]
   *
   * The derivative is evaluated using
   * \f[
   *	{P^{(\alpha,0)}_n}'(x) = \frac{\alpha+n+1}{2} P^{(\alpha+1,1)}_{n-1}(x)
   * \f]
   */
  static SCALAR derivative(unsigned n, double alpha, double x) {
    return eval_general(n, alpha + 1, 1, x) * (alpha + n + 2);
  }

  /**
   * @brief Evaluate the n-th degree Jacobi polynomial
   * @param n The degree of the polynomial
   * @param alpha The \f$ \alpha \f$ parameter to the Jacobi polynomial
   * @param beta The \f$ \beta \f$ parameter to the Jacobi polynomial
   * @param x The evaluation point in [0, 1]
   *
   * This method might be more unstable than the `JacobiPoly<SCALAR>::eval`
   * function which is specialized for \f$ \beta = 0 \f$.
   */
  static SCALAR eval_general(unsigned n, double alpha, double beta, double x) {
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
      const double alpha3 = (2 * j + alpha + beta) *
                            (2 * j + alpha + beta + 1) *
                            (2 * j + alpha + beta + 2);
      const double alpha4 =
          2 * (j + alpha) * (j + beta) * (2 * j + alpha + beta + 2);
      p2 = 1. / alpha1 * ((alpha2 + alpha3 * x) * p1 - alpha4 * p0);
      p0 = p1;
      p1 = p2;
    }
    return p2;
  }
};

/**
 * @brief Computes Chebyshev interpolation nodes in [0, 1]
 * @param n Degree of the Chebyshev interpolation nodes
 * @returns An Eigen vector containing the interpolation nodes on [0, 1]
 */
Eigen::VectorXd chebyshevNodes(unsigned n);

/**
 * @headerfile lf/fe/fe.h
 * @brief Linear finite element on a point
 *
 * This is a specialization of ScalarReferenceFiniteElement for an entity
 * of dimension 0, which is exactly one scalar value. It is an ingredient
 * of all hierarchic finite element spaces.
 */
template <class SCALAR>
class FeHierarchicPoint : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  /**
   * @brief Create a new FeHierarchicPoint by specifying the degree of the shape
   * functions.
   * @param degree The degree of the shape function.
   */
  explicit FeHierarchicPoint(unsigned degree) : degree_(degree) {}

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kPoint();
  }
  [[nodiscard]] unsigned Degree() const override { return degree_; }
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t subidx) const override {
    LF_ASSERT_MSG(codim == 0, "Codim out of bounds");
    LF_ASSERT_MSG(subidx == 0, "subidx out of bounds.");
    return 1;
  }
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 0, "refcoords has too many rows.");
    return Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>::Ones(1, refcoords.cols());
  }
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd & /*refcoords*/) const override {
    LF_VERIFY_MSG(false, "gradients not defined in points of mesh.");
  }
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return Eigen::MatrixXd(0, 1);
  }
  [[nodiscard]] size_type NumEvaluationNodes() const override { return 1; }

 private:
  unsigned degree_;
};

/**
 * @headerfile lf/fe/fe.h
 * @brief Hierarchic Finite Elements of arbitrary degree on segments
 *
 * The shape functions associated with the vertices are given by
 * \f[
 *  \begin{align*}
 *	\widehat{b^\cdot}^1(\widehat{x}) &= 1 - \widehat{x} \\
 *	\widehat{b^\cdot}^2(\widehat{x}) &= \widehat{x}
 *  \end{align*}
 * \f]
 * and the interior basis functions associated with the segment itself
 * are given by the integrated shifted Legendre polynomials
 * \f[
 *  \widehat{b^-}^i(\widehat{x}) = L_i(\widehat{x}) = \int_0^{\widehat{x}}\!
 *P_{i-1}(\xi) \,\mathrm{d}\xi \quad\mbox{ for }\quad i= 2, \cdots, p \f] where
 *\f$P_i : [0, 1] \to \mathbb{R}\f$ is the shifted Legendre polynomial of degree
 *\f$i\f$.
 *
 * To compute the basis function coefficients from point evaluations of a
 *function, we make use of the dual basis given by \f[ \lambda_i^-[f] =
 *\begin{cases}
 *	f(0) &\mbox{ for } i = 0 \\
 *	f(1) &\mbox{ for } i = 1 \\
 *	\frac{1}{2i - 1} \left [P_{i-1}(1)f(1) - P_{i-1}(0)f(0) - \int_0^1\!
 *P_{i-1}'(x)f(x) \,\mathrm{d}x \right] &\mbox{ for } i \geq 2 \end{cases} \f]
 *
 * @attention Note that the local coordinate \f$\widehat{x}\f$ may be flipped by
 *applying an affine transformation \f$\widehat{x} \mapsto 1 - \widehat{x}\f$.
 * This must be done if the relative orientation of the edge is
 *`lf::mesh::Orientation::negative` in order to keep a global ordering of DOFs
 *on the cell interfaces.
 *
 * @note The Shape Functions are taken from the following paper:
 * https://arxiv.org/pdf/1504.03025.pdf
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeHierarchicSegment final : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeHierarchicSegment(const FeHierarchicSegment &) = default;
  FeHierarchicSegment(FeHierarchicSegment &&) noexcept = default;
  FeHierarchicSegment &operator=(const FeHierarchicSegment &) = default;
  FeHierarchicSegment &operator=(FeHierarchicSegment &&) noexcept = default;
  ~FeHierarchicSegment() override = default;

  FeHierarchicSegment(unsigned degree,
                      nonstd::span<const lf::mesh::Orientation> rel_orient)
      : ScalarReferenceFiniteElement<SCALAR>(),
        degree_(degree),
        qr_dual_(lf::quad::make_QuadRule(RefEl(), 2 * (degree - 1))),
        rel_orient_(rel_orient) {}

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kSegment();
  }

  [[nodiscard]] unsigned Degree() const override { return degree_; }

  /**
   * @brief The local shape functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
    return degree_ + 1;
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions for the
   * segment
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim) const override {
    return codim == 0 ? degree_ - 1 : 1;
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions for the
   * segment
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,
   * sub_idx_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1, "refcoords must be a row vector");
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        NumRefShapeFunctions(), refcoords.cols());
    // Get the shape functions associated with the vertices
    result.row(0) =
        refcoords.unaryExpr([&](double x) -> SCALAR { return 1 - x; });
    result.row(1) = refcoords.unaryExpr([&](double x) -> SCALAR { return x; });
    // Get the shape functions associated with the interior of the segment
    for (int i = 0; i < degree_ - 1; ++i) {
      result.row(i + 2) = refcoords.unaryExpr([&](double x) -> SCALAR {
        return LegendrePoly<SCALAR>::integral(i + 2, x);
      });
    }
    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd &refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1, "refcoords must be a row vector");
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        NumRefShapeFunctions(), refcoords.cols());
    // Get the gradient of the shape functions associated with the vertices
    result.row(0) =
        refcoords.unaryExpr([&](double /*x*/) -> SCALAR { return -1; });
    result.row(1) =
        refcoords.unaryExpr([&](double /*x*/) -> SCALAR { return 1; });
    // Get the shape functions associated with the interior of the segment
    for (int i = 0; i < degree_ - 1; ++i) {
      result.row(i + 2) = refcoords.unaryExpr([&](double x) -> SCALAR {
        return LegendrePoly<SCALAR>::eval(i + 1, x);
      });
    }
    return result;
  }

  /**
   * @brief Evaluation nodes are the endpoints of the segment and the Chebyshev
   * nodes of degree p-1 on the segment
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    // First two nodes are vertices of the segment
    Eigen::MatrixXd nodes(1, NumEvaluationNodes());
    nodes(0, 0) = 0;
    nodes(0, 1) = 1;
    // The other nodes are quadrature points
    nodes.block(0, 2, 1, qr_dual_.NumPoints()) = qr_dual_.Points();
    return nodes;
  }

  /**
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    return qr_dual_.NumPoints() + 2;
  }

  /**
   * @brief Maps function evaluations to basis function coefficients
   * @param nodevals The value of the function at the evaluation nodes
   *
   * This function computes the basis function coefficients using the dual
   * basis of the basis functions on the segment. It is given by
   * \f[
   *  \lambda_i^-[f] = \begin{cases}
   *	f(0) &\mbox{ for } i = 0 \\
   *	f(1) &\mbox{ for } i = 1 \\
   *	\frac{1}{2i - 1} \left [P_{i-1}(1)f(1) - P_{i-1}(0)f(0) - \int_0^1\!
   *P_{i-1}'(x)f(x) \,\mathrm{d}x \right] &\mbox{ for } i \geq 2 \end{cases} \f]
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const override {
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(NumRefShapeFunctions());
    // Compute the first and second basis function coefficients
    dofs[0] = nodevals[0];
    dofs[1] = nodevals[1];
    // Compute the other basis function coefficients
    for (lf::base::size_type i = 2; i < NumRefShapeFunctions(); ++i) {
      const SCALAR P0 = LegendrePoly<SCALAR>::eval(i - 1, 0);
      const SCALAR P1 = LegendrePoly<SCALAR>::eval(i - 1, 1);
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
          qr_dual_.Points().unaryExpr([&](double x) -> SCALAR {
            return LegendrePoly<SCALAR>::derivative(i - 2, x);
          });
      // Evaluate the integral from the dual basis
      const SCALAR integ =
          (qr_dual_.Weights().transpose().array() * psidd.array() *
           nodevals.tail(qr_dual_.NumPoints()).array())
              .sum();
      // Compute the basis function coefficient
      dofs[i] = (P1 * nodevals[1] - P0 * nodevals[0] - integ) * (2 * i - 1);
    }
    return dofs;
  }

 private:
  unsigned degree_;
  lf::quad::QuadRule qr_dual_;
  nonstd::span<const lf::mesh::Orientation> rel_orient_;
};

/**
 * @headerfile lf/fe/fe.h
 * @brief Hierarchic Finite Elements of arbitrary degree on triangles
 *
 * The shape functions associated with the vertices are given by the barycentric
 *coordinates \f[ \begin{align*}
 *	\widehat{b^{\cdot}}^1(\widehat{x}, \widehat{y}) &= \lambda_1 = 1 -
 *\widehat{x} - \widehat{y} \\
 *	\widehat{b^{\cdot}}^2(\widehat{x}, \widehat{y}) &= \lambda_2 =
 *\widehat{x} \\ \widehat{b^{\cdot}}^3(\widehat{x}, \widehat{y}) &= \lambda_3 =
 *\widehat{y} \end{align*} \f]
 *
 * The basis functions associated with the triangle edges are given by the
 *homogenized integrated Legendre polynomials \f[ \begin{align*}
 *	\widehat{b^-}^i(\widehat{x}, \widehat{y}) &= (\lambda_1 + \lambda_2)^i
 *L_i\left(\frac{\lambda_2}{\lambda_1+\lambda_2}\right) &\mbox{ for edge 1 and }
 *i = 2, \cdots, p \\
 *	\widehat{b^-}^i(\widehat{x}, \widehat{y}) &= (\lambda_2 + \lambda_3)^i
 *L_i\left(\frac{\lambda_3}{\lambda_2+\lambda_3}\right) &\mbox{ for edge 2 and }
 *i = 2, \cdots, p \\ \widehat{b^-}^i(\widehat{x}, \widehat{y}) &= (\lambda_3 +
 *\lambda_1)^i L_i\left(\frac{\lambda_1}{\lambda_3+\lambda_1}\right) &\mbox{ for
 *edge 3 and } i = 2, \cdots, p \end{align*} \f] Note that the basis function on
 *a specific edge is always zero on the other two edges. This is needed to
 *guarantee continuity of the function space.
 *
 * The basis functions associated with the interior of the triangle are given by
 *the edge basis functions multiplied with an integrated Jacobi polynomial to
 *force the value of the basis function to be zero on all edges. \f[
 *  \widehat{b^{\triangle}}^{ij}(\widehat{x}, \widehat{y}) = (\lambda_1 +
 *\lambda_2)^i L_i\left(\frac{\lambda_2}{\lambda_1+\lambda_2}\right)
 *L_j^{2i}(\lambda_3) \quad\mbox{ for } i \geq 2, j \geq 1, i+j = 3, \cdots, p
 * \f]
 * where \f$L_i : [0, 1] \to \mathbb{R}\f$ and \f$L_i^{\alpha} : [0, 1] \to
 *\mathbb{R}\f$ are the integrated shifted Legendre and integrated shifted
 *Jacobi polynomials respectively: \f[ \begin{align*}
 *	L_i(x) &= \int_0^x\! P_{i-1}(\xi) \,\mathrm{d}\xi \\
 *	L_i^{\alpha}(x) &= \int_0^x\! P_{i-1}^{(\alpha, 0)}(\xi) \,\mathrm{d}\xi
 *  \end{align*}
 * \f]
 *
 * To compute the basis function coefficientsfrom point evaluations of a
 *function, we make use of the dual basis. For the vertices it is simply given
 *by \f[ \lambda_i^{\cdot}[f] = \begin{cases}
 *	f(0, 0) &\mbox{ for } i = 1 \\
 *	f(1, 0) &\mbox{ for } i = 2 \\
 *	f(0, 1) &\mbox{ for } i = 3
 *  \end{cases}
 * \f]
 * For the dual basis on the edges, we simply apply the segment dual basis along
 *the edges of the triangle \f[ \lambda_i^-[f] = \begin{cases}
 *	\frac{1}{2i-1} \left[ P_{i-1}(1)f(1, 0) - P_{i-1}(0)f(0, 0) - \int_0^1\!
 *P_{i-1}'(x)f(x, 0) \,\mathrm{d}x \right] &\mbox{ for edge 1} \\
 *	\frac{1}{2i-1} \left[ P_{i-1}(1)f(0, 1) - P_{i-1}(0)f(1, 0) - \int_0^1\!
 *P_{i-1}'(x)f(1-x, x) \,\mathrm{d}{x} \right] &\mbox{ for edge 2} \\
 *	\frac{1}{2i-1} \left[ P_{i-1}(1)f(0, 0) - P_{i-1}(0)f(0, 1) - \int_0^1\!
 *P_{i-1}'(x)f(0, 1-x) \,\mathrm{d}x \right] &\mbox{ for edge 2} \end{cases} \f]
 * The dual basis for the interior shape functions is quite a bit more involved.
 *It is given by \f[ \lambda_{ij}^{\triangle}[f] = \frac{1}{(2i-1)(2i+2j-1)}
 *\left[ \int_0^1\! \int_0^{1-y}\! f(x, y) \left( (1-y)^{i-1}{L_j^{2i}}''(y) -
 *2i(1-y)^{i-2}{L_j^{2i}}'(y) \right) L_{j+1}''\left(\frac{x}{1-y}\right)
 *\,\mathrm{d}x \,\mathrm{d}y \right] \f] and must additionally be
 *orthogonalized with respect to the dual basis on the vertices and edges of the
 *triangle by subtracting them accordingly.
 *
 * @attention Note that for the basis functions associated with the edges,
 *depending on the `lf::mesh::Orientation` of the according edge, the local
 *coordinate may be flipped to ensure continuity of the function space over the
 *cell interfaces of the mesh. The basis functions and the dual basis must be
 * adjusted accordingly in this case.
 *
 * @note The Shape Functions are taken from the following paper:
 * https://arxiv.org/pdf/1504.03025.pdf
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeHierarchicTria final : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeHierarchicTria(const FeHierarchicTria &) = default;
  FeHierarchicTria(FeHierarchicTria &&) noexcept = default;
  FeHierarchicTria &operator=(const FeHierarchicTria &) = default;
  FeHierarchicTria &operator=(FeHierarchicTria &&) noexcept = default;
  ~FeHierarchicTria() override = default;

  FeHierarchicTria(unsigned degree,
                   nonstd::span<const lf::mesh::Orientation> rel_orient)
      : ScalarReferenceFiniteElement<SCALAR>(),
        degree_(degree),
        qr_dual_segment_(lf::quad::make_QuadRule(lf::base::RefEl::kSegment(),
                                                 2 * (degree - 1))),
        qr_dual_tria_(lf::quad::make_QuadRule(lf::base::RefEl::kTria(),
                                              2 * (degree - 1))),
        rel_orient_(rel_orient) {}

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kTria();
  }

  [[nodiscard]] unsigned Degree() const override { return degree_; }

  /**
   * @brief The local shape functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
    return (degree_ + 1) * (degree_ + 2) / 2;
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and max(0, (p-2)*(p-1)/2) shape functions on the triangle
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim) const override {
    switch (codim) {
      case 0:
        if (degree_ <= 2) {
          return 0;
        } else {
          return (degree_ - 2) * (degree_ - 1) / 2;
        }
      case 1:
        return degree_ - 1;
      case 2:
        return 1;
      default:
        LF_ASSERT_MSG(false, "Illegal codim " << codim);
        // Silence compiler warnings
        return 0;
    }
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and max(0, (p-2)*(p-1)/2) shape functions on the triangle
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        NumRefShapeFunctions(), refcoords.cols());
    // Compute the barycentric coordinate functions
    const Eigen::RowVectorXd l1 = Eigen::RowVectorXd::Ones(refcoords.cols()) -
                                  refcoords.row(0) - refcoords.row(1);
    const Eigen::RowVectorXd l2 = refcoords.row(0);
    const Eigen::RowVectorXd l3 = refcoords.row(1);
    // Compute normalized coordinates for the scaled polynomials
    Eigen::RowVectorXd l121n(refcoords.cols());
    Eigen::RowVectorXd l122n(refcoords.cols());
    Eigen::RowVectorXd l232n(refcoords.cols());
    Eigen::RowVectorXd l233n(refcoords.cols());
    Eigen::RowVectorXd l313n(refcoords.cols());
    Eigen::RowVectorXd l311n(refcoords.cols());
    for (Eigen::Index i = 0; i < refcoords.cols(); ++i) {
      l121n[i] = l1[i] + l2[i] == 0 ? SCALAR(0) : l1[i] / (l1[i] + l2[i]);
      l122n[i] = l1[i] + l2[i] == 0 ? SCALAR(0) : l2[i] / (l1[i] + l2[i]);
      l232n[i] = l2[i] + l3[i] == 0 ? SCALAR(0) : l2[i] / (l2[i] + l3[i]);
      l233n[i] = l2[i] + l3[i] == 0 ? SCALAR(0) : l3[i] / (l2[i] + l3[i]);
      l313n[i] = l3[i] + l1[i] == 0 ? SCALAR(0) : l3[i] / (l3[i] + l1[i]);
      l311n[i] = l3[i] + l1[i] == 0 ? SCALAR(0) : l1[i] / (l3[i] + l1[i]);
    }
    // Get the basis functions associated with the vertices
    result.row(0) = l1.unaryExpr([&](double x) -> SCALAR { return x; });
    result.row(1) = l2.unaryExpr([&](double x) -> SCALAR { return x; });
    result.row(2) = l3.unaryExpr([&](double x) -> SCALAR { return x; });
    // Get the basis functions associated with the first edge
    for (unsigned i = 0; i < degree_ - 1; ++i) {
      if (rel_orient_[0] == lf::mesh::Orientation::positive) {
        // L_{i+2}(\lambda_2 ; \lambda_1+\lambda_2)
        result.row(3 + i) = ((l1 + l2)
                                 .unaryExpr([&](double x) -> SCALAR {
                                   return std::pow(x, i + 2);
                                 })
                                 .array() *
                             l122n.array().unaryExpr([&](double x) -> SCALAR {
                               return LegendrePoly<SCALAR>::integral(i + 2, x);
                             })).matrix();
      } else {
        // L_{i+2}(\lambda_1 ; \lambda_1+\lambda_2)
        result.row(degree_ + 1 - i) =
            ((l1 + l2)
                 .unaryExpr(
                     [&](double x) -> SCALAR { return std::pow(x, i + 2); })
                 .array() *
             l121n.array().unaryExpr([&](double x) -> SCALAR {
               return LegendrePoly<SCALAR>::integral(i + 2, x);
             })).matrix();
      }
    }
    // Get the basis functions associated with the second edge
    for (unsigned i = 0; i < degree_ - 1; ++i) {
      if (rel_orient_[1] == lf::mesh::Orientation::positive) {
        // L_{i+2}(\lambda_3 ; \lambda_2+\lambda_3)
        result.row(degree_ + 2 + i) =
            ((l2 + l3)
                 .unaryExpr(
                     [&](double x) -> SCALAR { return std::pow(x, i + 2); })
                 .array() *
             l233n.array().unaryExpr([&](double x) -> SCALAR {
               return LegendrePoly<SCALAR>::integral(i + 2, x);
             })).matrix();
      } else {
        // L_{i+2}(\lambda_2 ; \lambda_2+\lambda_3)
        result.row(2 * degree_ - i) =
            ((l2 + l3)
                 .unaryExpr(
                     [&](double x) -> SCALAR { return std::pow(x, i + 2); })
                 .array() *
             l232n.array().unaryExpr([&](double x) -> SCALAR {
               return LegendrePoly<SCALAR>::integral(i + 2, x);
             })).matrix();
      }
    }
    // Get the basis functions associated with the third edge
    for (unsigned i = 0; i < degree_ - 1; ++i) {
      if (rel_orient_[2] == lf::mesh::Orientation::positive) {
        // L_{i+2}(\lambda_1 ; \lambda_3+\lambda_1)
        result.row(2 * degree_ + 1 + i) =
            ((l3 + l1)
                 .unaryExpr(
                     [&](double x) -> SCALAR { return std::pow(x, i + 2); })
                 .array() *
             l311n.array().unaryExpr([&](double x) -> SCALAR {
               return LegendrePoly<SCALAR>::integral(i + 2, x);
             })).matrix();
      } else {
        // L_{i+2}(\lambda_3 ; \lambda_3+\lambda_1)
        result.row(3 * degree_ - 1 - i) =
            ((l3 + l1)
                 .unaryExpr(
                     [&](double x) -> SCALAR { return std::pow(x, i + 2); })
                 .array() *
             l313n.array().unaryExpr([&](double x) -> SCALAR {
               return LegendrePoly<SCALAR>::integral(i + 2, x);
             })).matrix();
      }
    }
    // Get the basis functions associated with the interior of the triangle
    if (degree_ > 2) {
      unsigned idx = 3 * degree_;
      // i is the degree of the edge function
      for (unsigned i = 0; i < degree_ - 2; ++i) {
        // j is the degree of the blending polynomial
        for (unsigned j = 0; j < degree_ - i - 2; ++j) {
          if (rel_orient_[0] == lf::mesh::Orientation::positive) {
            // L_{i+2}(\lambda_3 ; \lambda_2+\lambda_3) *
            // P_{j+1}^{2i+4}(\lambda_1)
            result.row(idx) =
                (result.row(3 + i).array() *
                 l3.array().unaryExpr([&](double x) -> SCALAR {
                   return JacobiPoly<SCALAR>::integral(j + 1, 2 * i + 4, x);
                 })).matrix();
          } else {
            // L_{i+2}(\lambda_2 ; \lambda_2+\lambda_3) *
            // P_{j+1}^{2i+4}(\lambda_1)
            result.row(idx) =
                (result.row(degree_ + 1 - i).array() *
                 l3.array().unaryExpr([&](double x) -> SCALAR {
                   return JacobiPoly<SCALAR>::integral(j + 1, 2 * i + 4, x);
                 })).matrix();
          }
          ++idx;
        }
      }
    }
    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        NumRefShapeFunctions(), 2 * refcoords.cols());
    // Compute the barycentric coordinate functions
    const Eigen::RowVectorXd l1 = Eigen::RowVectorXd::Ones(refcoords.cols()) -
                                  refcoords.row(0) - refcoords.row(1);
    const Eigen::RowVectorXd l2 = refcoords.row(0);
    const Eigen::RowVectorXd l3 = refcoords.row(1);
    // Comput the gradients of the barycentrc coordinate functions
    const Eigen::RowVectorXd l1_dx =
        Eigen::RowVectorXd::Constant(refcoords.cols(), -1);
    const Eigen::RowVectorXd l1_dy =
        Eigen::RowVectorXd::Constant(refcoords.cols(), -1);
    const Eigen::RowVectorXd l2_dx =
        Eigen::RowVectorXd::Constant(refcoords.cols(), 1);
    const Eigen::RowVectorXd l2_dy =
        Eigen::RowVectorXd::Constant(refcoords.cols(), 0);
    const Eigen::RowVectorXd l3_dx =
        Eigen::RowVectorXd::Constant(refcoords.cols(), 0);
    const Eigen::RowVectorXd l3_dy =
        Eigen::RowVectorXd::Constant(refcoords.cols(), 1);
    // Iterate over all refcoords
    for (Eigen::Index i = 0; i < refcoords.cols(); ++i) {
      // Compute normalized coordinates for edge 1
      const SCALAR l1p2 = l1[i] + l2[i];
      const SCALAR l1p2_dx = l1_dx[i] + l2_dx[i];
      const SCALAR l1p2_dy = l1_dy[i] + l2_dy[i];
      const SCALAR l121n = l1p2 == 0 ? SCALAR(0) : (l1[i] / l1p2);
      const SCALAR l121n_dx =
          l1p2 == 0 ? SCALAR(0)
                    : ((l1_dx[i] * l1p2 - l1[i] * l1p2_dx) / (l1p2 * l1p2));
      const SCALAR l121n_dy =
          l1p2 == 0 ? SCALAR(0)
                    : ((l1_dy[i] * l1p2 - l1[i] * l1p2_dy) / (l1p2 * l1p2));
      const SCALAR l122n = l1p2 == 0 ? SCALAR(0) : (l2[i] / l1p2);
      const SCALAR l122n_dx =
          l1p2 == 0 ? SCALAR(0)
                    : ((l2_dx[i] * l1p2 - l2[i] * l1p2_dx) / (l1p2 * l1p2));
      const SCALAR l122n_dy =
          l1p2 == 0 ? SCALAR(0)
                    : ((l2_dy[i] * l1p2 - l2[i] * l1p2_dy) / (l1p2 * l1p2));
      // Compute normalized coordinates for edge 2
      const SCALAR l2p3 = l2[i] + l3[i];
      const SCALAR l2p3_dx = l2_dx[i] + l3_dx[i];
      const SCALAR l2p3_dy = l2_dy[i] + l3_dy[i];
      const SCALAR l232n = l2p3 == 0 ? SCALAR(0) : (l2[i] / l2p3);
      const SCALAR l232n_dx =
          l2p3 == 0 ? SCALAR(0)
                    : ((l2_dx[i] * l2p3 - l2[i] * l2p3_dx) / (l2p3 * l2p3));
      const SCALAR l232n_dy =
          l2p3 == 0 ? SCALAR(0)
                    : ((l2_dy[i] * l2p3 - l2[i] * l2p3_dy) / (l2p3 * l2p3));
      const SCALAR l233n = l2p3 == 0 ? SCALAR(0) : (l3[i] / l2p3);
      const SCALAR l233n_dx =
          l2p3 == 0 ? SCALAR(0)
                    : ((l3_dx[i] * l2p3 - l3[i] * l2p3_dx) / (l2p3 * l2p3));
      const SCALAR l233n_dy =
          l2p3 == 0 ? SCALAR(0)
                    : ((l3_dy[i] * l2p3 - l3[i] * l2p3_dy) / (l2p3 * l2p3));
      // Compute normalized coordinates for edge 3
      const SCALAR l3p1 = l3[i] + l1[i];
      const SCALAR l3p1_dx = l3_dx[i] + l1_dx[i];
      const SCALAR l3p1_dy = l3_dy[i] + l1_dy[i];
      const SCALAR l313n = l3p1 == 0 ? SCALAR(0) : (l3[i] / l3p1);
      const SCALAR l313n_dx =
          l3p1 == 0 ? SCALAR(0)
                    : ((l3_dx[i] * l3p1 - l3[i] * l3p1_dx) / (l3p1 * l3p1));
      const SCALAR l313n_dy =
          l3p1 == 0 ? SCALAR(0)
                    : ((l3_dy[i] * l3p1 - l3[i] * l3p1_dy) / (l3p1 * l3p1));
      const SCALAR l311n = l3p1 == 0 ? SCALAR(0) : (l1[i] / l3p1);
      const SCALAR l311n_dx =
          l3p1 == 0 ? SCALAR(0)
                    : ((l1_dx[i] * l3p1 - l1[i] * l3p1_dx) / (l3p1 * l3p1));
      const SCALAR l311n_dy =
          l3p1 == 0 ? SCALAR(0)
                    : ((l1_dy[i] * l3p1 - l1[i] * l3p1_dy) / (l3p1 * l3p1));
      // Get the gradient of the basis functions associated with the vertices
      result(0, 2 * i + 0) = l1_dx[i];
      result(0, 2 * i + 1) = l1_dy[i];
      result(1, 2 * i + 0) = l2_dx[i];
      result(1, 2 * i + 1) = l2_dy[i];
      result(2, 2 * i + 0) = l3_dx[i];
      result(2, 2 * i + 1) = l3_dy[i];
      // Get the gradient of the basis functions associated with the first edge
      for (int j = 0; j < degree_ - 1; ++j) {
        if (rel_orient_[0] == lf::mesh::Orientation::positive) {
          const SCALAR leg2inte = LegendrePoly<SCALAR>::integral(j + 2, l122n);
          const SCALAR leg2eval = LegendrePoly<SCALAR>::eval(j + 1, l122n);
          result(3 + j, 2 * i + 0) =
              l1p2_dx * (j + 2) * std::pow(l1p2, j + 1) * leg2inte +
              std::pow(l1p2, j + 2) * l122n_dx * leg2eval;
          result(3 + j, 2 * i + 1) =
              l1p2_dy * (j + 2) * std::pow(l1p2, j + 1) * leg2inte +
              std::pow(l1p2, j + 2) * l122n_dy * leg2eval;
        } else {
          const SCALAR leg1inte = LegendrePoly<SCALAR>::integral(j + 2, l121n);
          const SCALAR leg1eval = LegendrePoly<SCALAR>::eval(j + 1, l121n);
          result(degree_ + 1 - j, 2 * i + 0) =
              l1p2_dx * (j + 2) * std::pow(l1p2, j + 1) * leg1inte +
              std::pow(l1p2, j + 2) * l121n_dx * leg1eval;
          result(degree_ + 1 - j, 2 * i + 1) =
              l1p2_dy * (j + 2) * std::pow(l1p2, j + 1) * leg1inte +
              std::pow(l1p2, j + 2) * l121n_dy * leg1eval;
        }
      }
      // Get the gradient of the basis functions associated with the second edge
      for (int j = 0; j < degree_ - 1; ++j) {
        if (rel_orient_[1] == lf::mesh::Orientation::positive) {
          const SCALAR leg3inte = LegendrePoly<SCALAR>::integral(j + 2, l233n);
          const SCALAR leg3eval = LegendrePoly<SCALAR>::eval(j + 1, l233n);
          result(2 + degree_ + j, 2 * i + 0) =
              l2p3_dx * (j + 2) * std::pow(l2p3, j + 1) * leg3inte +
              std::pow(l2p3, j + 2) * l233n_dx * leg3eval;
          result(2 + degree_ + j, 2 * i + 1) =
              l2p3_dy * (j + 2) * std::pow(l2p3, j + 1) * leg3inte +
              std::pow(l2p3, j + 2) * l233n_dy * leg3eval;
        } else {
          const SCALAR leg2inte = LegendrePoly<SCALAR>::integral(j + 2, l232n);
          const SCALAR leg2eval = LegendrePoly<SCALAR>::eval(j + 1, l232n);
          result(2 * degree_ - j, 2 * i + 0) =
              l2p3_dx * (j + 2) * std::pow(l2p3, j + 1) * leg2inte +
              std::pow(l2p3, j + 2) * l232n_dx * leg2eval;
          result(2 * degree_ - j, 2 * i + 1) =
              l2p3_dy * (j + 2) * std::pow(l2p3, j + 1) * leg2inte +
              std::pow(l2p3, j + 2) * l232n_dy * leg2eval;
        }
      }
      // Get the gradient of the basis functions associated with the third edge
      for (int j = 0; j < degree_ - 1; ++j) {
        if (rel_orient_[2] == lf::mesh::Orientation::positive) {
          const SCALAR leg1inte = LegendrePoly<SCALAR>::integral(j + 2, l311n);
          const SCALAR leg1eval = LegendrePoly<SCALAR>::eval(j + 1, l311n);
          result(1 + 2 * degree_ + j, 2 * i + 0) =
              l3p1_dx * (j + 2) * std::pow(l3p1, j + 1) * leg1inte +
              std::pow(l3p1, j + 2) * l311n_dx * leg1eval;
          result(1 + 2 * degree_ + j, 2 * i + 1) =
              l3p1_dy * (j + 2) * std::pow(l3p1, j + 1) * leg1inte +
              std::pow(l3p1, j + 2) * l311n_dy * leg1eval;
        } else {
          const SCALAR leg3inte = LegendrePoly<SCALAR>::integral(j + 2, l313n);
          const SCALAR leg3eval = LegendrePoly<SCALAR>::eval(j + 1, l313n);
          result(3 * degree_ - 1 - j, 2 * i + 0) =
              l3p1_dx * (j + 2) * std::pow(l3p1, j + 1) * leg3inte +
              std::pow(l3p1, j + 2) * l313n_dx * leg3eval;
          result(3 * degree_ - 1 - j, 2 * i + 1) =
              l3p1_dy * (j + 2) * std::pow(l3p1, j + 1) * leg3inte +
              std::pow(l3p1, j + 2) * l313n_dy * leg3eval;
        }
      }
      // Get the gradient of the basis functions associated with the interior of
      // the triangle
      if (degree_ > 2) {
        unsigned idx = 3 * degree_;
        for (unsigned j = 0; j < degree_ - 2; ++j) {
          SCALAR edge_eval;
          SCALAR edge_dx;
          SCALAR edge_dy;
          if (rel_orient_[0] == lf::mesh::Orientation::positive) {
            edge_eval = std::pow(l1p2, j + 2) *
                        LegendrePoly<SCALAR>::integral(j + 2, l122n);
            edge_dx = result(3 + j, 2 * i + 0);
            edge_dy = result(3 + j, 2 * i + 1);
          } else {
            edge_eval = std::pow(l1p2, j + 2) *
                        LegendrePoly<SCALAR>::integral(j + 2, l121n);
            edge_dx = result(degree_ + 1 - j, 2 * i + 0);
            edge_dy = result(degree_ + 1 - j, 2 * i + 1);
          }
          for (unsigned k = 0; k < degree_ - j - 2; ++k) {
            SCALAR jackinte =
                JacobiPoly<SCALAR>::integral(k + 1, 2 * j + 4, l3[i]);
            SCALAR jackeval = JacobiPoly<SCALAR>::eval(k, 2 * j + 4, l3[i]);
            result(idx, 2 * i + 0) =
                jackinte * edge_dx + edge_eval * jackeval * l3_dx[i];
            result(idx, 2 * i + 1) =
                jackinte * edge_dy + edge_eval * jackeval * l3_dy[i];
            ++idx;
          }
        }
      }
    }
    return result;
  }

  /**
   * @brief Evaluation nodes are the vertices, the points of a Gauss
   * quadrature rule for each edge and the points of a quadrature rule
   * on the interior of the triangle.
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    const lf::base::size_type Ns = qr_dual_segment_.NumPoints();
    const lf::base::size_type Nt = qr_dual_tria_.NumPoints();
    Eigen::MatrixXd nodes(2, 3 + 3 * Ns + Nt);
    // Add the vertices
    nodes(0, 0) = 0;
    nodes(1, 0) = 0;
    nodes(0, 1) = 1;
    nodes(1, 1) = 0;
    nodes(0, 2) = 0;
    nodes(1, 2) = 1;
    // Add the quadrature points on the first edge
    if (rel_orient_[0] == lf::mesh::Orientation::positive) {
      nodes.block(0, 3, 1, Ns) = qr_dual_segment_.Points();
    } else {
      nodes.block(0, 3, 1, Ns) =
          Eigen::RowVectorXd::Ones(Ns) - qr_dual_segment_.Points();
    }
    nodes.block(1, 3, 1, Ns).setZero();
    // Add the quadrature points on the second edge
    if (rel_orient_[1] == lf::mesh::Orientation::positive) {
      nodes.block(0, 3 + Ns, 1, Ns) =
          Eigen::RowVectorXd::Ones(Ns) - qr_dual_segment_.Points();
      nodes.block(1, 3 + Ns, 1, Ns) = qr_dual_segment_.Points();
    } else {
      nodes.block(0, 3 + Ns, 1, Ns) = qr_dual_segment_.Points();
      nodes.block(1, 3 + Ns, 1, Ns) =
          Eigen::RowVectorXd::Ones(Ns) - qr_dual_segment_.Points();
    }
    // Add the quadrature points on the third edge
    nodes.block(0, 3 + 2 * Ns, 1, Ns).setZero();
    if (rel_orient_[2] == lf::mesh::Orientation::positive) {
      nodes.block(1, 3 + 2 * Ns, 1, Ns) =
          Eigen::RowVectorXd::Ones(Ns) - qr_dual_segment_.Points();
    } else {
      nodes.block(1, 3 + 2 * Ns, 1, Ns) = qr_dual_segment_.Points();
    }
    // Add the quadrature points for the interior
    nodes.block(0, 3 + 3 * Ns, 2, Nt) = qr_dual_tria_.Points();
    return nodes;
  }

  /**
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    const lf::base::size_type Ns = qr_dual_segment_.NumPoints();
    const lf::base::size_type Nt = qr_dual_tria_.NumPoints();
    return 3 + 3 * Ns + Nt;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const override {
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(NumRefShapeFunctions());
    // Compute the basis function coefficients of the vertex basis functions
    dofs.segment(0, 3) = NodalValuesToVertexDofs(nodevals);
    // Compute the basis function coefficients on the edges
    dofs.segment(3, 3 * (Degree() - 1)) = NodalValuesToEdgeDofs(nodevals);
    // Compute the basis function coefficients for the face bubbles
    dofs.segment(3 * Degree(), NumRefShapeFunctions(0)) =
        NodalValuesToFaceDofs(nodevals);
    // We need to orthogonalize the face bubble dual basis w.r.t. the
    // dual basis on the vertices and edges
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
        basis_functions = EvalReferenceShapeFunctions(EvaluationNodes());
    const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> boundary_function =
        dofs.segment(0, 3 * Degree()) *
        basis_functions.block(0, 0, 3 * Degree(), basis_functions.cols());
    const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> boundary_face_dofs =
        NodalValuesToFaceDofs(boundary_function);
    dofs.segment(3 * Degree(), NumRefShapeFunctions(0)) -= boundary_face_dofs;
    return dofs;
  }

 private:
  unsigned degree_;
  lf::quad::QuadRule qr_dual_segment_;
  lf::quad::QuadRule qr_dual_tria_;
  nonstd::span<const lf::mesh::Orientation> rel_orient_;

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>
  NodalValuesToVertexDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const {
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(3);
    dofs[0] = nodevals[0];
    dofs[1] = nodevals[1];
    dofs[2] = nodevals[2];
    return dofs;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToEdgeDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const {
    const lf::base::size_type Ns = qr_dual_segment_.NumPoints();
    const lf::base::size_type Nt = qr_dual_tria_.NumPoints();
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(3 * (Degree() - 1));
    // Compute the basis function coefficients on the edges
    // by applying the dual basis of the segment
    for (lf::base::size_type i = 2; i < Degree() + 1; ++i) {
      const SCALAR P0 = LegendrePoly<SCALAR>::eval(i - 1, 0);
      const SCALAR P1 = LegendrePoly<SCALAR>::eval(i - 1, 1);
      Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
          qr_dual_segment_.Points().unaryExpr([&](double x) -> SCALAR {
            return LegendrePoly<SCALAR>::derivative(i - 2, x);
          });
      // Compute the basis function coefficients for the first edge
      const SCALAR integ1 = (qr_dual_segment_.Weights().transpose().array() *
                             psidd.array() * nodevals.segment(3, Ns).array())
                                .sum();
      if (rel_orient_[0] == lf::mesh::Orientation::positive) {
        dofs[1 + i - 3] =
            (P1 * nodevals[1] - P0 * nodevals[0] - integ1) * (2 * i - 1);
      } else {
        dofs[4 - i + (Degree() - 1) - 3] =
            (P1 * nodevals[0] - P0 * nodevals[1] - integ1) * (2 * i - 1);
      }
      // Compute the basis function coefficients for the second edge
      const SCALAR integ2 =
          (qr_dual_segment_.Weights().transpose().array() * psidd.array() *
           nodevals.segment(3 + Ns, Ns).array())
              .sum();
      if (rel_orient_[1] == lf::mesh::Orientation::positive) {
        dofs[1 + i + (Degree() - 1) - 3] =
            (P1 * nodevals[2] - P0 * nodevals[1] - integ2) * (2 * i - 1);
      } else {
        dofs[4 - i + 2 * (Degree() - 1) - 3] =
            (P1 * nodevals[1] - P0 * nodevals[2] - integ2) * (2 * i - 1);
      }
      // Compute the basis function coefficients for the third edge
      const SCALAR integ3 =
          (qr_dual_segment_.Weights().transpose().array() * psidd.array() *
           nodevals.segment(3 + 2 * Ns, Ns).array())
              .sum();
      if (rel_orient_[2] == lf::mesh::Orientation::positive) {
        dofs[1 + i + 2 * (Degree() - 1) - 3] =
            (P1 * nodevals[0] - P0 * nodevals[2] - integ3) * (2 * i - 1);
      } else {
        dofs[4 - i + 3 * (Degree() - 1) - 3] =
            (P1 * nodevals[2] - P0 * nodevals[0] - integ3) * (2 * i - 1);
      }
    }
    return dofs;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToFaceDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const {
    const lf::base::size_type Ns = qr_dual_segment_.NumPoints();
    const lf::base::size_type Nt = qr_dual_tria_.NumPoints();
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(NumRefShapeFunctions(0));
    if (Degree() > 2) {
      unsigned idx = 0;
      const Eigen::Matrix<double, 1, Eigen::Dynamic> xnorm =
          (qr_dual_tria_.Points().row(0).array() /
           qr_dual_tria_.Points().row(1).array().unaryExpr([&](double y) {
             return 1 - y;
           })).matrix();
      const Eigen::Matrix<double, 1, Eigen::Dynamic> xnorm_adj =
          rel_orient_[0] == lf::mesh::Orientation::positive
              ? xnorm
              : Eigen::RowVectorXd::Ones(Nt) - xnorm;
      const Eigen::Matrix<double, 1, Eigen::Dynamic> y =
          qr_dual_tria_.Points().row(1);
      // i is the degree of the edge function
      for (unsigned i = 0; i < Degree() - 2; ++i) {
        const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> ypow =
            y.array()
                .unaryExpr(
                    [&](double y) -> SCALAR { return std::pow(1 - y, i); })
                .matrix();
        const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> ypowp1 =
            y.array()
                .unaryExpr(
                    [&](double y) -> SCALAR { return std::pow(1 - y, i + 1); })
                .matrix();
        const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
            xnorm_adj.unaryExpr([&](double x) -> SCALAR {
              return LegendrePoly<SCALAR>::derivative(i, x);
            });
        // j is the degree of the blending Jacobi polynomial
        for (unsigned j = 0; j < Degree() - i - 2; ++j) {
          const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> jacdd =
              qr_dual_tria_.Points().row(1).unaryExpr([&](double y) -> SCALAR {
                return j == 0 ? SCALAR(0)
                              : JacobiPoly<SCALAR>::derivative(j - 1, 2 * i + 4,
                                                               y);
              });
          const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> jacd =
              qr_dual_tria_.Points().row(1).unaryExpr([&](double y) -> SCALAR {
                return JacobiPoly<SCALAR>::eval(j, 2 * i + 4, y);
              });
          dofs[idx] =
              (qr_dual_tria_.Weights().transpose().array() *
               nodevals.block(0, 3 + 3 * Ns, 1, Nt).array() * psidd.array() *
               (ypowp1.array() * jacdd.array() -
                (2 * i + 4) * ypow.array() * jacd.array()))
                  .sum();
          dofs[idx] *=
              2 * i + 3;  // Normalization factor for the Legendre polynomial
          dofs[idx] *= 2 * j + 2 * i +
                       5;  // Normalization factor for the Jacobi Polynomial
          ++idx;
        }
      }
    }
    return dofs;
  }
};

/**
 * @headerfile lf/fe/fe.h
 * @brief Hierarchic Finite Elements of arbitrary degree on quadrilaterals
 *
 * The basis functions on the quadrilateral has a tensor product structure and
 * can thus be represented by products of basis functions on segments.
 * The basis functions are therefore given by
 * \f[
 *  \widehat{b^{\square}}^{ij}(\widehat{x}, \widehat{y}) =
 * \widehat{b}^i(\widehat{x})\widehat{b}^j(\widehat{y}) \f] where
 * \f$\widehat{b}^i(\widehat{x})\f$ and \f$\widehat{b}^j(\widehat{y})\f$ are the
 * basis functions on a segment.
 *
 * The dual basis is therefore also quite simple, as we can recycle the one from
 * the segments by first applying the dual basis along the \f$x\f$-axis and then
 * apply the dual basis to the resulting 1d function.
 *
 * @note The Shape Functions are taken from the following paper:
 * https://arxiv.org/pdf/1504.03025.pdf
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeHierarchicQuad final : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeHierarchicQuad(const FeHierarchicQuad &) = default;
  FeHierarchicQuad(FeHierarchicQuad &&) noexcept = default;
  FeHierarchicQuad &operator=(const FeHierarchicQuad &) = default;
  FeHierarchicQuad &operator=(FeHierarchicQuad &&) noexcept = default;
  ~FeHierarchicQuad() override = default;

  FeHierarchicQuad(unsigned degree,
                   nonstd::span<const lf::mesh::Orientation> rel_orient)
      : ScalarReferenceFiniteElement<SCALAR>(),
        degree_(degree),
        qr_dual_(lf::quad::make_QuadRule(lf::base::RefEl::kSegment(),
                                         2 * (degree - 1))),
        rel_orient_(rel_orient),
        fe1d_(degree, rel_orient) {}

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kQuad();
  }

  [[nodiscard]] unsigned Degree() const override { return degree_; }

  /**
   * @brief The local shape functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
    return (degree_ + 1) * (degree_ + 1);
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and (p-1)^2 shape functions on the quadrilateral
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim) const override {
    switch (codim) {
      case 0:
        return (degree_ - 1) * (degree_ - 1);
      case 1:
        return degree_ - 1;
      case 2:
        return 1;
      default:
        LF_ASSERT_MSG(false, "Illegal codim " << codim);
        // Silence compiler warnings
        return 0;
    }
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and (p-1)^2 shape functions on the quadrilateral
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        NumRefShapeFunctions(), refcoords.cols());
    // Compute the 1D shape functions at the x and y coordinates
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_x =
        fe1d_.EvalReferenceShapeFunctions(refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_y =
        fe1d_.EvalReferenceShapeFunctions(refcoords.row(1));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_x =
        fe1d_.EvalReferenceShapeFunctions(
            Eigen::RowVectorXd::Constant(refcoords.cols(), 1) -
            refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_y =
        fe1d_.EvalReferenceShapeFunctions(
            Eigen::RowVectorXd::Constant(refcoords.cols(), 1) -
            refcoords.row(1));
    // Get the basis functions associated with the vertices
    result.row(0) = (sf1d_x.row(0).array() * sf1d_y.row(0).array()).matrix();
    result.row(1) = (sf1d_x.row(1).array() * sf1d_y.row(0).array()).matrix();
    result.row(2) = (sf1d_x.row(1).array() * sf1d_y.row(1).array()).matrix();
    result.row(3) = (sf1d_x.row(0).array() * sf1d_y.row(1).array()).matrix();
    // Get the basis functions associated with the first edge
    // as a tensor product of 1D basis functions
    for (int i = 0; i < degree_ - 1; ++i) {
      if (rel_orient_[0] == lf::mesh::Orientation::positive) {
        result.row(4 + i) =
            (sf1d_x.row(2 + i).array() * sf1d_y.row(0).array()).matrix();
      } else {
        result.row(2 + degree_ - i) =
            (sf1df_x.row(2 + i).array() * sf1d_y.row(0).array()).matrix();
      }
    }
    // Get the basis functions associated with the second edge
    // as a tensor product of 1D basis functions
    for (int i = 0; i < degree_ - 1; ++i) {
      if (rel_orient_[1] == lf::mesh::Orientation::positive) {
        result.row(3 + degree_ + i) =
            (sf1d_x.row(1).array() * sf1d_y.row(2 + i).array()).matrix();
      } else {
        result.row(1 + 2 * degree_ - i) =
            (sf1d_x.row(1).array() * sf1df_y.row(2 + i).array()).matrix();
      }
    }
    // Get the basis functions associated with the third edge
    // as a tensor product of 1D basis functions
    for (int i = 0; i < degree_ - 1; ++i) {
      if (rel_orient_[2] == lf::mesh::Orientation::positive) {
        result.row(2 + 2 * degree_ + i) =
            (sf1df_x.row(2 + i).array() * sf1d_y.row(1).array()).matrix();
      } else {
        result.row(3 * degree_ - i) =
            (sf1d_x.row(2 + i).array() * sf1d_y.row(1).array()).matrix();
      }
    }
    // Get the basis functions associated with the fourth edge
    // as a tensor product of 1D basis functions
    for (int i = 0; i < degree_ - 1; ++i) {
      if (rel_orient_[3] == lf::mesh::Orientation::positive) {
        result.row(1 + 3 * degree_ + i) =
            (sf1d_x.row(0).array() * sf1df_y.row(2 + i).array()).matrix();
      } else {
        result.row(4 * degree_ - 1 - i) =
            (sf1d_x.row(0).array() * sf1d_y.row(2 + i).array()).matrix();
      }
    }
    // Get the basis functions associated with the interior of the quad
    // as a tensor product of 1D basis functions
    for (int i = 0; i < degree_ - 1; ++i) {
      for (int j = 0; j < degree_ - 1; ++j) {
        result.row(4 * degree_ + (degree_ - 1) * i + j) =
            (sf1d_x.row(j + 2).array() * sf1d_y.row(i + 2).array()).matrix();
      }
    }
    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        NumRefShapeFunctions(), 2 * refcoords.cols());
    // Compute the gradient of the 1D shape functions at the x and y coordinates
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_x =
        fe1d_.EvalReferenceShapeFunctions(refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_y =
        fe1d_.EvalReferenceShapeFunctions(refcoords.row(1));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_dx =
        fe1d_.GradientsReferenceShapeFunctions(refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1d_dy =
        fe1d_.GradientsReferenceShapeFunctions(refcoords.row(1));
    // Compute the gradient of the flipped 1D shape functions at the x and y
    // coordinates
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_x =
        fe1d_.EvalReferenceShapeFunctions(
            Eigen::RowVectorXd::Constant(refcoords.cols(), 1) -
            refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_y =
        fe1d_.EvalReferenceShapeFunctions(
            Eigen::RowVectorXd::Constant(refcoords.cols(), 1) -
            refcoords.row(1));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_dx =
        fe1d_.GradientsReferenceShapeFunctions(
            Eigen::RowVectorXd::Constant(refcoords.cols(), 1) -
            refcoords.row(0));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> sf1df_dy =
        fe1d_.GradientsReferenceShapeFunctions(
            Eigen::RowVectorXd::Constant(refcoords.cols(), 1) -
            refcoords.row(1));
    for (Eigen::Index i = 0; i < refcoords.cols(); ++i) {
      // Get the gradient of the basis functions associated with the vertices
      result(0, 2 * i + 0) = sf1d_dx(0, i) * sf1d_y(0, i);
      result(0, 2 * i + 1) = sf1d_x(0, i) * sf1d_dy(0, i);
      result(1, 2 * i + 0) = sf1d_dx(1, i) * sf1d_y(0, i);
      result(1, 2 * i + 1) = sf1d_x(1, i) * sf1d_dy(0, i);
      result(2, 2 * i + 0) = sf1d_dx(1, i) * sf1d_y(1, i);
      result(2, 2 * i + 1) = sf1d_x(1, i) * sf1d_dy(1, i);
      result(3, 2 * i + 0) = sf1d_dx(0, i) * sf1d_y(1, i);
      result(3, 2 * i + 1) = sf1d_x(0, i) * sf1d_dy(1, i);
      // Get the basis functions associated with the first edge
      for (int j = 0; j < degree_ - 1; ++j) {
        if (rel_orient_[0] == lf::mesh::Orientation::positive) {
          result(4 + j, 2 * i + 0) = sf1d_dx(2 + j, i) * sf1d_y(0, i);
          result(4 + j, 2 * i + 1) = sf1d_x(2 + j, i) * sf1d_dy(0, i);
        } else {
          result(2 + degree_ - j, 2 * i + 0) =
              -sf1df_dx(2 + j, i) * sf1d_y(0, i);
          result(2 + degree_ - j, 2 * i + 1) =
              sf1df_x(2 + j, i) * sf1d_dy(0, i);
        }
      }
      // Get the basis functions associated with the second edge
      for (int j = 0; j < degree_ - 1; ++j) {
        if (rel_orient_[1] == lf::mesh::Orientation::positive) {
          result(3 + degree_ + j, 2 * i + 0) = sf1d_dx(1, i) * sf1d_y(2 + j, i);
          result(3 + degree_ + j, 2 * i + 1) = sf1d_x(1, i) * sf1d_dy(2 + j, i);
        } else {
          result(1 + 2 * degree_ - j, 2 * i + 0) =
              sf1d_dx(1, i) * sf1df_y(2 + j, i);
          result(1 + 2 * degree_ - j, 2 * i + 1) =
              sf1d_x(1, i) * -sf1df_dy(2 + j, i);
        }
      }
      // Get the basis functions associated with the third edge
      for (int j = 0; j < degree_ - 1; ++j) {
        if (rel_orient_[2] == lf::mesh::Orientation::positive) {
          result(2 + 2 * degree_ + j, 2 * i + 0) =
              -sf1df_dx(2 + j, i) * sf1d_y(1, i);
          result(2 + 2 * degree_ + j, 2 * i + 1) =
              sf1df_x(2 + j, i) * sf1d_dy(1, i);
        } else {
          result(3 * degree_ - j, 2 * i + 0) = sf1d_dx(2 + j, i) * sf1d_y(1, i);
          result(3 * degree_ - j, 2 * i + 1) = sf1d_x(2 + j, i) * sf1d_dy(1, i);
        }
      }
      // Get the basis functions associated with the fourth edge
      for (int j = 0; j < degree_ - 1; ++j) {
        if (rel_orient_[3] == lf::mesh::Orientation::positive) {
          result(1 + 3 * degree_ + j, 2 * i + 0) =
              sf1d_dx(0, i) * sf1df_y(2 + j, i);
          result(1 + 3 * degree_ + j, 2 * i + 1) =
              sf1d_x(0, i) * -sf1df_dy(2 + j, i);
        } else {
          result(4 * degree_ - 1 - j, 2 * i + 0) =
              sf1d_dx(0, i) * sf1d_y(2 + j, i);
          result(4 * degree_ - 1 - j, 2 * i + 1) =
              sf1d_x(0, i) * sf1d_dy(2 + j, i);
        }
      }
      // Get the basis functions associated with the interior of the quad
      for (int j = 0; j < degree_ - 1; ++j) {
        for (int k = 0; k < degree_ - 1; ++k) {
          result(4 * degree_ + (degree_ - 1) * j + k, 2 * i + 0) =
              sf1d_dx(k + 2, i) * sf1d_y(j + 2, i);
          result(4 * degree_ + (degree_ - 1) * j + k, 2 * i + 1) =
              sf1d_x(k + 2, i) * sf1d_dy(j + 2, i);
        }
      }
    }
    return result;
  }

  /**
   * @brief Evaluation nodes are the vertices, the points of
   * a quadrature rule on each edge and the points of a quadrature
   * rule on the interior of the quadrilateral.
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    const lf::base::size_type N = qr_dual_.NumPoints();
    Eigen::MatrixXd nodes(2, 4 + 4 * N + N * N);
    // Add the vertices
    nodes(0, 0) = 0;
    nodes(1, 0) = 0;
    nodes(0, 1) = 1;
    nodes(1, 1) = 0;
    nodes(0, 2) = 1;
    nodes(1, 2) = 1;
    nodes(0, 3) = 0;
    nodes(1, 3) = 1;
    // Add the quadrature points on the first edge
    if (rel_orient_[0] == lf::mesh::Orientation::positive) {
      nodes.block(0, 4, 1, N) = qr_dual_.Points();
    } else {
      nodes.block(0, 4, 1, N) = Eigen::RowVectorXd::Ones(N) - qr_dual_.Points();
    }
    nodes.block(1, 4, 1, N).setZero();
    // Add the quadrature points on the second edge
    nodes.block(0, 4 + N, 1, N).setOnes();
    if (rel_orient_[1] == lf::mesh::Orientation::positive) {
      nodes.block(1, 4 + N, 1, N) = qr_dual_.Points();
    } else {
      nodes.block(1, 4 + N, 1, N) =
          Eigen::RowVectorXd::Ones(N) - qr_dual_.Points();
    }
    // Add the quadrature points on the third edge
    if (rel_orient_[2] == lf::mesh::Orientation::positive) {
      nodes.block(0, 4 + 2 * N, 1, N) =
          Eigen::RowVectorXd::Ones(N) - qr_dual_.Points();
    } else {
      nodes.block(0, 4 + 2 * N, 1, N) = qr_dual_.Points();
    }
    nodes.block(1, 4 + 2 * N, 1, N).setOnes();
    // Add the quadrature points on the fourth edge
    nodes.block(0, 4 + 3 * N, 1, N).setZero();
    if (rel_orient_[3] == lf::mesh::Orientation::positive) {
      nodes.block(1, 4 + 3 * N, 1, N) =
          Eigen::RowVectorXd::Ones(N) - qr_dual_.Points();
    } else {
      nodes.block(1, 4 + 3 * N, 1, N) = qr_dual_.Points();
    }
    // Add the quadrature points for the face
    for (lf::base::size_type i = 0; i < N; ++i) {
      nodes.block(0, 4 + (4 + i) * N, 1, N) = qr_dual_.Points();
      nodes.block(1, 4 + (4 + i) * N, 1, N)
          .setConstant(qr_dual_.Points()(0, i));
    }
    return nodes;
  }

  /**
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    const lf::base::size_type N = qr_dual_.NumPoints();
    return (N + 2) * (N + 2);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const override {
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(NumRefShapeFunctions());
    const lf::base::size_type N = qr_dual_.NumPoints();
    // Compute the basis function coefficients of the vertex basis functions
    dofs[0] = nodevals[0];
    dofs[1] = nodevals[1];
    dofs[2] = nodevals[2];
    dofs[3] = nodevals[3];
    // Compute the basis function coefficients on the edges
    // by applying the dual basis of the segment
    for (lf::base::size_type i = 2; i < Degree() + 1; ++i) {
      const SCALAR P0 = LegendrePoly<SCALAR>::eval(i - 1, 0);
      const SCALAR P1 = LegendrePoly<SCALAR>::eval(i - 1, 1);
      Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
          qr_dual_.Points().unaryExpr([&](double x) -> SCALAR {
            return LegendrePoly<SCALAR>::derivative(i - 2, x);
          });
      // Compute the basis function coefficients for the first edge
      const SCALAR integ1 = (qr_dual_.Weights().transpose().array() *
                             psidd.array() * nodevals.segment(4, N).array())
                                .sum();
      if (rel_orient_[0] == lf::mesh::Orientation::positive) {
        dofs[2 + i] =
            (P1 * nodevals[1] - P0 * nodevals[0] - integ1) * (2 * i - 1);
      } else {
        dofs[5 - i + (Degree() - 1)] =
            (P1 * nodevals[0] - P0 * nodevals[1] - integ1) * (2 * i - 1);
      }
      // Compute the basis function coefficients for the second edge
      const SCALAR integ2 = (qr_dual_.Weights().transpose().array() *
                             psidd.array() * nodevals.segment(4 + N, N).array())
                                .sum();
      if (rel_orient_[1] == lf::mesh::Orientation::positive) {
        dofs[2 + i + (Degree() - 1)] =
            (P1 * nodevals[2] - P0 * nodevals[1] - integ2) * (2 * i - 1);
      } else {
        dofs[5 - i + 2 * (Degree() - 1)] =
            (P1 * nodevals[1] - P0 * nodevals[2] - integ2) * (2 * i - 1);
      }
      // Compute the basis function coefficients for the third edge
      const SCALAR integ3 =
          (qr_dual_.Weights().transpose().array() * psidd.array() *
           nodevals.segment(4 + 2 * N, N).array())
              .sum();
      if (rel_orient_[2] == lf::mesh::Orientation::positive) {
        dofs[2 + i + 2 * (Degree() - 1)] =
            (P1 * nodevals[3] - P0 * nodevals[2] - integ3) * (2 * i - 1);
      } else {
        dofs[5 - i + 3 * (Degree() - 1)] =
            (P1 * nodevals[2] - P0 * nodevals[3] - integ3) * (2 * i - 1);
      }
      // Compute the basis function coefficients for the fourth edge
      const SCALAR integ4 =
          (qr_dual_.Weights().transpose().array() * psidd.array() *
           nodevals.segment(4 + 3 * N, N).array())
              .sum();
      if (rel_orient_[3] == lf::mesh::Orientation::positive) {
        dofs[2 + i + 3 * (Degree() - 1)] =
            (P1 * nodevals[0] - P0 * nodevals[3] - integ4) * (2 * i - 1);
      } else {
        dofs[5 - i + 4 * (Degree() - 1)] =
            (P1 * nodevals[3] - P0 * nodevals[0] - integ4) * (2 * i - 1);
      }
    }
    // Compute the basis function coefficients for the face bubbles
    for (unsigned i = 0; i < Degree() - 1; ++i) {
      // Apply the 1d dual along the x-axis
      const SCALAR P0i = LegendrePoly<SCALAR>::eval(i + 1, 0);
      const SCALAR P1i = LegendrePoly<SCALAR>::eval(i + 1, 1);
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psiddi =
          qr_dual_.Points().unaryExpr([&](double x) -> SCALAR {
            return LegendrePoly<SCALAR>::derivative(i, x);
          });
      Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> partial(N);
      for (lf::base::size_type j = 0; j < N; ++j) {
        const SCALAR integ =
            (qr_dual_.Weights().transpose().array() * psiddi.array() *
             nodevals.segment(4 + (4 + j) * N, N).array())
                .sum();
        SCALAR n0;
        if (rel_orient_[3] == lf::mesh::Orientation::positive) {
          n0 = nodevals[3 + 4 * N - j];
        } else {
          n0 = nodevals[4 + 3 * N + j];
        }
        SCALAR n1;
        if (rel_orient_[1] == lf::mesh::Orientation::positive) {
          n1 = nodevals[4 + N + j];
        } else {
          n1 = nodevals[3 + 2 * N - j];
        }
        partial[j] = (P1i * n1 - P0i * n0 - integ) * (2 * i + 3);
      }
      SCALAR f0_integ;
      if (rel_orient_[0] == lf::mesh::Orientation::positive) {
        f0_integ = (qr_dual_.Weights().transpose().array() * psiddi.array() *
                    nodevals.segment(4, N).array())
                       .sum();
      } else {
        f0_integ = (qr_dual_.Weights().transpose().array() * psiddi.array() *
                    nodevals.segment(4, N).array().reverse())
                       .sum();
      }
      const SCALAR f0 =
          (P1i * nodevals[1] - P0i * nodevals[0] - f0_integ) * (2 * i + 3);
      SCALAR f1_integ;
      if (rel_orient_[2] == lf::mesh::Orientation::positive) {
        f1_integ = (qr_dual_.Weights().transpose().array() * psiddi.array() *
                    nodevals.segment(4 + 2 * N, N).array().reverse())
                       .sum();
      } else {
        f1_integ = (qr_dual_.Weights().transpose().array() * psiddi.array() *
                    nodevals.segment(4 + 2 * N, N).array())
                       .sum();
      }
      const SCALAR f1 =
          (P1i * nodevals[2] - P0i * nodevals[3] - f1_integ) * (2 * i + 3);
      for (unsigned j = 0; j < Degree() - 1; ++j) {
        // Apply the 1d dual along the y-axis
        const SCALAR P0j = LegendrePoly<SCALAR>::eval(j + 1, 0);
        const SCALAR P1j = LegendrePoly<SCALAR>::eval(j + 1, 1);
        const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psiddj =
            qr_dual_.Points().unaryExpr([&](double x) -> SCALAR {
              return LegendrePoly<SCALAR>::derivative(j, x);
            });
        const SCALAR integ = (qr_dual_.Weights().transpose().array() *
                              psiddj.array() * partial.array())
                                 .sum();
        dofs[4 * Degree() + j * (Degree() - 1) + i] =
            (P1j * f1 - P0j * f0 - integ) * (2 * j + 3);
      }
    }
    return dofs;
  }

 private:
  unsigned degree_;
  lf::quad::QuadRule qr_dual_;
  FeHierarchicSegment<SCALAR> fe1d_;
  nonstd::span<const lf::mesh::Orientation> rel_orient_;
};

}  // end namespace lf::fe

#endif  // LF_USCALFE_HP_FE_H_
