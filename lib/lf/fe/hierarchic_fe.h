#ifndef LF_USCALFE_HP_FE_H_
#define LF_USCALFE_HP_FE_H_

/**
 * @file
 * @brief Data structures representing HP finite elements
 * @author Tobias Rohner
 * @date May 2020
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>
#include <lf/quad/quad.h>

#include <cmath>
#include <memory>
#include <vector>

#include "fe_point.h"
#include "scalar_reference_finite_element.h"

namespace lf::fe {

// clang-format off
/**
 * @brief computes the `n`-th degree scaled Legendre Polynomial \f$ P_n(x;t)
 * \f$
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
// clang-format on
double Legendre(unsigned n, double x, double t = 1);

// clang-format off
/**
 * @brief computes the integral of the (n-1)-th degree scaled Legendre
 * Polynomial
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
// clang-format on
double ILegendre(unsigned n, double x, double t = 1);

// clang-format off
/**
 * @brief Computes \f$ \frac{\partial}{\partial x} L_n(x;t) \f$
 * @param n The degree of the integrated scaled Legendre polynomial
 * @param x The evaluation coordinate
 * @param t The scaling parameter
 *
 * The derivative is simply given by \f$ \frac{\partial}{\partial x} L_n(x;t) =
 * P_{n-1}(x;t) \f$
 */
// clang-format on
double ILegendreDx(unsigned n, double x, double t = 1);

// clang-format off
/**
 * @brief Computes \f$ \frac{\partial}{\partial t} L_n(x;t) \f$
 * @param n The degree of the integrated scaled Legendre polynomial
 * @param x The evaluation coordinate
 * @param t The scaling parameter
 *
 * The derivative is given by
 * \f[
 *	\begin{aligned}
 *	    \frac{\partial}{\partial t} L_1(x;t) &= 0 \\
 *	    \frac{\partial}{\partial t} L_n(x;t) &= -\frac{1}{2} \left(
 * P_{n-1}(x;t) + tP_{n-2}(x;t) \right) \end{aligned} \f]
 */
// clang-format on
double ILegendreDt(unsigned n, double x, double t = 1);

// clang-format off
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
 * (2x-t)\frac{\partial}{\partial x}P_{n-1}(x;t) \\ \end{aligned} \f]
 */
// clang-format on
double LegendreDx(unsigned n, double x, double t = 1);

// clang-format off
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
 * (\alpha + \beta + 2)x \right) \\ P_{n+1}^{(\alpha,\beta)}(x) &= \frac{1}{a_n}
 * \left( (b_n+c_nx)P_n^{(\alpha,\beta)}(x) - d_nP_{n-1}^{(\alpha,\beta)}(x)
 * \right) \end{aligned} \f] where \f[ \begin{aligned}
 *	    a_n &= 2(n+1)(n+\alpha+\beta+1)(2n+\alpha+\beta) \\
 *	    b_n &= (2n+\alpha+\beta+1)(\alpha^2-\beta^2) \\
 *	    c_n &= (2n+\alpha+\beta)(2n+\alpha+\beta+1)(2n+\alpha+\beta+2) \\
 *	    d_n &= 2(n+\alpha)(n+\beta)(2n+\alpha+\beta+2)
 *	\end{aligned}
 * \f]
 */
// clang-format on
double Jacobi(unsigned n, double alpha, double beta, double x);

// clang-format off
/**
 * @brief Computes the n-th degree shifted Jacobi polynomial for \f$ \beta = 0
 * \f$
 * @param n The degree of the polynomial
 * @param alpha The \f$ \alpha \f$ parameter of the Jacobi polynomial
 * @param x The evaluation coordinate
 */
// clang-format on
double Jacobi(unsigned n, double alpha, double x);

// clang-format off
/**
 * @brief Evaluate the integral of the (n-1)-th degree Jacobi Polynomial for \f$
 * \beta = 0 \f$
 * @param n The degree of the integrated polynomial
 * @param alpha The \f$ \alpha \f$ parameter of the Jacobi polynomial
 * @param x The evaluation coordinate
 *
 * The integral is evaluated using
 * \f[
 *	\begin{aligned}
 *	    L_1^\alpha(x) &= x \\
 *	    L_p^\alpha(x) &= a_pP_p^\alpha(x) + b_pP_{p-1}^\alpha(x) -
 * c_pP_{p-2}^\alpha(x) \end{aligned} \f] where the coefficients are defined as
 * \f[
 *	\begin{aligned}
 *	    a_p &= \frac{p+\alpha}{(2p+\alpha-1)(2p+\alpha)} \\
 *	    b_p &= \frac{\alpha}{(2p+\alpha-1)(2p+\alpha)} \\
 *	    c_p &= \frac{p-1}{(2p+\alpha-2)(2p+\alpha-1)}
 *	\end{aligned}
 * \f]
 */
// clang-format on
double IJacobi(unsigned n, double alpha, double x);

// clang-format off
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
// clang-format on
double IJacobiDx(unsigned n, double alpha, double x);

// clang-format off
/**
 * @brief Computes the derivative of the n-th degree Jacobi Polynomial for \f$
 * \beta = 0 \f$
 * @param n The degree of the differentiated polynomial
 * @param alpha The \f$ \alpha \f$ parameter of the Jacobi Polynomial
 * @param x The evaluation coordinate
 *
 * The derivative is evaluated using
 * \f[
 *	{P^{(\alpha,0)}_n}'(x) = \frac{\alpha+n+1}{2} P^{(\alpha+1,1)}_{n-1}(x)
 * \f]
 */
// clang-format on
double JacobiDx(unsigned n, double alpha, double x);

// clang-format off
/**
 * @headerfile lf/fe/fe.h
 * @brief Hierarchic Finite Elements of arbitrary degree on segments
 *
 * The shape functions associated with the vertices are given by
 * \f[
 *  \begin{align*}
 *	\widehat{b^\cdot}^1(x) &:= 1 - x \\
 *	\widehat{b^\cdot}^2(x) &:= x
 *  \end{align*}
 * \f]
 * and the interior basis functions associated with the segment itself
 * are given by the integrated shifted Legendre polynomials
 * \f[
 *  \widehat{b^-}^i(x) = L_i(x) = \int_0^{x}\!
 *  P_{i-1}(\xi) \,\mathrm{d}\xi \quad\mbox{ for }\quad i= 2, \cdots, p \f] where
 * \f$P_i : [0, 1] \to \mathbb{R}\f$ is the shifted Legendre polynomial of degree
 * \f$i\f$.
 *
 * To compute the basis function coefficients from point evaluations of a
 * function, we make use of the dual basis given by \f[ \lambda_i^-[f] =
 * \begin{cases}
 *	f(0) &\mbox{ for } i = 0 \\
 *	f(1) &\mbox{ for } i = 1 \\
 *	\frac{1}{2i - 1} \left [P_{i-1}(1)f(1) - P_{i-1}(0)f(0) - \int_0^1\!
 *	P_{i-1}'(x)f(x) \,\mathrm{d}x \right] &\mbox{ for } i \geq 2 \end{cases} \f]
 *
 * @attention Note that the local coordinate \f$\widehat{x}\f$ may be flipped by
 * applying an affine transformation \f$\widehat{x} \mapsto 1 - \widehat{x}\f$.
 * This must be done if the relative orientation of the edge is
 *`lf::mesh::Orientation::negative` in order to keep a global ordering of DOFs
 * on the cell interfaces.
 *
 * A complete description of the basis functions and dual basis can be found
 * <a href="https://raw.githubusercontent.com/craffael/lehrfempp/master/doc/pfem/hierarchical_basis.pdf" target="_blank"><b>here</b></a>.
 *
 * @see ScalarReferenceFiniteElement
 */
// clang-format on
template <typename SCALAR>
class FeHierarchicSegment final : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeHierarchicSegment(const FeHierarchicSegment &) = default;
  FeHierarchicSegment(FeHierarchicSegment &&) noexcept = default;
  FeHierarchicSegment &operator=(const FeHierarchicSegment &) = default;
  FeHierarchicSegment &operator=(FeHierarchicSegment &&) noexcept = default;
  ~FeHierarchicSegment() override = default;

  FeHierarchicSegment(unsigned degree, const quad::QuadRuleCache &qr_cache)
      : ScalarReferenceFiniteElement<SCALAR>(),
        degree_(degree),
        qr_dual_(&qr_cache.Get(RefEl(), 2 * (degree - 1))) {}

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
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t) const
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim) const override {
    LF_ASSERT_MSG(codim >= 0 && codim <= 1, "codim out of range.");
    return codim == 0 ? degree_ - 1 : 2;
  }

  // clang-format off
  /**
   * @brief One shape function for each vertex, p-1 shape functions for the
   * segment
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t, sub_idx_t) const
   */
  // clang-format on
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t subidx) const override {
    LF_ASSERT_MSG((codim == 0 && subidx == 0) || (codim == 1 && subidx < 2),
                  "codim/subidx out of range.");
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
    result.row(1) = refcoords;
    // Get the shape functions associated with the interior of the segment
    for (int i = 0; i < degree_ - 1; ++i) {
      result.row(i + 2) = refcoords.unaryExpr(
          [&](double x) -> SCALAR { return ILegendre(i + 2, x); });
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
      result.row(i + 2) = refcoords.unaryExpr(
          [&](double x) -> SCALAR { return ILegendreDx(i + 2, x); });
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
    nodes.block(0, 2, 1, qr_dual_->NumPoints()) = qr_dual_->Points();
    return nodes;
  }

  /**
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    return qr_dual_->NumPoints() + 2;
  }

  // clang-format off
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
   * P_{i-1}'(x)f(x) \,\mathrm{d}x \right] &\mbox{ for } i \geq 2 \end{cases} \f]
   */
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const override {
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(NumRefShapeFunctions());
    // Compute the first and second basis function coefficients
    dofs[0] = nodevals[0];
    dofs[1] = nodevals[1];
    // Compute the other basis function coefficients
    for (lf::base::size_type i = 2; i < NumRefShapeFunctions(); ++i) {
      const SCALAR P0 = ILegendreDx(i, 0);
      const SCALAR P1 = ILegendreDx(i, 1);
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
          qr_dual_->Points().unaryExpr(
              [&](double x) -> SCALAR { return LegendreDx(i - 1, x); });
      // Evaluate the integral from the dual basis
      const SCALAR integ =
          (qr_dual_->Weights().transpose().array() * psidd.array() *
           nodevals.tail(qr_dual_->NumPoints()).array())
              .sum();
      // Compute the basis function coefficient
      dofs[i] = (P1 * nodevals[1] - P0 * nodevals[0] - integ) * (2 * i - 1);
    }
    return dofs;
  }

 private:
  unsigned degree_;
  const lf::quad::QuadRule *qr_dual_;
};

// clang-format off
/**
 * @headerfile lf/fe/fe.h
 * @brief Hierarchic Finite Elements of arbitrary degree on triangles
 *
 * The shape functions associated with the vertices are given by the barycentric
 * coordinates \f[ \begin{align*}
 *	\widehat{b^{\cdot}}^1(\widehat{x}, \widehat{y}) &= \lambda_1 = 1 -
 * \widehat{x} - \widehat{y} \\
 *	\widehat{b^{\cdot}}^2(\widehat{x}, \widehat{y}) &= \lambda_2 =
 * \widehat{x} \\ \widehat{b^{\cdot}}^3(\widehat{x}, \widehat{y}) &= \lambda_3 =
 * \widehat{y} \end{align*} \f]
 *
 * The basis functions associated with the triangle edges are given by the
 * homogenized integrated Legendre polynomials \f[ \begin{align*}
 *	\widehat{b^-}^i(\widehat{x}, \widehat{y}) &= (\lambda_1 + \lambda_2)^i
 * L_i\left(\frac{\lambda_2}{\lambda_1+\lambda_2}\right) &\mbox{ for edge 1 and }
 * i = 2, \cdots, p \\
 *	\widehat{b^-}^i(\widehat{x}, \widehat{y}) &= (\lambda_2 + \lambda_3)^i
 * L_i\left(\frac{\lambda_3}{\lambda_2+\lambda_3}\right) &\mbox{ for edge 2 and }
 * i = 2, \cdots, p \\ \widehat{b^-}^i(\widehat{x}, \widehat{y}) &= (\lambda_3 +
 * \lambda_1)^i L_i\left(\frac{\lambda_1}{\lambda_3+\lambda_1}\right) &\mbox{ for
 * edge 3 and } i = 2, \cdots, p \end{align*} \f] Note that the basis function on
 * a specific edge is always zero on the other two edges. This is needed to
 * guarantee continuity of the function space.
 *
 * The basis functions associated with the interior of the triangle are given by
 * the edge basis functions multiplied with an integrated Jacobi polynomial to
 * force the value of the basis function to be zero on all edges. \f[
 *  \widehat{b^{\triangle}}^{ij}(\widehat{x}, \widehat{y}) = (\lambda_1 +
 * \lambda_2)^i L_i\left(\frac{\lambda_2}{\lambda_1+\lambda_2}\right)
 * L_j^{2i}(\lambda_3) \quad\mbox{ for } i \geq 2, j \geq 1, i+j = 3, \cdots, p
 * \f]
 * where \f$L_i : [0, 1] \to \mathbb{R}\f$ and \f$L_i^{\alpha} : [0, 1] \to
 * \mathbb{R}\f$ are the integrated shifted Legendre and integrated shifted
 * Jacobi polynomials respectively: \f[ \begin{align*}
 *	L_i(x) &= \int_0^x\! P_{i-1}(\xi) \,\mathrm{d}\xi \\
 *	L_i^{\alpha}(x) &= \int_0^x\! P_{i-1}^{(\alpha, 0)}(\xi) \,\mathrm{d}\xi
 *  \end{align*}
 * \f]
 *
 * To compute the basis function coefficients from point evaluations of a
 * function, we make use of the dual basis. For the vertices it is simply given
 * by \f[ \lambda_i^{\cdot}[f] = \begin{cases}
 *	f(0, 0) &\mbox{ for } i = 1 \\
 *	f(1, 0) &\mbox{ for } i = 2 \\
 *	f(0, 1) &\mbox{ for } i = 3
 *  \end{cases}
 * \f]
 * For the dual basis on the edges, we simply apply the segment dual basis along
 * the edges of the triangle \f[ \lambda_i^-[f] = \begin{cases}
 *	\frac{1}{2i-1} \left[ P_{i-1}(1)f(1, 0) - P_{i-1}(0)f(0, 0) - \int_0^1\!
 * P_{i-1}'(x)f(x, 0) \,\mathrm{d}x \right] &\mbox{ for edge 1} \\
 *	\frac{1}{2i-1} \left[ P_{i-1}(1)f(0, 1) - P_{i-1}(0)f(1, 0) - \int_0^1\!
 * P_{i-1}'(x)f(1-x, x) \,\mathrm{d}{x} \right] &\mbox{ for edge 2} \\
 *	\frac{1}{2i-1} \left[ P_{i-1}(1)f(0, 0) - P_{i-1}(0)f(0, 1) - \int_0^1\!
 * P_{i-1}'(x)f(0, 1-x) \,\mathrm{d}x \right] &\mbox{ for edge 2} \end{cases} \f]
 * The dual basis for the interior shape functions is quite a bit more involved.
 * It is given by \f[ \lambda_{ij}^{\triangle}[f] = \frac{1}{(2i-1)(2i+2j-1)}
 * \left[ \int_0^1\! \int_0^{1-y}\! f(x, y) \left( (1-y)^{i-1}{L_j^{2i}}''(y) -
 * 2i(1-y)^{i-2}{L_j^{2i}}'(y) \right) L_{j+1}''\left(\frac{x}{1-y}\right)
 * \,\mathrm{d}x \,\mathrm{d}y \right] \f] and must additionally be
 * orthogonalized with respect to the dual basis on the vertices and edges of the
 * triangle by Gram Schmidt.
 *
 * @attention Note that for the basis functions associated with the edges,
 * depending on the `lf::mesh::Orientation` of the according edge, the local
 * coordinate may be flipped to ensure continuity of the function space over the
 * cell interfaces of the mesh. The basis functions and the dual basis must be
 * adjusted accordingly in this case.
 *
 * A complete description of the basis functions and dual basis can be found
 * <a href="https://raw.githubusercontent.com/craffael/lehrfempp/master/doc/pfem/hierarchical_basis.pdf" target="_blank"><b>here</b></a>.
 *
 * @see ScalarReferenceFiniteElement
 */
// clang-format on
template <typename SCALAR>
class FeHierarchicTria final : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeHierarchicTria(const FeHierarchicTria &) = default;
  FeHierarchicTria(FeHierarchicTria &&) noexcept = default;
  FeHierarchicTria &operator=(const FeHierarchicTria &) = default;
  FeHierarchicTria &operator=(FeHierarchicTria &&) noexcept = default;
  ~FeHierarchicTria() override = default;

  FeHierarchicTria(unsigned interior_degree,
                   std::array<unsigned, 3> edge_degrees,
                   const quad::QuadRuleCache &qr_cache,
                   std::span<const lf::mesh::Orientation> rel_orient)
      : ScalarReferenceFiniteElement<SCALAR>(),
        interior_degree_(interior_degree),
        edge_degrees_(edge_degrees),
        qr_dual_edge_(
            {&qr_cache.Get(base::RefEl::kSegment(), 2 * (edge_degrees_[0] - 1)),
             &qr_cache.Get(base::RefEl::kSegment(), 2 * (edge_degrees_[1] - 1)),
             &qr_cache.Get(base::RefEl::kSegment(),
                           2 * (edge_degrees_[2] - 1))}),
        qr_dual_tria_(&qr_cache.Get(lf::base::RefEl::kTria(),
                                    2 * (interior_degree_ - 1))),
        rel_orient_(rel_orient) {
    LF_ASSERT_MSG(interior_degree_ >= 0, "illegal interior degree.");
    LF_ASSERT_MSG(edge_degrees_[0] >= 0, "illegal degree for edge 0");
    LF_ASSERT_MSG(edge_degrees_[1] >= 0, "illegal degree for edge 1");
    LF_ASSERT_MSG(edge_degrees_[2] >= 0, "illegal degree for edge 2");
  }

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kTria();
  }

  [[nodiscard]] unsigned Degree() const override {
    return std::max({interior_degree_, edge_degrees_[0], edge_degrees_[1],
                     edge_degrees_[2], 1U});
  }

  /**
   * @brief The local shape functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
    return NumRefShapeFunctions(0) + NumRefShapeFunctions(1, 0) +
           NumRefShapeFunctions(1, 1) + NumRefShapeFunctions(1, 2) + 3;
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and max(0, (p-2)*(p-1)/2) shape functions on the triangle
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t) const
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim) const override {
    switch (codim) {
      case 0:
        if (interior_degree_ <= 2) {
          return 0;
        } else {
          return (interior_degree_ - 2) * (interior_degree_ - 1) / 2;
        }
      case 1:
        LF_VERIFY_MSG(
            edge_degrees_[0] == edge_degrees_[1] &&
                edge_degrees_[0] == edge_degrees_[2],
            "You cannot call this method with codim=1 if the edges of the "
            "triangle have a differing number of shape functions.");
        return edge_degrees_[0] - 1;
      case 2:
        return 1;
      default:
        LF_ASSERT_MSG(false, "Illegal codim " << codim);
        // Silence compiler warnings
        return 0;
    }
  }

  // clang-format off
  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and max(0, (p-2)*(p-1)/2) shape functions on the triangle
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t, sub_idx_t) const
   */
  // clang-format on
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t subidx) const override {
    switch (codim) {
      case 0:
        LF_ASSERT_MSG(subidx == 0, "illegal codim and subidx.");
        return NumRefShapeFunctions(0);
      case 1:
        LF_ASSERT_MSG(subidx >= 0 && subidx <= 2,
                      "illegal codim/subidx combination.");
        return edge_degrees_[subidx] - 1;
      case 2:
        LF_ASSERT_MSG(subidx >= 0 && subidx <= 2,
                      "illegal codim/subidx combination.");
        return 1;
      default:
        LF_VERIFY_MSG(false, "Illegal codim " << codim);
        // Silence compiler warnings
        return 0;
    }
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        NumRefShapeFunctions(), refcoords.cols());
    // Compute the barycentric coordinate functions
    const Eigen::RowVectorXd l1 = Eigen::RowVectorXd::Ones(refcoords.cols()) -
                                  refcoords.row(0) - refcoords.row(1);
    const Eigen::RowVectorXd l2 = refcoords.row(0);
    const Eigen::RowVectorXd l3 = refcoords.row(1);
    // Get the basis functions associated with the vertices
    result.row(0) = l1;
    result.row(1) = l2;
    result.row(2) = l3;
    int current_dof = 3;  // counter for the next row in the result.

    // Get the basis functions associated with the first edge
    for (unsigned i = 0; i < edge_degrees_[0] - 1; ++i) {
      if (rel_orient_[0] == lf::mesh::Orientation::positive) {
        // L_{i+2}(\lambda_2 ; \lambda_1+\lambda_2)
        for (long j = 0; j < refcoords.cols(); ++j) {
          result(3 + i, j) = ILegendre(i + 2, l2[j], l1[j] + l2[j]);
        }
      } else {
        // L_{i+2}(\lambda_1 ; \lambda_1+\lambda_2)
        for (long j = 0; j < refcoords.cols(); ++j) {
          result(edge_degrees_[0] + 1 - i, j) =
              ILegendre(i + 2, l1[j], l1[j] + l2[j]);
        }
      }
    }
    current_dof += edge_degrees_[0] - 1;

    // Get the basis functions associated with the second edge
    for (unsigned i = 0; i < edge_degrees_[1] - 1; ++i) {
      if (rel_orient_[1] == lf::mesh::Orientation::positive) {
        // L_{i+2}(\lambda_3 ; \lambda_2+\lambda_3)
        for (long j = 0; j < refcoords.cols(); ++j) {
          result(current_dof + i, j) = ILegendre(i + 2, l3[j], l2[j] + l3[j]);
        }
      } else {
        // L_{i+2}(\lambda_2 ; \lambda_2+\lambda_3)
        for (long j = 0; j < refcoords.cols(); ++j) {
          result(current_dof + edge_degrees_[1] - 2 - i, j) =
              ILegendre(i + 2, l2[j], l2[j] + l3[j]);
        }
      }
    }
    current_dof += edge_degrees_[1] - 1;

    // Get the basis functions associated with the third edge
    for (unsigned i = 0; i < edge_degrees_[2] - 1; ++i) {
      if (rel_orient_[2] == lf::mesh::Orientation::positive) {
        // L_{i+2}(\lambda_1 ; \lambda_3+\lambda_1)
        for (long j = 0; j < refcoords.cols(); ++j) {
          result(current_dof + i, j) = ILegendre(i + 2, l1[j], l3[j] + l1[j]);
        }
      } else {
        // L_{i+2}(\lambda_3 ; \lambda_3+\lambda_1)
        for (long j = 0; j < refcoords.cols(); ++j) {
          result(current_dof + edge_degrees_[2] - 2 - i, j) =
              ILegendre(i + 2, l3[j], l3[j] + l1[j]);
        }
      }
    }
    current_dof += edge_degrees_[2] - 1;

    // Get the basis functions associated with the interior of the triangle
    if (interior_degree_ > 2) {
      // i is the degree of the edge function
      for (unsigned i = 0; i < interior_degree_ - 2; ++i) {
        // value of the edge function (agrees with the values computed above,
        // but since the degrees of the edge and interior are not linked at all,
        // we cannot make any assumptions.
        Eigen::Array<SCALAR, 1, Eigen::Dynamic> edge(refcoords.cols());
        for (Eigen::Index j = 0; j < refcoords.cols(); ++j) {
          if (rel_orient_[0] == mesh::Orientation::positive) {
            edge(j) = ILegendre(i + 2, l2[j], l1[j] + l2[j]);
          } else {
            edge(j) = ILegendre(i + 2, l1[j], l1[j] + l2[j]);
          }
        }

        // j is the degree of the blending
        // polynomial
        for (unsigned j = 0; j < interior_degree_ - i - 2; ++j) {
          if (rel_orient_[0] == lf::mesh::Orientation::positive) {
            // L_{i+2}(\lambda_3 ; \lambda_2+\lambda_3) *
            // P_{j+1}^{2i+4}(\lambda_1)
            result.row(current_dof++) =
                (edge * l3.array().unaryExpr([&](double x) -> SCALAR {
                  return IJacobi(j + 1, 2 * i + 4, x);
                })).matrix();
          } else {
            // L_{i+2}(\lambda_2 ; \lambda_2+\lambda_3) *
            // P_{j+1}^{2i+4}(\lambda_1)
            result.row(current_dof++) =
                (edge * l3.array().unaryExpr([&](double x) -> SCALAR {
                  return IJacobi(j + 1, 2 * i + 4, x);
                })).matrix();
          }
        }
      }
    }
    LF_ASSERT_MSG(current_dof == result.rows(),
                  "Something's wrong, not all rows have been filled.");
    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
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
      int current_dof = 3;
      // Get the gradient of the basis functions associated with the vertices
      result(0, 2 * i + 0) = l1_dx[i];
      result(0, 2 * i + 1) = l1_dy[i];
      result(1, 2 * i + 0) = l2_dx[i];
      result(1, 2 * i + 1) = l2_dy[i];
      result(2, 2 * i + 0) = l3_dx[i];
      result(2, 2 * i + 1) = l3_dy[i];
      // Get the gradient of the basis functions associated with the first edge
      for (int j = 0; j < edge_degrees_[0] - 1; ++j) {
        if (rel_orient_[0] == lf::mesh::Orientation::positive) {
          result(current_dof + j, 2 * i + 0) =
              ILegendreDx(j + 2, l2[i], l1[i] + l2[i]) * l2_dx[i] +
              ILegendreDt(j + 2, l2[i], l1[i] + l2[i]) * (l1_dx[i] + l2_dx[i]);
          result(current_dof + j, 2 * i + 1) =
              ILegendreDx(j + 2, l2[i], l1[i] + l2[i]) * l2_dy[i] +
              ILegendreDt(j + 2, l2[i], l1[i] + l2[i]) * (l1_dy[i] + l2_dy[i]);
        } else {
          result(current_dof + edge_degrees_[0] - 2 - j, 2 * i + 0) =
              ILegendreDx(j + 2, l1[i], l1[i] + l2[i]) * l1_dx[i] +
              ILegendreDt(j + 2, l1[i], l1[i] + l2[i]) * (l1_dx[i] + l2_dx[i]);
          result(current_dof + edge_degrees_[0] - 2 - j, 2 * i + 1) =
              ILegendreDx(j + 2, l1[i], l1[i] + l2[i]) * l1_dy[i] +
              ILegendreDt(j + 2, l1[i], l1[i] + l2[i]) * (l1_dy[i] + l2_dy[i]);
        }
      }
      current_dof += edge_degrees_[0] - 1;

      // Get the gradient of the basis functions associated with the second edge
      for (int j = 0; j < edge_degrees_[1] - 1; ++j) {
        if (rel_orient_[1] == lf::mesh::Orientation::positive) {
          result(current_dof + j, 2 * i + 0) =
              ILegendreDx(j + 2, l3[i], l2[i] + l3[i]) * l3_dx[i] +
              ILegendreDt(j + 2, l3[i], l2[i] + l3[i]) * (l2_dx[i] + l3_dx[i]);
          result(current_dof + j, 2 * i + 1) =
              ILegendreDx(j + 2, l3[i], l2[i] + l3[i]) * l3_dy[i] +
              ILegendreDt(j + 2, l3[i], l2[i] + l3[i]) * (l2_dy[i] + l3_dy[i]);
        } else {
          result(current_dof + edge_degrees_[1] - 2 - j, 2 * i + 0) =
              ILegendreDx(j + 2, l2[i], l2[i] + l3[i]) * l2_dx[i] +
              ILegendreDt(j + 2, l2[i], l2[i] + l3[i]) * (l2_dx[i] + l3_dx[i]);
          result(current_dof + edge_degrees_[1] - 2 - j, 2 * i + 1) =
              ILegendreDx(j + 2, l2[i], l2[i] + l3[i]) * l2_dy[i] +
              ILegendreDt(j + 2, l2[i], l2[i] + l3[i]) * (l2_dy[i] + l3_dy[i]);
        }
      }
      current_dof += edge_degrees_[1] - 1;

      // Get the gradient of the basis functions associated with the third edge
      for (int j = 0; j < edge_degrees_[2] - 1; ++j) {
        if (rel_orient_[2] == lf::mesh::Orientation::positive) {
          result(current_dof + j, 2 * i + 0) =
              ILegendreDx(j + 2, l1[i], l3[i] + l1[i]) * l1_dx[i] +
              ILegendreDt(j + 2, l1[i], l3[i] + l1[i]) * (l3_dx[i] + l1_dx[i]);
          result(current_dof + j, 2 * i + 1) =
              ILegendreDx(j + 2, l1[i], l3[i] + l1[i]) * l1_dy[i] +
              ILegendreDt(j + 2, l1[i], l3[i] + l1[i]) * (l3_dy[i] + l1_dy[i]);
        } else {
          result(current_dof + edge_degrees_[2] - 2 - j, 2 * i + 0) =
              ILegendreDx(j + 2, l3[i], l3[i] + l1[i]) * l3_dx[i] +
              ILegendreDt(j + 2, l3[i], l3[i] + l1[i]) * (l3_dx[i] + l1_dx[i]);
          result(current_dof + edge_degrees_[2] - 2 - j, 2 * i + 1) =
              ILegendreDx(j + 2, l3[i], l3[i] + l1[i]) * l3_dy[i] +
              ILegendreDt(j + 2, l3[i], l3[i] + l1[i]) * (l3_dy[i] + l1_dy[i]);
        }
      }
      current_dof += edge_degrees_[2] - 1;

      // Get the gradient of the basis functions associated with the interior of
      // the triangle
      if (interior_degree_ > 2) {
        for (unsigned j = 0; j < interior_degree_ - 2; ++j) {
          SCALAR edge_eval;
          SCALAR edge_dx;
          SCALAR edge_dy;
          if (rel_orient_[0] == lf::mesh::Orientation::positive) {
            edge_eval = ILegendre(j + 2, l2[i], l1[i] + l2[i]);
            edge_dx = ILegendreDx(j + 2, l2[i], l1[i] + l2[i]) * l2_dx[i] +
                      ILegendreDt(j + 2, l2[i], l1[i] + l2[i]) *
                          (l1_dx[i] + l2_dx[i]);
            edge_dy = ILegendreDx(j + 2, l2[i], l1[i] + l2[i]) * l2_dy[i] +
                      ILegendreDt(j + 2, l2[i], l1[i] + l2[i]) *
                          (l1_dy[i] + l2_dy[i]);
          } else {
            edge_eval = ILegendre(j + 2, l1[i], l1[i] + l2[i]);
            edge_dx = ILegendreDx(j + 2, l1[i], l1[i] + l2[i]) * l1_dx[i] +
                      ILegendreDt(j + 2, l1[i], l1[i] + l2[i]) *
                          (l1_dx[i] + l2_dx[i]);
            edge_dy = ILegendreDx(j + 2, l1[i], l1[i] + l2[i]) * l1_dy[i] +
                      ILegendreDt(j + 2, l1[i], l1[i] + l2[i]) *
                          (l1_dy[i] + l2_dy[i]);
          }
          for (unsigned k = 0; k < interior_degree_ - j - 2; ++k) {
            SCALAR jackinte = IJacobi(k + 1, 2 * j + 4, l3[i]);
            SCALAR jackeval = IJacobiDx(k + 1, 2 * j + 4, l3[i]);
            result(current_dof, 2 * i + 0) =
                jackinte * edge_dx + edge_eval * jackeval * l3_dx[i];
            result(current_dof++, 2 * i + 1) =
                jackinte * edge_dy + edge_eval * jackeval * l3_dy[i];
          }
        }
      }
      LF_ASSERT_MSG(current_dof == NumRefShapeFunctions(),
                    "internal error, not all rows have been filled.");
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
    const auto Ns0 = edge_degrees_[0] > 1 ? qr_dual_edge_[0]->NumPoints() : 0;
    const auto Ns1 = edge_degrees_[1] > 1 ? qr_dual_edge_[1]->NumPoints() : 0;
    const auto Ns2 = edge_degrees_[2] > 1 ? qr_dual_edge_[2]->NumPoints() : 0;
    const lf::base::size_type Nt = qr_dual_tria_->NumPoints();
    Eigen::MatrixXd nodes(2, 3 + Ns0 + Ns1 + Ns2 + Nt);
    // Add the vertices
    nodes(0, 0) = 0;
    nodes(1, 0) = 0;
    nodes(0, 1) = 1;
    nodes(1, 1) = 0;
    nodes(0, 2) = 0;
    nodes(1, 2) = 1;
    // Add the quadrature points on the first edge
    if (Ns0 > 0) {
      if (rel_orient_[0] == lf::mesh::Orientation::positive) {
        nodes.block(0, 3, 1, Ns0) = qr_dual_edge_[0]->Points();
      } else {
        nodes.block(0, 3, 1, Ns0) =
            Eigen::RowVectorXd::Ones(Ns0) - qr_dual_edge_[0]->Points();
      }
      nodes.block(1, 3, 1, Ns0).setZero();
    }
    // Add the quadrature points on the second edge
    if (Ns1 > 0) {
      if (rel_orient_[1] == lf::mesh::Orientation::positive) {
        nodes.block(0, 3 + Ns0, 1, Ns1) =
            Eigen::RowVectorXd::Ones(Ns1) - qr_dual_edge_[1]->Points();
        nodes.block(1, 3 + Ns0, 1, Ns1) = qr_dual_edge_[1]->Points();
      } else {
        nodes.block(0, 3 + Ns0, 1, Ns1) = qr_dual_edge_[1]->Points();
        nodes.block(1, 3 + Ns0, 1, Ns1) =
            Eigen::RowVectorXd::Ones(Ns1) - qr_dual_edge_[1]->Points();
      }
    }
    // Add the quadrature points on the third edge
    if (Ns2 > 0) {
      nodes.block(0, 3 + Ns0 + Ns1, 1, Ns2).setZero();
      if (rel_orient_[2] == lf::mesh::Orientation::positive) {
        nodes.block(1, 3 + Ns0 + Ns1, 1, Ns2) =
            Eigen::RowVectorXd::Ones(Ns2) - qr_dual_edge_[2]->Points();
      } else {
        nodes.block(1, 3 + Ns0 + Ns1, 1, Ns2) = qr_dual_edge_[2]->Points();
      }
    }
    if (Nt > 0) {
      // Add the quadrature points for the interior
      nodes.block(0, 3 + Ns0 + Ns1 + Ns2, 2, Nt) = qr_dual_tria_->Points();
    }
    return nodes;
  }

  /**
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    const auto Ns0 = edge_degrees_[0] > 1 ? qr_dual_edge_[0]->NumPoints() : 0;
    const auto Ns1 = edge_degrees_[1] > 1 ? qr_dual_edge_[1]->NumPoints() : 0;
    const auto Ns2 = edge_degrees_[2] > 1 ? qr_dual_edge_[2]->NumPoints() : 0;
    const lf::base::size_type Nt = qr_dual_tria_->NumPoints();
    return 3 + Ns0 + Ns1 + Ns2 + Nt;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const override {
    const auto d0 = edge_degrees_[0] - 1;
    const auto d1 = edge_degrees_[1] - 1;
    const auto d2 = edge_degrees_[2] - 1;

    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(NumRefShapeFunctions());
    // Compute the basis function coefficients of the vertex basis functions
    dofs.segment(0, 3) = NodalValuesToVertexDofs(nodevals);
    // Compute the basis function coefficients on the edges
    dofs.segment(3, d0 + d1 + d2) = NodalValuesToEdgeDofs(nodevals);
    // Compute the basis function coefficients for the face bubbles
    dofs.segment(3 + d0 + d1 + d2, NumRefShapeFunctions(0)) =
        NodalValuesToFaceDofs(nodevals);
    // We need to orthogonalize the face bubble dual basis w.r.t. the
    // dual basis on the vertices and edges
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
        basis_functions = EvalReferenceShapeFunctions(EvaluationNodes());
    const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> boundary_function =
        dofs.segment(0, 3 + d0 + d1 + d2) *
        basis_functions.block(0, 0, 3 + d0 + d1 + d2, basis_functions.cols());
    const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> boundary_face_dofs =
        NodalValuesToFaceDofs(boundary_function);
    dofs.segment(3 + d0 + d1 + d2, NumRefShapeFunctions(0)) -=
        boundary_face_dofs;
    return dofs;
  }

 private:
  unsigned interior_degree_;              // degree for inner shape functions
  std::array<unsigned, 3> edge_degrees_;  // degree of the edges.
  std::array<const quad::QuadRule *, 3>
      qr_dual_edge_;  // Quadrature rules for the edges.
  const lf::quad::QuadRule *qr_dual_tria_;
  std::span<const lf::mesh::Orientation> rel_orient_;

  /*
   * @brief Compute the DOFs of the vertex functions from some
   * given nodal evaluations
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>
  NodalValuesToVertexDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const {
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(3);
    dofs[0] = nodevals[0];
    dofs[1] = nodevals[1];
    dofs[2] = nodevals[2];
    return dofs;
  }

  /*
   * @brief Compute the DOFs of the edge functions from some
   * given nodal evaluations
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToEdgeDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const {
    const auto Ns0 = edge_degrees_[0] > 1 ? qr_dual_edge_[0]->NumPoints() : 0;
    const auto Ns1 = edge_degrees_[1] > 1 ? qr_dual_edge_[1]->NumPoints() : 0;
    const auto Ns2 = edge_degrees_[2] > 1 ? qr_dual_edge_[2]->NumPoints() : 0;
    const lf::base::size_type Nt = qr_dual_tria_->NumPoints();
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(
        edge_degrees_[0] + edge_degrees_[1] + edge_degrees_[2] - 3);
    // Compute the basis function coefficients on the edges
    // by applying the dual basis of the segment

    // first edge:
    for (base::size_type i = 2; i < edge_degrees_[0] + 1; ++i) {
      const SCALAR P0 = ILegendreDx(i, 0);
      const SCALAR P1 = ILegendreDx(i, 1);
      Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
          qr_dual_edge_[0]->Points().unaryExpr(
              [&](double x) -> SCALAR { return LegendreDx(i - 1, x); });

      const SCALAR integ1 = (qr_dual_edge_[0]->Weights().transpose().array() *
                             psidd.array() * nodevals.segment(3, Ns0).array())
                                .sum();
      if (rel_orient_[0] == lf::mesh::Orientation::positive) {
        dofs[1 + i - 3] =
            (P1 * nodevals[1] - P0 * nodevals[0] - integ1) * (2 * i - 1);
      } else {
        dofs[4 - i + (edge_degrees_[0] - 1) - 3] =
            (P1 * nodevals[0] - P0 * nodevals[1] - integ1) * (2 * i - 1);
      }
    }

    // Compute the basis function coefficients for the second edge
    for (base::size_type i = 2; i < edge_degrees_[1] + 1; ++i) {
      const SCALAR P0 = ILegendreDx(i, 0);
      const SCALAR P1 = ILegendreDx(i, 1);
      Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
          qr_dual_edge_[1]->Points().unaryExpr(
              [&](double x) -> SCALAR { return LegendreDx(i - 1, x); });

      const SCALAR integ2 =
          (qr_dual_edge_[1]->Weights().transpose().array() * psidd.array() *
           nodevals.segment(3 + Ns0, Ns1).array())
              .sum();
      if (rel_orient_[1] == lf::mesh::Orientation::positive) {
        dofs[1 + i + (edge_degrees_[0] - 1) - 3] =
            (P1 * nodevals[2] - P0 * nodevals[1] - integ2) * (2 * i - 1);
      } else {
        dofs[4 - i + (edge_degrees_[0] + edge_degrees_[1] - 2) - 3] =
            (P1 * nodevals[1] - P0 * nodevals[2] - integ2) * (2 * i - 1);
      }
    }

    // Compute the basis function coefficients for the second edge
    for (base::size_type i = 2; i < edge_degrees_[2] + 1; ++i) {
      const SCALAR P0 = ILegendreDx(i, 0);
      const SCALAR P1 = ILegendreDx(i, 1);
      Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
          qr_dual_edge_[2]->Points().unaryExpr(
              [&](double x) -> SCALAR { return LegendreDx(i - 1, x); });

      const SCALAR integ3 =
          (qr_dual_edge_[2]->Weights().transpose().array() * psidd.array() *
           nodevals.segment(3 + Ns0 + Ns1, Ns2).array())
              .sum();
      if (rel_orient_[2] == lf::mesh::Orientation::positive) {
        dofs[1 + i + (edge_degrees_[0] + edge_degrees_[1] - 2) - 3] =
            (P1 * nodevals[0] - P0 * nodevals[2] - integ3) * (2 * i - 1);
      } else {
        dofs[4 - i +
             (edge_degrees_[0] + edge_degrees_[1] + edge_degrees_[2] - 3) - 3] =
            (P1 * nodevals[2] - P0 * nodevals[0] - integ3) * (2 * i - 1);
      }
    }

    return dofs;
  }

  /*
   * @brief Compute the non-orthogonalized DOFs of the face bubbles from some
   * given nodal evaluations
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToFaceDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const {
    const auto Ns0 = edge_degrees_[0] > 1 ? qr_dual_edge_[0]->NumPoints() : 0;
    const auto Ns1 = edge_degrees_[1] > 1 ? qr_dual_edge_[1]->NumPoints() : 0;
    const auto Ns2 = edge_degrees_[2] > 1 ? qr_dual_edge_[2]->NumPoints() : 0;
    const lf::base::size_type Ns = Ns0 + Ns1 + Ns2;
    const lf::base::size_type Nt = qr_dual_tria_->NumPoints();
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(NumRefShapeFunctions(0));
    if (interior_degree_ > 2) {
      unsigned idx = 0;
      const Eigen::Matrix<double, 1, Eigen::Dynamic> xnorm =
          (qr_dual_tria_->Points().row(0).array() /
           qr_dual_tria_->Points().row(1).array().unaryExpr([&](double y) {
             return 1 - y;
           })).matrix();
      // We need to flip the local coordinates in case the
      // relative orientation is negative
      const Eigen::Matrix<double, 1, Eigen::Dynamic> xnorm_adj =
          rel_orient_[0] == lf::mesh::Orientation::positive
              ? xnorm
              : Eigen::RowVectorXd::Ones(Nt) - xnorm;
      const Eigen::Matrix<double, 1, Eigen::Dynamic> y =
          qr_dual_tria_->Points().row(1);
      // i is the degree of the edge function
      for (unsigned i = 0; i < interior_degree_ - 2; ++i) {
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
            xnorm_adj.unaryExpr(
                [&](double x) -> SCALAR { return LegendreDx(i + 1, x); });
        // j is the degree of the blending Jacobi polynomial
        for (unsigned j = 0; j < interior_degree_ - i - 2; ++j) {
          const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> jacdd =
              qr_dual_tria_->Points().row(1).unaryExpr([&](double y) -> SCALAR {
                return JacobiDx(j, 2 * i + 4, y);
              });
          const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> jacd =
              qr_dual_tria_->Points().row(1).unaryExpr([&](double y) -> SCALAR {
                return IJacobiDx(j + 1, 2 * i + 4, y);
              });
          dofs[idx] =
              (qr_dual_tria_->Weights().transpose().array() *
               nodevals.block(0, 3 + Ns, 1, Nt).array() * psidd.array() *
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

// clang-format off
/**
 * @headerfile lf/fe/fe.h
 * @brief Hierarchic Finite Elements of arbitrary degree on quadrilaterals
 *
 * The basis functions on the quadrilateral has a tensor product structure and
 * can thus be represented by products of basis functions on segments.
 * The vertex basis functions on the reference quad are therefore given by
 * \f[
 * \begin{align*}
 *  \widehat{b^{\cdot}}^0(x, y) &:= (1 - x)(1 - y) \\
 *  \widehat{b^{\cdot}}^1(x, y) &:= x(1 - y) \\
 *  \widehat{b^{\cdot}}^2(x, y) &:= xy \\
 *  \widehat{b^{\cdot}}^3(x, y) &:= (1 - x)y.
 * \end{align*}
 * \f]
 * The edge basis functions can be written as
 * \f[
 *  \begin{align*}
 *	\widehat{b^{-}}^{0,n} &:= (1-y)L_n(x) \\
 *	\widehat{b^{-}}^{1,n} &:= xL_n(y) \\
 *	\widehat{b^{-}}^{2,n} &:= yL_n(1-x) \\
 *	\widehat{b^{-}}^{3,n} &:= (1-x)L_n(1-y)
 *  \end{align*}
 * \f]
 * where \f$ n \geq 2 \f$ is the degree of the basis function.
 * Finally, the face bubbles are given by
 * \f[
 *  \widehat{b^{\square}}^{n,m}(x, y) := L_n(x)L_m(y)
 * \f]
 * where \f$ n \geq 2 \f$, \f$ m \geq 2 \f$.
 *
 * The dual basis is therefore also quite simple, as we can recycle the one from
 * the segments by first applying the dual basis along the \f$x\f$-axis and then
 * apply the dual basis to the resulting 1d function.
 *
 * @attention Note that for the basis functions associated with the edges,
 * depending on the `lf::mesh::Orientation` of the according edge, the local
 * coordinate may be flipped to ensure continuity of the function space over the
 * cell interfaces of the mesh. The basis functions and the dual basis must be
 * adjusted accordingly in this case.
 *
 * A complete description of the basis functions and dual basis can be found
 * <a href="https://raw.githubusercontent.com/craffael/lehrfempp/master/doc/pfem/hierarchical_basis.pdf" target="_blank"><b>here</b></a>.
 *
 * @see ScalarReferenceFiniteElement
 */
// clang-format on
template <typename SCALAR>
class FeHierarchicQuad final : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeHierarchicQuad(const FeHierarchicQuad &) = default;
  FeHierarchicQuad(FeHierarchicQuad &&) noexcept = default;
  FeHierarchicQuad &operator=(const FeHierarchicQuad &) = default;
  FeHierarchicQuad &operator=(FeHierarchicQuad &&) noexcept = default;
  ~FeHierarchicQuad() override = default;

  FeHierarchicQuad(unsigned interior_degree,
                   std::array<unsigned, 4> edge_degrees,
                   const quad::QuadRuleCache &qr_cache,
                   std::span<const lf::mesh::Orientation> rel_orient)
      : ScalarReferenceFiniteElement<SCALAR>(),
        interior_degree_(interior_degree),
        edge_degrees_(edge_degrees),
        qr_dual_edge_({&qr_cache.Get(lf::base::RefEl::kSegment(),
                                     2 * (edge_degrees_[0] - 1)),
                       &qr_cache.Get(lf::base::RefEl::kSegment(),
                                     2 * (edge_degrees_[1] - 1)),
                       &qr_cache.Get(lf::base::RefEl::kSegment(),
                                     2 * (edge_degrees_[2] - 1)),
                       &qr_cache.Get(lf::base::RefEl::kSegment(),
                                     2 * (edge_degrees_[3] - 1))}),
        rel_orient_(rel_orient),
        fe1d_(std::max({interior_degree_, edge_degrees_[0], edge_degrees_[1],
                        edge_degrees_[2], edge_degrees_[3]}),
              qr_cache) {
    LF_ASSERT_MSG(interior_degree_ >= 0, "illegal interior degree.");
    LF_ASSERT_MSG(edge_degrees_[0] >= 0, "illegal degree for edge 0");
    LF_ASSERT_MSG(edge_degrees_[1] >= 0, "illegal degree for edge 1");
    LF_ASSERT_MSG(edge_degrees_[2] >= 0, "illegal degree for edge 2");
    LF_ASSERT_MSG(edge_degrees_[3] >= 0, "illegal degree for edge 3");
  }

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kQuad();
  }

  [[nodiscard]] unsigned Degree() const override {
    return std::max({interior_degree_, edge_degrees_[0], edge_degrees_[1],
                     edge_degrees_[2], edge_degrees_[3]});
  }

  /**
   * @brief The local shape functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions() const override {
    return 4 + NumRefShapeFunctions(0) + NumRefShapeFunctions(1, 0) +
           NumRefShapeFunctions(1, 1) + NumRefShapeFunctions(1, 2) +
           NumRefShapeFunctions(1, 3);
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and (p-1)^2 shape functions on the quadrilateral
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t) const
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim) const override {
    switch (codim) {
      case 0:
        return (interior_degree_ - 1) * (interior_degree_ - 1);
      case 1:
        LF_VERIFY_MSG(
            edge_degrees_[0] == edge_degrees_[1] &&
                edge_degrees_[0] == edge_degrees_[2] &&
                edge_degrees_[0] == edge_degrees_[3],
            "You cannot call this method with codim=1 if the edges of the "
            "quad have a differing number of shape functions.");
        return edge_degrees_[0] - 1;
      case 2:
        return 1;
      default:
        LF_ASSERT_MSG(false, "Illegal codim " << codim);
        // Silence compiler warnings
        return 0;
    }
  }

  // clang-format off
  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and (p-1)^2 shape functions on the quadrilateral
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t) const
   */
  // clang-format on
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t subidx) const override {
    switch (codim) {
      case 0:
        LF_ASSERT_MSG(subidx == 0, "illegal codim and subidex.");
        return NumRefShapeFunctions(0);
      case 1:
        LF_ASSERT_MSG(subidx >= 0 && subidx < 4, "illegal codim and subidx.");
        return edge_degrees_[subidx] - 1;
      case 2:
        LF_ASSERT_MSG(subidx >= 0 && subidx < 4, "illegal codim and subidx.");
        return 1;
      default:
        LF_VERIFY_MSG(false, "Illegal codim " << codim);
        // Silence compiler warnings
        return 0;
    }
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

    int current_dof = 4;

    // Get the basis functions associated with the first edge
    // as a tensor product of 1D basis functions
    for (int i = 0; i < edge_degrees_[0] - 1; ++i) {
      if (rel_orient_[0] == lf::mesh::Orientation::positive) {
        result.row(current_dof + i) =
            (sf1d_x.row(2 + i).array() * sf1d_y.row(0).array()).matrix();
      } else {
        result.row(current_dof - 2 + edge_degrees_[0] - i) =
            (sf1df_x.row(2 + i).array() * sf1d_y.row(0).array()).matrix();
      }
    }
    current_dof += edge_degrees_[0] - 1;
    // Get the basis functions associated with the second edge
    // as a tensor product of 1D basis functions
    for (int i = 0; i < edge_degrees_[1] - 1; ++i) {
      if (rel_orient_[1] == lf::mesh::Orientation::positive) {
        result.row(current_dof + i) =
            (sf1d_x.row(1).array() * sf1d_y.row(2 + i).array()).matrix();
      } else {
        result.row(current_dof - 2 + edge_degrees_[1] - i) =
            (sf1d_x.row(1).array() * sf1df_y.row(2 + i).array()).matrix();
      }
    }
    current_dof += edge_degrees_[1] - 1;
    // Get the basis functions associated with the third edge
    // as a tensor product of 1D basis functions
    for (int i = 0; i < edge_degrees_[2] - 1; ++i) {
      if (rel_orient_[2] == lf::mesh::Orientation::positive) {
        result.row(current_dof + i) =
            (sf1df_x.row(2 + i).array() * sf1d_y.row(1).array()).matrix();
      } else {
        result.row(current_dof + edge_degrees_[2] - 2 - i) =
            (sf1d_x.row(2 + i).array() * sf1d_y.row(1).array()).matrix();
      }
    }
    current_dof += edge_degrees_[2] - 1;
    // Get the basis functions associated with the fourth edge
    // as a tensor product of 1D basis functions
    for (int i = 0; i < edge_degrees_[3] - 1; ++i) {
      if (rel_orient_[3] == lf::mesh::Orientation::positive) {
        result.row(current_dof + i) =
            (sf1d_x.row(0).array() * sf1df_y.row(2 + i).array()).matrix();
      } else {
        result.row(current_dof + edge_degrees_[3] - 2 - i) =
            (sf1d_x.row(0).array() * sf1d_y.row(2 + i).array()).matrix();
      }
    }
    current_dof += edge_degrees_[3] - 1;
    // Get the basis functions associated with the interior of the quad
    // as a tensor product of 1D basis functions
    for (int i = 0; i < interior_degree_ - 1; ++i) {
      for (int j = 0; j < interior_degree_ - 1; ++j) {
        result.row(current_dof++) =
            (sf1d_x.row(j + 2).array() * sf1d_y.row(i + 2).array()).matrix();
      }
    }
    LF_ASSERT_MSG(current_dof == result.rows(),
                  "Not all rows have been filled.");
    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
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

      int current_dof = 4;
      // Get the basis functions associated with the first edge
      for (int j = 0; j < edge_degrees_[0] - 1; ++j) {
        if (rel_orient_[0] == lf::mesh::Orientation::positive) {
          result(current_dof + j, 2 * i + 0) = sf1d_dx(2 + j, i) * sf1d_y(0, i);
          result(current_dof + j, 2 * i + 1) = sf1d_x(2 + j, i) * sf1d_dy(0, i);
        } else {
          result(current_dof + edge_degrees_[0] - 2 - j, 2 * i + 0) =
              -sf1df_dx(2 + j, i) * sf1d_y(0, i);
          result(current_dof + edge_degrees_[0] - 2 - j, 2 * i + 1) =
              sf1df_x(2 + j, i) * sf1d_dy(0, i);
        }
      }
      current_dof += edge_degrees_[0] - 1;
      // Get the basis functions associated with the second edge
      for (int j = 0; j < edge_degrees_[1] - 1; ++j) {
        if (rel_orient_[1] == lf::mesh::Orientation::positive) {
          result(current_dof + j, 2 * i + 0) = sf1d_dx(1, i) * sf1d_y(2 + j, i);
          result(current_dof + j, 2 * i + 1) = sf1d_x(1, i) * sf1d_dy(2 + j, i);
        } else {
          result(current_dof + edge_degrees_[1] - 2 - j, 2 * i + 0) =
              sf1d_dx(1, i) * sf1df_y(2 + j, i);
          result(current_dof + edge_degrees_[1] - 2 - j, 2 * i + 1) =
              sf1d_x(1, i) * -sf1df_dy(2 + j, i);
        }
      }
      current_dof += edge_degrees_[1] - 1;
      // Get the basis functions associated with the third edge
      for (int j = 0; j < edge_degrees_[2] - 1; ++j) {
        if (rel_orient_[2] == lf::mesh::Orientation::positive) {
          result(current_dof + j, 2 * i + 0) =
              -sf1df_dx(2 + j, i) * sf1d_y(1, i);
          result(current_dof + j, 2 * i + 1) =
              sf1df_x(2 + j, i) * sf1d_dy(1, i);
        } else {
          result(current_dof + edge_degrees_[2] - 2 - j, 2 * i + 0) =
              sf1d_dx(2 + j, i) * sf1d_y(1, i);
          result(current_dof + edge_degrees_[2] - 2 - j, 2 * i + 1) =
              sf1d_x(2 + j, i) * sf1d_dy(1, i);
        }
      }
      current_dof += edge_degrees_[2] - 1;
      // Get the basis functions associated with the fourth edge
      for (int j = 0; j < edge_degrees_[3] - 1; ++j) {
        if (rel_orient_[3] == lf::mesh::Orientation::positive) {
          result(current_dof + j, 2 * i + 0) =
              sf1d_dx(0, i) * sf1df_y(2 + j, i);
          result(current_dof + j, 2 * i + 1) =
              sf1d_x(0, i) * -sf1df_dy(2 + j, i);
        } else {
          result(current_dof + edge_degrees_[3] - 2 - j, 2 * i + 0) =
              sf1d_dx(0, i) * sf1d_y(2 + j, i);
          result(current_dof + edge_degrees_[3] - 2 - j, 2 * i + 1) =
              sf1d_x(0, i) * sf1d_dy(2 + j, i);
        }
      }
      current_dof += edge_degrees_[3] - 1;
      // Get the basis functions associated with the interior of the quad
      for (int j = 0; j < interior_degree_ - 1; ++j) {
        for (int k = 0; k < interior_degree_ - 1; ++k) {
          result(current_dof, 2 * i + 0) = sf1d_dx(k + 2, i) * sf1d_y(j + 2, i);
          result(current_dof++, 2 * i + 1) =
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
    const base::size_type Ne0 =
        edge_degrees_[0] > 1 ? qr_dual_edge_[0]->NumPoints() : 0;
    const base::size_type Ne1 =
        edge_degrees_[1] > 1 ? qr_dual_edge_[1]->NumPoints() : 0;
    const base::size_type Ne2 =
        edge_degrees_[2] > 1 ? qr_dual_edge_[2]->NumPoints() : 0;
    const base::size_type Ne3 =
        edge_degrees_[3] > 1 ? qr_dual_edge_[3]->NumPoints() : 0;
    auto Nq = interior_degree_ > 1 ? fe1d_.NumEvaluationNodes() : 0;

    Eigen::MatrixXd nodes(2, 4 + Ne0 + Ne1 + Ne2 + Ne3 + Nq * Nq);
    // Add the vertices
    nodes(0, 0) = 0;
    nodes(1, 0) = 0;
    nodes(0, 1) = 1;
    nodes(1, 1) = 0;
    nodes(0, 2) = 1;
    nodes(1, 2) = 1;
    nodes(0, 3) = 0;
    nodes(1, 3) = 1;
    if (edge_degrees_[0] > 1) {
      // Add the quadrature points on the first edge
      if (rel_orient_[0] == lf::mesh::Orientation::positive) {
        nodes.block(0, 4, 1, Ne0) = qr_dual_edge_[0]->Points();
      } else {
        nodes.block(0, 4, 1, Ne0) =
            Eigen::RowVectorXd::Ones(Ne0) - qr_dual_edge_[0]->Points();
      }
      nodes.block(1, 4, 1, Ne0).setZero();
    }

    if (edge_degrees_[1] > 1) {
      // Add the quadrature points on the second edge
      nodes.block(0, 4 + Ne0, 1, Ne1).setOnes();
      if (rel_orient_[1] == lf::mesh::Orientation::positive) {
        nodes.block(1, 4 + Ne0, 1, Ne1) = qr_dual_edge_[1]->Points();
      } else {
        nodes.block(1, 4 + Ne0, 1, Ne1) =
            Eigen::RowVectorXd::Ones(Ne1) - qr_dual_edge_[1]->Points();
      }
    }
    if (edge_degrees_[2] > 1) {
      // Add the quadrature points on the third edge
      if (rel_orient_[2] == lf::mesh::Orientation::positive) {
        nodes.block(0, 4 + Ne0 + Ne1, 1, Ne2) =
            Eigen::RowVectorXd::Ones(Ne2) - qr_dual_edge_[2]->Points();
      } else {
        nodes.block(0, 4 + Ne0 + Ne1, 1, Ne2) = qr_dual_edge_[2]->Points();
      }
      nodes.block(1, 4 + Ne0 + Ne1, 1, Ne2).setOnes();
    }
    if (edge_degrees_[3] > 1) {
      // Add the quadrature points on the fourth edge
      nodes.block(0, 4 + Ne0 + Ne1 + Ne2, 1, Ne3).setZero();
      if (rel_orient_[3] == lf::mesh::Orientation::positive) {
        nodes.block(1, 4 + Ne0 + Ne1 + Ne2, 1, Ne3) =
            Eigen::RowVectorXd::Ones(Ne3) - qr_dual_edge_[3]->Points();
      } else {
        nodes.block(1, 4 + Ne0 + Ne1 + Ne2, 1, Ne3) =
            qr_dual_edge_[3]->Points();
      }
    }
    if (interior_degree_ > 1) {
      // Add the quadrature points for the face
      auto points = fe1d_.EvaluationNodes();
      for (lf::base::size_type i = 0; i < Nq; ++i) {
        nodes.block(0, 4 + Ne0 + Ne1 + Ne2 + Ne3 + i * Nq, 1, Nq) = points;
        nodes.block(1, 4 + Ne0 + Ne1 + Ne2 + Ne3 + i * Nq, 1, Nq)
            .setConstant(points(0, i));
      }
    }
    return nodes;
  }

  /**
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    const base::size_type Ne0 =
        edge_degrees_[0] > 1 ? qr_dual_edge_[0]->NumPoints() : 0;
    const base::size_type Ne1 =
        edge_degrees_[1] > 1 ? qr_dual_edge_[1]->NumPoints() : 0;
    const base::size_type Ne2 =
        edge_degrees_[2] > 1 ? qr_dual_edge_[2]->NumPoints() : 0;
    const base::size_type Ne3 =
        edge_degrees_[3] > 1 ? qr_dual_edge_[3]->NumPoints() : 0;
    const lf::base::size_type Nq =
        interior_degree_ > 1 ? fe1d_.NumEvaluationNodes() : 0;
    return 4 + Ne0 + Ne1 + Ne2 + Ne3 + Nq * Nq;
  }

  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> &nodevals) const override {
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> dofs(NumRefShapeFunctions());
    const base::size_type Ne0 =
        edge_degrees_[0] > 1 ? qr_dual_edge_[0]->NumPoints() : 0;
    const base::size_type Ne1 =
        edge_degrees_[1] > 1 ? qr_dual_edge_[1]->NumPoints() : 0;
    const base::size_type Ne2 =
        edge_degrees_[2] > 1 ? qr_dual_edge_[2]->NumPoints() : 0;
    const base::size_type Ne3 =
        edge_degrees_[3] > 1 ? qr_dual_edge_[3]->NumPoints() : 0;

    // Compute the basis function coefficients of the vertex basis functions
    dofs[0] = nodevals[0];
    dofs[1] = nodevals[1];
    dofs[2] = nodevals[2];
    dofs[3] = nodevals[3];
    // Compute the basis function coefficients on the edges
    // by applying the dual basis of the segment
    for (lf::base::size_type i = 2; i < Degree() + 1; ++i) {
      const SCALAR P0 = ILegendreDx(i, 0);
      const SCALAR P1 = ILegendreDx(i, 1);

      // Compute the basis function coefficients for the first edge
      if (i <= edge_degrees_[0]) {
        Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
            qr_dual_edge_[0]->Points().unaryExpr(
                [&](double x) -> SCALAR { return LegendreDx(i - 1, x); });
        const SCALAR integ1 = (qr_dual_edge_[0]->Weights().transpose().array() *
                               psidd.array() * nodevals.segment(4, Ne0).array())
                                  .sum();
        if (rel_orient_[0] == lf::mesh::Orientation::positive) {
          dofs[2 + i] =
              (P1 * nodevals[1] - P0 * nodevals[0] - integ1) * (2 * i - 1);
        } else {
          dofs[5 - i + (edge_degrees_[0] - 1)] =
              (P1 * nodevals[0] - P0 * nodevals[1] - integ1) * (2 * i - 1);
        }
      }
      if (i <= edge_degrees_[1]) {
        // Compute the basis function coefficients for the second edge
        Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
            qr_dual_edge_[1]->Points().unaryExpr(
                [&](double x) -> SCALAR { return LegendreDx(i - 1, x); });
        const SCALAR integ2 =
            (qr_dual_edge_[1]->Weights().transpose().array() * psidd.array() *
             nodevals.segment(4 + Ne0, Ne1).array())
                .sum();
        if (rel_orient_[1] == lf::mesh::Orientation::positive) {
          dofs[2 + i + (edge_degrees_[0] - 1)] =
              (P1 * nodevals[2] - P0 * nodevals[1] - integ2) * (2 * i - 1);
        } else {
          dofs[5 - i + (edge_degrees_[0] + edge_degrees_[1] - 2)] =
              (P1 * nodevals[1] - P0 * nodevals[2] - integ2) * (2 * i - 1);
        }
      }
      if (i <= edge_degrees_[2]) {
        // Compute the basis function coefficients for the third edge
        Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
            qr_dual_edge_[2]->Points().unaryExpr(
                [&](double x) -> SCALAR { return LegendreDx(i - 1, x); });
        const SCALAR integ3 =
            (qr_dual_edge_[2]->Weights().transpose().array() * psidd.array() *
             nodevals.segment(4 + Ne0 + Ne1, Ne2).array())
                .sum();
        if (rel_orient_[2] == lf::mesh::Orientation::positive) {
          dofs[2 + i + (edge_degrees_[0] + edge_degrees_[1] - 2)] =
              (P1 * nodevals[3] - P0 * nodevals[2] - integ3) * (2 * i - 1);
        } else {
          dofs[5 - i +
               (edge_degrees_[0] + edge_degrees_[1] + edge_degrees_[2] - 3)] =
              (P1 * nodevals[2] - P0 * nodevals[3] - integ3) * (2 * i - 1);
        }
      }
      if (i <= edge_degrees_[3]) {
        // Compute the basis function coefficients for the fourth edge
        Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> psidd =
            qr_dual_edge_[3]->Points().unaryExpr(
                [&](double x) -> SCALAR { return LegendreDx(i - 1, x); });
        const SCALAR integ4 =
            (qr_dual_edge_[3]->Weights().transpose().array() * psidd.array() *
             nodevals.segment(4 + Ne0 + Ne1 + Ne2, Ne3).array())
                .sum();
        if (rel_orient_[3] == lf::mesh::Orientation::positive) {
          dofs[2 + i +
               (edge_degrees_[0] + edge_degrees_[1] + edge_degrees_[2] - 3)] =
              (P1 * nodevals[0] - P0 * nodevals[3] - integ4) * (2 * i - 1);
        } else {
          dofs[5 - i +
               (edge_degrees_[0] + edge_degrees_[1] + edge_degrees_[2] +
                edge_degrees_[3] - 4)] =
              (P1 * nodevals[3] - P0 * nodevals[0] - integ4) * (2 * i - 1);
        }
      }
    }
    if (interior_degree_ > 1) {
      // Compute the basis function coefficients for the face bubbles
      auto Nq = fe1d_.NumEvaluationNodes();
      Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dof_temp(
          Nq, interior_degree_ - 1);
      for (unsigned i = 0; i < Nq; ++i) {
        dof_temp.row(i) = fe1d_
                              .NodalValuesToDofs(nodevals.segment(
                                  4 + Ne0 + Ne1 + Ne2 + Ne3 + Nq * i, Nq))
                              .segment(2, interior_degree_ - 1);
      }
      Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dof_temp2(
          interior_degree_ - 1, interior_degree_ - 1);
      for (unsigned i = 0; i < interior_degree_ - 1; ++i) {
        dofs.segment(edge_degrees_[0] + edge_degrees_[1] + edge_degrees_[2] +
                         edge_degrees_[3] + i * (interior_degree_ - 1),
                     interior_degree_ - 1) = dof_temp2.row(i) =
            fe1d_.NodalValuesToDofs(dof_temp.col(i))
                .segment(2, interior_degree_ - 1);
      }
      dofs.segment(edge_degrees_[0] + edge_degrees_[1] + edge_degrees_[2] +
                       edge_degrees_[3],
                   (interior_degree_ - 1) * (interior_degree_ - 1)) =
          Eigen::Map<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>>(
              dof_temp2.data(), 1,
              (interior_degree_ - 1) * (interior_degree_ - 1));
    }
    return dofs;
  }

 private:
  unsigned interior_degree_;
  std::array<unsigned, 4> edge_degrees_;
  std::array<const quad::QuadRule *, 4> qr_dual_edge_;
  FeHierarchicSegment<SCALAR>
      fe1d_;  // degree = max(interior_degree_, edge_degrees_)
  std::span<const lf::mesh::Orientation> rel_orient_;
};

}  // end namespace lf::fe

#endif  // LF_USCALFE_HP_FE_H_
