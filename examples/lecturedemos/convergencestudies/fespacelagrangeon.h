#ifndef EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FESPACELAGRANGEON_H_
#define EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FESPACELAGRANGEON_H_

#define _USE_MATH_DEFINES

#include <lf/assemble/assemble.h>
#include <lf/mesh/mesh.h>
#include <lf/uscalfe/uscalfe.h>

#include <cmath>
#include <memory>

/**
 * @brief Computes Chebyshev interpolation nodes in [0, 1]
 * @param n Degree of the Chebyshev interpolation nodes
 * @returns An Eigen vector containing the interpolation nodes on [0, 1]
 */
Eigen::VectorXd chebyshevNodes(unsigned n) {
  // Compute the chebyshev nodes in the interval [0, 1]
  const auto cosine = [](double x) -> double { return std::cos(x); };
  return (Eigen::VectorXd::Ones(n) +
          Eigen::VectorXd::LinSpaced(n, M_PI - M_PI / (2 * n), M_PI / (2 * n))
              .unaryExpr(cosine)) /
         2;
}

/**
 * @brief Lagrangian Finite Elements of arbitrary degreen on segments
 *
 * The evaluation nodes are given by the Chebyshev nodes and include the
 * endpoints of the segment.
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeONSegment final
    : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeONSegment(const FeLagrangeONSegment &) = default;
  FeLagrangeONSegment(FeLagrangeONSegment &&) = default;
  FeLagrangeONSegment &operator=(const FeLagrangeONSegment &) = default;
  FeLagrangeONSegment &operator=(FeLagrangeONSegment &&) = default;
  ~FeLagrangeONSegment() = default;

  FeLagrangeONSegment(unsigned degree) : degree_(degree), ref_func_coeffs_() {
    eval_nodes_ = ComputeEvaluationNodes(degree);
    ref_func_coeffs_ = ComputePolyBasis(eval_nodes_).inverse().transpose();
  }

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
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim) const override {
    return codim == 0 ? degree_ - 1 : 1;
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions for the
   * segment
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t,
   * lf::base::sub_idx_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim,
      lf::base::sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    const auto poly_basis = ComputePolyBasis(refcoords);
    return ref_func_coeffs_ * poly_basis.transpose();
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd &refcoords) const override {
    const auto dx = ComputePolyBasisDerivative(refcoords);
    return ref_func_coeffs_ * dx.transpose();
  }

  /**
   * @brief Evaluation nodes are the endpoints of the segment and the Chebyshev
   * nodes of degree p-1 on the segment
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return eval_nodes_;
  }

  /**
   * @brief p+1 shape functions
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }

 private:
  unsigned degree_;
  Eigen::MatrixXd eval_nodes_;
  Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ref_func_coeffs_;

  Eigen::MatrixXd ComputeEvaluationNodes(unsigned p) const {
    Eigen::MatrixXd nodes(1, p + 1);
    // Add the evaluation nodes corresponding to the vertices of the segment
    nodes(0, 0) = 0;
    nodes(0, 1) = 1;
    // Add the evaluation nodes corresponding to the interior of the segment
    if (p > 1) {
      nodes.block(0, 2, 1, p - 1) = chebyshevNodes(p - 1).transpose();
    }
    return nodes;
  }

  Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ComputePolyBasis(
      const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> &refcoords)
      const {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> poly_basis(
        refcoords.cols(), degree_ + 1);
    for (unsigned exponent = 0; exponent <= degree_; ++exponent) {
      for (unsigned node_idx = 0; node_idx < refcoords.cols(); ++node_idx) {
        poly_basis(node_idx, exponent) =
            std::pow(refcoords(0, node_idx), exponent);
      }
    }
    return poly_basis;
  }

  Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  ComputePolyBasisDerivative(
      const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> &refcoords)
      const {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dx(refcoords.cols(),
                                                             degree_ + 1);
    dx.col(0).setZero();
    for (unsigned exponent = 1; exponent <= degree_; ++exponent) {
      for (unsigned node_idx = 0; node_idx < refcoords.cols(); ++node_idx) {
        dx(node_idx, exponent) =
            exponent * std::pow(refcoords(0, node_idx), exponent - 1);
      }
    }
    return dx;
  }
};

/**
 * @brief Lagrangian Finite Elements of arbitrary degreen on triangles
 *
 * The evaluation nodes are derived from the Chebyshev interpolation nodes
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeONTria final
    : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeONTria(const FeLagrangeONTria &) = default;
  FeLagrangeONTria(FeLagrangeONTria &&) = default;
  FeLagrangeONTria &operator=(const FeLagrangeONTria &) = default;
  FeLagrangeONTria &operator=(FeLagrangeONTria &&) = default;
  ~FeLagrangeONTria() = default;

  FeLagrangeONTria(unsigned degree)
      : degree_(degree), eval_nodes_(), ref_func_coeffs_() {
    eval_nodes_ = ComputeEvaluationNodes(degree);
    ref_func_coeffs_ = ComputePolyBasis(eval_nodes_).inverse().transpose();
  }

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
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim) const override {
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
    }
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and max(0, (p-2)*(p-1)/2) shape functions on the triangle
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim,
      lf::base::sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    const auto poly_basis = ComputePolyBasis(refcoords);
    return ref_func_coeffs_ * poly_basis.transpose();
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> grads(
        NumRefShapeFunctions(), 2 * refcoords.cols());
    const auto [basis_dx, basis_dy] = ComputePolyBasisDerivative(refcoords);
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dx =
        ref_func_coeffs_ * basis_dx.transpose();
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dy =
        ref_func_coeffs_ * basis_dy.transpose();
    for (int refcoord_idx = 0; refcoord_idx < refcoords.cols();
         ++refcoord_idx) {
      grads.col(2 * refcoord_idx + 0) = dx.col(refcoord_idx);
      grads.col(2 * refcoord_idx + 1) = dy.col(refcoord_idx);
    }
    return grads;
  }

  /**
   * @brief Evaluation nodes are the vertices, the Chebyshev nodes of degree p-1
   * on the edges and the corresponding nodes on the triangle
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return eval_nodes_;
  }

  /**
   * @brief (p+1)*(p+2)/2 evaluation nodes
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }

 private:
  unsigned degree_;
  Eigen::MatrixXd eval_nodes_;
  Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ref_func_coeffs_;

  Eigen::MatrixXd ComputeEvaluationNodes(unsigned p) const {
    Eigen::MatrixXd eval_nodes(2, (p + 1) * (p + 2) / 2);
    // Add the evaluation nodes corresponding to the vertices of the triangle
    eval_nodes(0, 0) = 0;
    eval_nodes(1, 0) = 0;
    eval_nodes(0, 1) = 1;
    eval_nodes(1, 1) = 0;
    eval_nodes(0, 2) = 0;
    eval_nodes(1, 2) = 1;
    // Add the evaluation nodes corresponding to the edges of the triangle
    const Eigen::VectorXd cheb = chebyshevNodes(p - 1);
    for (int i = 0; i < p - 1; ++i) {
      eval_nodes(0, 3 + i) = cheb[i];
      eval_nodes(1, 3 + i) = 0;
    }
    for (int i = 0; i < p - 1; ++i) {
      eval_nodes(0, 2 + p + i) = 1. - cheb[i];
      eval_nodes(1, 2 + p + i) = cheb[i];
    }
    for (int i = 0; i < p - 1; ++i) {
      eval_nodes(0, 1 + p + p + i) = 0;
      eval_nodes(1, 1 + p + p + i) = 1. - cheb[i];
    }
    // Add the evaluation nodes corresponding to the interior of the triangle
    if (p > 2) {
      for (int i = 0; i < p - 2; ++i) {
        for (int j = 0; j <= i; ++j) {
          eval_nodes(0, (3 * p) + (i * (i + 1) / 2) + j) = cheb[p - 3 - i];
          eval_nodes(1, (3 * p) + (i * (i + 1) / 2) + j) = cheb[j];
        }
      }
    }
    return eval_nodes;
  }

  Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ComputePolyBasis(
      const Eigen::MatrixXd &refcoords) const {
    // The coefficients are ordered x^0y^0, x^0y^1, ..., x^0y^p, x^1y^0, x^1y^1,
    // ..., x^1, y^(p-1), ... x^(p-1)y^0, x^(p-1)y^1, x^py^0
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> coeffs(
        refcoords.cols(), (degree_ + 1) * (degree_ + 2) / 2);
    unsigned coeff_idx = 0;
    for (unsigned powx = 0; powx <= degree_; ++powx) {
      for (unsigned powy = 0; powy <= degree_ - powx; ++powy) {
        for (unsigned node_idx = 0; node_idx < refcoords.cols(); ++node_idx) {
          coeffs(node_idx, coeff_idx) = std::pow(refcoords(0, node_idx), powx) *
                                        std::pow(refcoords(1, node_idx), powy);
        }
        ++coeff_idx;
      }
    }
    return coeffs;
  }

  std::pair<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
            Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>>
  ComputePolyBasisDerivative(const Eigen::MatrixXd &refcoords) const {
    // The coefficients are ordered x^0y^0, x^0y^1, ..., x^0y^p, x^1y^0, x^1y^1,
    // ..., x^1, y^(p-1), ... x^(p-1)y^0, x^(p-1)y^1, x^py^0
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> coeffs_dx(
        refcoords.cols(), (degree_ + 1) * (degree_ + 2) / 2);
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> coeffs_dy(
        refcoords.cols(), (degree_ + 1) * (degree_ + 2) / 2);
    unsigned coeff_idx = 0;
    for (unsigned powx = 0; powx <= degree_; ++powx) {
      for (unsigned powy = 0; powy <= degree_ - powx; ++powy) {
        for (unsigned node_idx = 0; node_idx < refcoords.cols(); ++node_idx) {
          if (powx > 0) {
            coeffs_dx(node_idx, coeff_idx) =
                powx * std::pow(refcoords(0, node_idx), powx - 1) *
                std::pow(refcoords(1, node_idx), powy);
          } else {
            coeffs_dx(node_idx, coeff_idx) = 0;
          }
          if (powy > 0) {
            coeffs_dy(node_idx, coeff_idx) =
                powy * std::pow(refcoords(0, node_idx), powx) *
                std::pow(refcoords(1, node_idx), powy - 1);
          } else {
            coeffs_dy(node_idx, coeff_idx) = 0;
          }
        }
        ++coeff_idx;
      }
    }
    return {coeffs_dx, coeffs_dy};
  }
};

/**
 * @brief Lagrangian Finite Elements of arbitrary degreen on quadrilaterals
 *
 * The evaluation nodes are derived from the Chebyshev interpolation nodes
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeONQuad final
    : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeONQuad(const FeLagrangeONQuad &) = default;
  FeLagrangeONQuad(FeLagrangeONQuad &&) = default;
  FeLagrangeONQuad &operator=(const FeLagrangeONQuad &) = default;
  FeLagrangeONQuad &operator=(FeLagrangeONQuad &&) = default;
  ~FeLagrangeONQuad() = default;

  FeLagrangeONQuad(unsigned degree)
      : degree_(degree), eval_nodes_(), ref_func_coeffs_() {
    eval_nodes_ = ComputeEvaluationNodes(degree);
    ref_func_coeffs_ = ComputePolyBasis(eval_nodes_).inverse().transpose();
  }

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
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim) const override {
    switch (codim) {
      case 0:
        return (degree_ - 1) * (degree_ - 1);
      case 1:
        return degree_ - 1;
      case 2:
        return 1;
      default:
        LF_ASSERT_MSG(false, "Illegal codim " << codim);
    }
  }

  /**
   * @brief One shape function for each vertex, p-1 shape functions on the edges
   * and (p-1)^2 shape functions on the quadrilateral
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(lf::assemble::dim_t)
   */
  [[nodiscard]] lf::base::size_type NumRefShapeFunctions(
      lf::assemble::dim_t codim,
      lf::base::sub_idx_t /*subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    const auto poly_basis = ComputePolyBasis(refcoords);
    return ref_func_coeffs_ * poly_basis.transpose();
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd &refcoords) const override {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> grads(
        NumRefShapeFunctions(), 2 * refcoords.cols());
    const auto [basis_dx, basis_dy] = ComputePolyBasisDerivative(refcoords);
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dx =
        ref_func_coeffs_ * basis_dx.transpose();
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> dy =
        ref_func_coeffs_ * basis_dy.transpose();
    for (int refcoord_idx = 0; refcoord_idx < refcoords.cols();
         ++refcoord_idx) {
      grads.col(2 * refcoord_idx + 0) = dx.col(refcoord_idx);
      grads.col(2 * refcoord_idx + 1) = dy.col(refcoord_idx);
    }
    return grads;
  }

  /**
   * @brief Evaluation nodes are the vertices, the Chebyshev nodes of degree p-1
   * on the edges and the corresponding nodes on the quadrilateral
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return eval_nodes_;
  }

  /**
   * @brief (p+1)^2 evaluation nodes
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] lf::base::size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }

 private:
  unsigned degree_;
  Eigen::MatrixXd eval_nodes_;
  Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ref_func_coeffs_;

  Eigen::MatrixXd ComputeEvaluationNodes(unsigned p) const {
    Eigen::MatrixXd nodes(2, (p + 1) * (p + 1));
    // Add the evaluation nodes corresponding to the vertices
    nodes(0, 0) = 0;
    nodes(1, 0) = 0;
    nodes(0, 1) = 1;
    nodes(1, 1) = 0;
    nodes(0, 2) = 1;
    nodes(1, 2) = 1;
    nodes(0, 3) = 0;
    nodes(1, 3) = 1;
    // Add the evaluation nodes corresponding to the edges of the quad
    const Eigen::VectorXd cheb = chebyshevNodes(p - 1);
    for (int i = 0; i < p - 1; ++i) {
      nodes(0, 4 + i) = cheb[i];
      nodes(1, 4 + i) = 0;
    }
    for (int i = 0; i < p - 1; ++i) {
      nodes(0, 3 + p + i) = 1;
      nodes(1, 3 + p + i) = cheb[i];
    }
    for (int i = 0; i < p - 1; ++i) {
      nodes(0, 2 + p + p + i) = 1. - cheb[i];
      nodes(1, 2 + p + p + i) = 1;
    }
    for (int i = 0; i < p - 1; ++i) {
      nodes(0, 1 + p + p + p + i) = 0;
      nodes(1, 1 + p + p + p + i) = 1. - cheb[i];
    }
    // Add the evaluation nodes corresponding to the interior of the quad
    for (int i = 0; i < p - 1; ++i) {
      for (int j = 0; j < p - 1; ++j) {
        nodes(0, 4 * p + (p - 1) * i + j) = cheb[j];
        nodes(1, 4 * p + (p - 1) * i + j) = cheb[i];
      }
    }
    return nodes;
  }

  Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ComputePolyBasis(
      const Eigen::MatrixXd &refcoords) const {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> poly(
        refcoords.cols(), (degree_ + 1) * (degree_ + 1));
    for (unsigned powx = 0; powx <= degree_; ++powx) {
      for (unsigned powy = 0; powy <= degree_; ++powy) {
        for (unsigned node_idx = 0; node_idx < refcoords.cols(); ++node_idx) {
          poly(node_idx, (degree_ + 1) * powx + powy) =
              std::pow(refcoords(0, node_idx), powx) *
              std::pow(refcoords(1, node_idx), powy);
        }
      }
    }
    return poly;
  }

  std::pair<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
            Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>>
  ComputePolyBasisDerivative(const Eigen::MatrixXd &refcoords) const {
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> coeffs_dx(
        refcoords.cols(), (degree_ + 1) * (degree_ + 1));
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> coeffs_dy(
        refcoords.cols(), (degree_ + 1) * (degree_ + 1));
    for (unsigned powx = 0; powx <= degree_; ++powx) {
      for (unsigned powy = 0; powy <= degree_; ++powy) {
        for (unsigned node_idx = 0; node_idx < refcoords.cols(); ++node_idx) {
          if (powx > 0) {
            coeffs_dx(node_idx, (degree_ + 1) * powx + powy) =
                powx * std::pow(refcoords(0, node_idx), powx - 1) *
                std::pow(refcoords(1, node_idx), powy);
          } else {
            coeffs_dx(node_idx, (degree_ + 1) * powx + powy) = 0;
          }
          if (powy > 0) {
            coeffs_dy(node_idx, (degree_ + 1) * powx + powy) =
                powy * std::pow(refcoords(0, node_idx), powx) *
                std::pow(refcoords(1, node_idx), powy - 1);
          } else {
            coeffs_dy(node_idx, (degree_ + 1) * powx + powy) = 0;
          }
        }
      }
    }
    return {coeffs_dx, coeffs_dy};
  }
};

/**
 * @brief Lagrangian Finite Element Space of arbitrary degree
 */
template <typename SCALAR>
class FeSpaceLagrangeON : public lf::uscalfe::UniformScalarFESpace<SCALAR> {
 public:
  using Scalar = SCALAR;

  FeSpaceLagrangeON() = delete;
  FeSpaceLagrangeON(const FeSpaceLagrangeON &) = delete;
  FeSpaceLagrangeON(FeSpaceLagrangeON &&) noexcept = default;
  FeSpaceLagrangeON &operator=(const FeSpaceLagrangeON &) = delete;
  FeSpaceLagrangeON &operator=(FeSpaceLagrangeON &&) noexcept = default;

  /**
   * @brief Constructor: Sets up the dof handler
   * @param mesh_p A shared pointer to the underlying mesh (immutable)
   * @param N The polynomial degree of the Finite Element Space
   */
  explicit FeSpaceLagrangeON(
      const std::shared_ptr<const lf::mesh::Mesh> &mesh_p, unsigned N)
      : lf::uscalfe::UniformScalarFESpace<SCALAR>(
            mesh_p, std::make_shared<FeLagrangeONTria<SCALAR>>(N),
            std::make_shared<FeLagrangeONQuad<SCALAR>>(N),
            std::make_shared<FeLagrangeONSegment<SCALAR>>(N),
            std::make_shared<lf::uscalfe::FeLagrangePoint<SCALAR>>(N)) {}
  ~FeSpaceLagrangeON() override = default;
};

#endif  // EXAMPLES_LECTUREDEMOS_CONVERGENCESTUDIES_FESPACELAGRANGEON_H_
