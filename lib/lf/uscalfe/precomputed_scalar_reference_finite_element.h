/**
 * @file
 * @brief A wrapper around a ScalarReferenceFiniteElement that additionally
 * precomputes the values of shape function at a set of quadrature points.
 * @author Raffael Casagrande
 * @date   2018-12-29 03:30:59
 * @copyright MIT License
 */

#ifndef __e281cd0ab7fb476e9315a3dda7f45ffe
#define __e281cd0ab7fb476e9315a3dda7f45ffe

#include <Eigen/src/Core/util/ForwardDeclarations.h>
#include <lf/quad/quad.h>
#include <memory>
#include "uniform_scalar_fe_space.h"

namespace lf::uscalfe {

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Helper class which stores a ScalarReferenceFiniteElement with
 * precomputed values at the nodes of a quadrature rule.
 * @tparam SCALAR The scalar type of the shape functions, e.g. `double`
 *
 * This class does essentially three things:
 * -# It wraps any ScalarReferenceFiniteElement and forwards all calls to the
 * wrapped instance.
 * -# It provides access to a quad::QuadratureRule (passed in constructor)
 * -# It provides additional member functions to access the precomputed values
 * -# the shape functions/gradients at the nodes of the quadrature rule.
 *
 * Precomputing entity-independent quantities boost efficiency in the context of
 * parametric finite element methods provided that a uniform quadrature rule is
 * used for local computations on all mesh entities of the same topological
 * type.
 *
 * Detailed explanations can be found in @\lref{par:locparm},
 * @lref{{par:locparm2}.
 */
template <class SCALAR>
class PrecomputedScalarReferenceFiniteElement
    : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  /**
   * @brief Default constructor which does not initialize this class at all
   * (invalid state). If any method is called upon it, an error is thrown.
   */
  PrecomputedScalarReferenceFiniteElement() = default;

  PrecomputedScalarReferenceFiniteElement(
      const PrecomputedScalarReferenceFiniteElement&) = delete;

  PrecomputedScalarReferenceFiniteElement(
      PrecomputedScalarReferenceFiniteElement&&) noexcept = default;

  PrecomputedScalarReferenceFiniteElement& operator=(
      const PrecomputedScalarReferenceFiniteElement&) = delete;

  PrecomputedScalarReferenceFiniteElement& operator=(
      PrecomputedScalarReferenceFiniteElement&&) noexcept = default;

  PrecomputedScalarReferenceFiniteElement(
      std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>> fe,
      quad::QuadRule qr)
      : ScalarReferenceFiniteElement<SCALAR>(),
        fe_(std::move(fe)),
        qr_(std::move(qr)),
        shap_fun_(fe_->EvalReferenceShapeFunctions(qr_.Points())),
        grad_shape_fun_(fe_->GradientsReferenceShapeFunctions(qr_.Points())) {}

  /**
   * @brief Tells initialization status of object
   *
   * An object is in an undefined state when built by the default constructor
   */
  [[nodiscard]] bool isInitialized() const { return (fe_ != nullptr); }

  [[nodiscard]] base::RefEl RefEl() const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->RefEl();
  }

  [[nodiscard]] unsigned Degree() const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->Degree();
  }

  [[nodiscard]] size_type NumRefShapeFunctions() const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->NumRefShapeFunctions();
  }

  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->NumRefShapeFunctions(codim);
  }

  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t subidx) const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->NumRefShapeFunctions(codim, subidx);
  }

  [[nodiscard]] Eigen::MatrixXd EvalReferenceShapeFunctions(
      const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->EvalReferenceShapeFunctions(local);
  }
  [[nodiscard]] Eigen::MatrixXd GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->GradientsReferenceShapeFunctions(local);
  }
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->EvaluationNodes();
  }
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->NumEvaluationNodes();
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>& nodvals) const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->NodalValuesToDofs(nodvals);
  }

  std::ostream& print(std::ostream& o) const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->print(o);
  }

  /**
   * @brief Return the Quadrature rule at which the shape functions (and their
   * gradients) have been precomputed.
   */
  [[nodiscard]] const quad::QuadRule& Qr() const {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return qr_;
  }

  /**
   * @brief Value of `EvalReferenceShapeFunctions(Qr().Points())`
   */
  [[nodiscard]] const Eigen::MatrixXd& PrecompReferenceShapeFunctions() const {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return shap_fun_;
  }

  /**
   * @brief Value of `EvalGradientsReferenceShapeFunctions(Qr().Weights())`
   *
   * See @ref ScalarReferenceFiniteElement::EvalGradientsReferenceShapeFunctions
   * for the packed format in which the gradients are returned.
   */
  [[nodiscard]] const Eigen::MatrixXd& PrecompGradientsReferenceShapeFunctions()
      const {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return grad_shape_fun_;
  }

  /** virtual destructor */
  ~PrecomputedScalarReferenceFiniteElement() override = default;

 private:
  /** The underlying scalar-valued parametric finite element */
  std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>> fe_;
  /** Uniform parametric quadrature rule for the associated type of reference
   * element */
  quad::QuadRule qr_;
  /** Holds values of reference shape functions at reference quadrature nodes */
  Eigen::MatrixXd shap_fun_;
  /** Holds gradients of reference shape functions at quadrature nodes */
  Eigen::MatrixXd grad_shape_fun_;
};

}  // namespace lf::uscalfe

#endif  // __e281cd0ab7fb476e9315a3dda7f45ffe
