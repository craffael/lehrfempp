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
#include "fe_space_uniform_scalar.h"

namespace lf::lagrfe {

/**
 * @brief Helper class which stores a ScalarReferenceFiniteElement with
 * precomputed values at the nodes of a quadrature rule.
 * @tparam SCALAR The scalar type of the shape functions.
 *
 * This class does essentially three things:
 * 1. It wraps any ScalarReferenceFiniteElement and forwards all calls to the
 * wrapped instance.
 * 2. It provides access to a quad::QuadratureRule (passed in constructor)
 * 3. It provides additional member functions to access the precomputed values
 * of the shape functions/gradients at the nodes of the quadrature rule.
 */
template <class SCALAR>
class PrecomputedScalarReferenceFiniteElement
    : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  /**
   * @brief Default constructor which does not initialize this class at all
   * (invalid state).
   */
  PrecomputedScalarReferenceFiniteElement() = default;

  PrecomputedScalarReferenceFiniteElement(
      std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>> fe,
      quad::QuadRule qr)
      : ScalarReferenceFiniteElement<SCALAR>(),
        fe_(std::move(fe)),
        qr_(std::move(qr)),
        shap_fun_(fe_->EvalReferenceShapeFunctions(qr_.Points())),
        grad_shape_fun_(fe_->GradientsReferenceShapeFunctions(qr_.Points())) {}

  base::RefEl RefEl() const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->RefEl();
  }

  unsigned Degree() const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->Degree();
  }

  size_type NumRefShapeFunctions() const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->NumRefShapeFunctions();
  }

  size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->NumRefShapeFunctions(codim);
  }

  size_type NumRefShapeFunctions(dim_t codim, sub_idx_t subidx) const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->NumRefShapeFunctions(codim, subidx);
  }

  Eigen::MatrixXd EvalReferenceShapeFunctions(
      const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->EvalReferenceShapeFunctions(local);
  }
  Eigen::MatrixXd GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->GradientsReferenceShapeFunctions(local);
  }
  Eigen::MatrixXd EvaluationNodes() const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->EvaluationNodes();
  }
  size_type NumEvaluationNodes() const override {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return fe_->NumEvaluationNodes();
  }

  Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
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
  const quad::QuadRule& Qr() const {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return qr_;
  }

  /**
   * @brief Value of `EvalReferenceShapeFunctions(Qr().Points())`
   */
  const Eigen::MatrixXd& PrecompReferenceShapeFunctions() const {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return shap_fun_;
  }

  /**
   * @brief Value of `EvalGradientsReferenceShapeFunctions(Qr().Weights())`
   */
  const Eigen::MatrixXd& PrecompGradientsReferenceShapeFunctions() const {
    LF_ASSERT_MSG(fe_ != nullptr, "Not initialized.");
    return grad_shape_fun_;
  }

 private:
  std::shared_ptr<const ScalarReferenceFiniteElement<SCALAR>> fe_;
  quad::QuadRule qr_;
  Eigen::MatrixXd shap_fun_;
  Eigen::MatrixXd grad_shape_fun_;
};

}  // namespace lf::lagrfe

#endif  // __e281cd0ab7fb476e9315a3dda7f45ffe
