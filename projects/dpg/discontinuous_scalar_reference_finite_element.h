#ifndef PROJECTS_DPG_DISCONTINUOUS_SCALAR_REFERENCE_FINITE_ELEMENT
#define PROJECTS_DPG_DISCONTINUOUS_SCALAR_REFERENCE_FINITE_ELEMENT

/**
 * @file
 * @brief A Decorator around a ScalarReferenceFiniteElement to represent
 * discontinuous shape functions
 * @author Philippe Peter
 * @date   June 2019
 * @copyright MIT License
 */

#include <iostream>
#include <typeinfo>

#include <lf/fe/fe.h>
#include <lf/uscalfe/uscalfe.h>

#include "dpg.h"

namespace projects::dpg {

/**
 * @headerfile projecte/dpg/discontinuous_scalar_reference_finite_element.h
 * @brief Decorator class around a  ScalarReferenceFiniteElement to represent
 * discontinuous shape functions.
 * @tparam SCALAR The scalar type of the shape functions e.g. 'double'
 *
 * The class decorates any lf::fe::ScalarReferenceFiniteElement and
 * forwards most calls to the decorated instance. The exception are methods
 * requesting the number of shape functions associated with certain codimensions
 * or subentities. Here the class changes the underlying implementation and
 * associates all shape functions to the underlying entity of codimension 0.
 *
 * In particular this class is used to represent \f$L^2(\Omega) \f$ conforming
 * finite elements. Standard \f$ H^1(\Omega) \f$ Lagrangian finite elements
 * fullfill some continuity constraints, since certain shape functions are
 * associated with vertices or edges of the cells. These continuity constraints
 * are broken, by considering all shape functions as interior.
 */
template <typename SCALAR>
class DiscontinuousScalarReferenceFiniteElement
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  /**
   * Default constructor, does not initialize this class (invalid state).
   * If any method is called upon it, an error is thrown.
   */
  DiscontinuousScalarReferenceFiniteElement() = default;

  DiscontinuousScalarReferenceFiniteElement(
      const DiscontinuousScalarReferenceFiniteElement&) = delete;
  DiscontinuousScalarReferenceFiniteElement(
      DiscontinuousScalarReferenceFiniteElement&&) noexcept = default;
  DiscontinuousScalarReferenceFiniteElement& operator=(
      const DiscontinuousScalarReferenceFiniteElement&) = delete;
  DiscontinuousScalarReferenceFiniteElement& operator=(
      DiscontinuousScalarReferenceFiniteElement&&) noexcept = default;

  explicit DiscontinuousScalarReferenceFiniteElement(
      std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>> cfe)
      : lf::fe::ScalarReferenceFiniteElement<SCALAR>(), cfe_(std::move(cfe)) {}

  /**
   * @brief Reports initialization status of the object
   *
   * Objects built by the default constructor are undefined
   */
  [[nodiscard]] inline bool isInitialized() const { return (cfe_ != nullptr); }

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized!");
    return cfe_->RefEl();
  }

  [[nodiscard]] unsigned Degree() const override {
    LF_ASSERT_MSG(isInitialized(), "Not Initialized");
    return cfe_->Degree();
  }

  [[nodiscard]] size_type NumRefShapeFunctions() const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return cfe_->NumRefShapeFunctions();
  }

  /**  Associates all shape functions to the underlying entity of codim
   * 0*/
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return (codim == 0) ? cfe_->NumRefShapeFunctions() : 0;
  }

  /** Associates all shape functions to the underlying entity of codim
   * 0*/
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return (codim == 0) ? cfe_->NumRefShapeFunctions() : 0;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return cfe_->EvalReferenceShapeFunctions(local);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return cfe_->GradientsReferenceShapeFunctions(local);
  }

  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return cfe_->EvaluationNodes();
  }

  [[nodiscard]] size_type NumEvaluationNodes() const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return cfe_->NumEvaluationNodes();
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>& nodvals) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return cfe_->NodalValuesToDofs(nodvals);
  }

  /** virtual destructor*/
  ~DiscontinuousScalarReferenceFiniteElement() override = default;

 private:
  /** The underlying (continuous) scalar-valued paramteric finite element */
  std::shared_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>> cfe_;
};

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_DISCONTINUOUS_SCALAR_REFERENCE_FINITE_ELEMENT
