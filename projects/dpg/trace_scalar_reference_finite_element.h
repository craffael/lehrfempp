#ifndef PROJECTS_DPG_TRACE_SCALAR_REFERENCE_FINITE_ELEMENT
#define PROJECTS_DPG_TRACE_SCALAR_REFERENCE_FINITE_ELEMENT

/**
 * @file
 * @brief A decorator around a ScalarReferenceFiniteElement, that
 * drops all interior shape functions.
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <iostream>
#include <typeinfo>

#include <lf/uscalfe/uscalfe.h>

#include "dpg.h"

namespace projects::dpg {

/**
 * @headerfile projects/dpg/trace_scalar_reference_finite_element.h
 * @brief Decorator class around a ScalarReferenceFiniteElement to represent
 *  traces on the cell boundary.
 * @tparam SCALAR The scalar type of the shape functions e.g.'double'
 *
 * The class decorates any lf::uscalfe::ScalarReferenceFiniteElement and
 * drops all its interior shape functions.
 *
 *
 * This class is used to represent finite elements discretizing functions in
 *  \f$ H^{1/2}(\mathcal S) \f$, the space of traces of functions in \f$
 * H^1(\Omega) \f$ on the mesh skeleton. Traces are evaluated/integrated only
 * on/over the boundary of cells. If a cardinal basis is used in the underlying
 * finite element, interior shape functions evalute to zero on any point of the
 * cell boundary. To construct a basis corresponding to the traces they thus
 * have to be dropped.
 */
template <typename SCALAR>
class TraceScalarReferenceFiniteElement
    : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  /**
   * @brief Default constructor, does not initialize this class (invalid state).
   *
   * If any method is called upon it, an error is thrown.
   */
  TraceScalarReferenceFiniteElement() = default;

  TraceScalarReferenceFiniteElement(const TraceScalarReferenceFiniteElement&) =
      delete;
  TraceScalarReferenceFiniteElement(
      TraceScalarReferenceFiniteElement&&) noexcept = default;
  TraceScalarReferenceFiniteElement& operator=(
      const TraceScalarReferenceFiniteElement&) = delete;
  TraceScalarReferenceFiniteElement& operator=(
      TraceScalarReferenceFiniteElement&&) noexcept = default;

  explicit TraceScalarReferenceFiniteElement(
      std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>>
          cfe)
      : lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>(),
        cfe_(std::move(cfe)) {}

  /**
   * @brief Report the initialization status of the object
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
    return cfe_->NumRefShapeFunctions() - cfe_->NumRefShapeFunctions(0);
  }

  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return (codim == 0) ? 0 : cfe_->NumRefShapeFunctions(codim);
  }

  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t subidx) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return (codim == 0) ? 0 : cfe_->NumRefShapeFunctions(codim, subidx);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return cfe_->EvalReferenceShapeFunctions(local).topRows(
        NumRefShapeFunctions());
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& local) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return cfe_->GradientsReferenceShapeFunctions(local).topRows(
        NumRefShapeFunctions());
  }

  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return cfe_->EvaluationNodes().leftCols(NumRefShapeFunctions());
  }

  [[nodiscard]] lf::uscalfe::size_type NumEvaluationNodes() const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    return cfe_->NumEvaluationNodes() - cfe_->NumRefShapeFunctions(0);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>& nodvals) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    LF_ASSERT_MSG(nodvals.cols() == NumEvaluationNodes(),
                  "nodvals = " << nodvals << " <-> " << NumEvaluationNodes());

    // the inner element expects more nodal values, than the trace element
    // since the internal shape functions were dropped, so
    // the nodal values are extended by zeroes
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> extendedNodvals(
        cfe_->NumEvaluationNodes());
    Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> zeros(cfe_->NumEvaluationNodes() -
                                                   NumEvaluationNodes());
    zeros.setZero();
    extendedNodvals << nodvals, zeros;
    return cfe_->NodalValuesToDofs(extendedNodvals)
        .leftCols(NumRefShapeFunctions());
  }

  std::ostream& PrintInfo(std::ostream& o, unsigned int = 0) const override {
    LF_ASSERT_MSG(isInitialized(), "Not initialized");
    o << typeid(*this).name() << "trace wrapper around \n";
    cfe_->PrintInfo(o);
    return o;
  }

  /** virtual destructor*/
  ~TraceScalarReferenceFiniteElement() override = default;

 private:
  /** The underlying scalar-valued paramteric finite element */
  std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<SCALAR>> cfe_;
};

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_TRACE_SCALAR_REFERENCE_FINITE_ELEMENT
