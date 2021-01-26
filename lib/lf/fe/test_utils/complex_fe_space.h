/**
 * @file
 * @brief Provides methods to create complex value FESpaces for testing
 * purposes.
 * @author Raffael Casagrande
 * @date   2021-01-22 05:41:14
 * @copyright MIT License
 */

#ifndef __b9e80d63ee55424493f538a7971df592
#define __b9e80d63ee55424493f538a7971df592

#include <lf/fe/fe.h>
#include <lf/uscalfe/uscalfe.h>

namespace lf::fe::test_utils {
/**
 * @brief Wraps another ScalarReferenceFiniteElement and multiplies the shape
 * functions with the imaginary unit to create complex valued shape functions.
 * @tparam SCALAR Scalar type of the wrapped FiniteElement.
 */
template <class SCALAR>
class ComplexScalarReferenceFiniteElement
    : public lf::fe::ScalarReferenceFiniteElement<std::complex<SCALAR>> {
 public:
  ComplexScalarReferenceFiniteElement(
      std::unique_ptr<lf::fe::ScalarReferenceFiniteElement<SCALAR>> fe)
      : inner_(std::move(fe)) {}

  [[nodiscard]] base::RefEl RefEl() const override { return inner_->RefEl(); }
  [[nodiscard]] unsigned Degree() const override { return inner_->Degree(); }
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t subidx) const override {
    return inner_->NumRefShapeFunctions(codim, subidx);
  }
  [[nodiscard]] Eigen::Matrix<std::complex<SCALAR>, Eigen::Dynamic,
                              Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd &refcoords) const override {
    return std::complex<SCALAR>(0, 1) *
           inner_->EvalReferenceShapeFunctions(refcoords);
  }
  [[nodiscard]] Eigen::Matrix<std::complex<SCALAR>, Eigen::Dynamic,
                              Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd &refcoords) const override {
    return std::complex<SCALAR>(0, 1) *
           inner_->GradientsReferenceShapeFunctions(refcoords);
  }
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return inner_->EvaluationNodes();
  }
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return inner_->NumEvaluationNodes();
  }

  [[nodiscard]] Eigen::Matrix<std::complex<SCALAR>, 1, Eigen::Dynamic>
  NodalValuesToDofs(const Eigen::Matrix<std::complex<SCALAR>, 1, Eigen::Dynamic>
                        &nodvals) const override {
    return inner_->NodalValuesToDofs(
        (nodvals / std::complex<SCALAR>(0, 1)).real());
  }

  [[nodiscard]] size_type NumRefShapeFunctions(dim_t dim) const override {
    return inner_->NumRefShapeFunctions(dim);
  }

 private:
  std::unique_ptr<const lf::fe::ScalarReferenceFiniteElement<SCALAR>> inner_;
};

/**
 * @brief Returns a UniformScalarFESpace that is made up of "complexified" (via
 * ComplexScalarReferenceFiniteElement) Hierarchic order 1 finite elements.
 */
inline std::shared_ptr<uscalfe::UniformScalarFESpace<std::complex<double>>>
MakeComplexLagrangeO1FeSpace(std::shared_ptr<const mesh::Mesh> mesh_p) {
  return std::make_shared<uscalfe::UniformScalarFESpace<std::complex<double>>>(
      mesh_p,
      std::make_shared<ComplexScalarReferenceFiniteElement<double>>(
          std::make_unique<uscalfe::FeLagrangeO1Tria<double>>()),
      std::make_shared<ComplexScalarReferenceFiniteElement<double>>(
          std::make_unique<uscalfe::FeLagrangeO1Quad<double>>()),
      std::make_shared<ComplexScalarReferenceFiniteElement<double>>(
          std::make_unique<uscalfe::FeLagrangeO1Segment<double>>()),
      nullptr);
}
}  // namespace lf::fe::test_utils

#endif  // __b9e80d63ee55424493f538a7971df592
