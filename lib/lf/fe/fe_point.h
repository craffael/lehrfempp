#ifndef LF_FE_FE_POINT_H_
#define LF_FE_FE_POINT_H_

/**
 * @file
 * @brief Finite Element on vertices for general scalar FE spaces
 * @author Tobias Rohner
 * @date January 2021
 * @copyright MIT License
 */

#include "scalar_reference_finite_element.h"

namespace lf::fe {

/**
 * @headerfile lf/fe/fe.h
 * @brief Linear finite element on a point
 *
 * This is a specialization of ScalarReferenceFiniteElement for an entity
 * of dimension 0, which is exactly one scalar value. It is an ingredient
 * of all scalar finite element spaces.
 */
template <class SCALAR>
class FePoint : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  /**
   * @brief Create a new FePoint by specifying the degree of the shape
   * functions.
   * @param degree The degree of the shape function.
   */
  explicit FePoint(unsigned degree) : degree_(degree) {}

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
    return {0, 1};
  }
  [[nodiscard]] size_type NumEvaluationNodes() const override { return 1; }

 private:
  unsigned degree_;
};

}  // end namespace lf::fe

#endif  // LF_FE_FE_POINT_H_
