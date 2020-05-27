/**
 * @file
 * @brief Data structures representing discontinuous constant finite elments
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#ifndef PROJECTS_DPG_DISCONTINUOUS_FE_CONSTANT_H
#define PROJECTS_DPG_DISCONTINUOUS_FE_CONSTANT_H

#include <lf/uscalfe/uscalfe.h>
#include <lf/fe/fe.h>
#include <typeinfo>
#include "dpg.h"

namespace projects::dpg {

/**
 * @headerfile projects/dpg/lagr_fe_constant.h
 * @brief Discontinuous constant finite element on a triangular reference
 * element.
 *
 * This is a specialization of lf::fe::ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The only reference shape function is constant  and  associated with
 * the barycenter of the reference triangle, it is an interior shape function.
 */
template <typename SCALAR>
class FeDiscontinuousO0Tria final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeDiscontinuousO0Tria(const FeDiscontinuousO0Tria&) = default;
  FeDiscontinuousO0Tria(FeDiscontinuousO0Tria&&) noexcept = default;
  FeDiscontinuousO0Tria& operator=(const FeDiscontinuousO0Tria&) = default;
  FeDiscontinuousO0Tria& operator=(FeDiscontinuousO0Tria&&) noexcept = default;
  FeDiscontinuousO0Tria() = default;
  ~FeDiscontinuousO0Tria() override = default;

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kTria();
  }

  [[nodiscard]] unsigned Degree() const override { return 0; }

  /**
   * @copydoc lf::fe::ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 1; }

  /**   Only one interior shape function
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 1 : 0;
  }
  /**  Only one interior shape function
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 1 : 0;
  }

  // clang-format off
  /** @copydoc lf::fe::ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(1, n_pts);
    result.row(0) = Eigen::RowVectorXd::Constant(n_pts, 1);
    return result;
  }

  // clang-format off
  /** @copydoc lf::fe::ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(1, 2 * n_pts);
    result.row(0) = Eigen::RowVectorXd::Constant(2 * n_pts, 0);
    return result;
  }

  /** @brief Only evaluation node is the barycenter of the reference triangle.
   * @copydoc lf::fe::ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::MatrixXd nodes(2, 1);
    nodes << 1.0 / 3.0, 1.0 / 3.0;
    return nodes;
  }

  /** @brief One evaluation node
   * @copydoc lf::fe::ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

/**
 * @headerfile projects/dpg/lagr_fe_constant.h
 * @brief Discontinuous constant finite element on a quadrilateral reference
 * element.
 *
 * This is a specialization of lf::fe::ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The only reference shape function is constant and associated with
 *  the center of the reference quadrilateral.
 */
template <typename SCALAR>
class FeDiscontinuousO0Quad final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeDiscontinuousO0Quad(const FeDiscontinuousO0Quad&) = default;
  FeDiscontinuousO0Quad(FeDiscontinuousO0Quad&&) noexcept = default;
  FeDiscontinuousO0Quad& operator=(const FeDiscontinuousO0Quad&) = default;
  FeDiscontinuousO0Quad& operator=(FeDiscontinuousO0Quad&&) noexcept = default;
  FeDiscontinuousO0Quad() = default;
  ~FeDiscontinuousO0Quad() override = default;

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kQuad();
  }

  [[nodiscard]] unsigned Degree() const override { return 0; }

  /** @copydoc lf::fe::ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 1; }

  /**   Only one interior shape function
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 1 : 0;
  }
  /**  Only one interior shape function
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 1 : 0;
  }

  // clang-format off
  /** @copydoc lf::fe::ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(1, n_pts);
    result.row(0) = Eigen::RowVectorXd::Constant(n_pts, 1);
    return result;
  }

  // clang-format off
  /** @copydoc lf::fe::ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(1, 2 * n_pts);
    result.row(0) = Eigen::RowVectorXd::Constant(2 * n_pts, 0);
    return result;
  }

  /** @brief Only evaluation node is the center of the reference quadrilateral
   * @copydoc lf::fe::ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::MatrixXd nodes(2, 1);
    nodes << 0.5, 0.5;
    return nodes;
  }
  /** @brief Only one evaluation node
   * @copydoc lf::fe::ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

/**
 * @headerfile projects/dpg/lagr_fe_constant.h
 * @brief Discontinuous constant finite element on a line segment.
 *
 * This is a specialization of lf::fe::ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The only reference shape function is constant  and associated with
 * the  center of the reference line segment
 */
template <typename SCALAR>
class FeDiscontinuousO0Segment final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeDiscontinuousO0Segment(const FeDiscontinuousO0Segment&) = default;
  FeDiscontinuousO0Segment(FeDiscontinuousO0Segment&&) noexcept = default;
  FeDiscontinuousO0Segment& operator=(const FeDiscontinuousO0Segment&) =
      default;
  FeDiscontinuousO0Segment& operator=(FeDiscontinuousO0Segment&&) noexcept =
      default;
  FeDiscontinuousO0Segment() = default;
  ~FeDiscontinuousO0Segment() override = default;

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kSegment();
  }

  [[nodiscard]] unsigned Degree() const override { return 0; }

  /** @copydoc lf::fe::ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 1; }

  /**  Only one interior shape function
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 1 : 0;
  }
  /**  Only one interior shape function
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 1 : 0;
  }

  // clang-format off
  /** @copydoc lf::fe::ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(1, n_pts);
    result.row(0) = Eigen::RowVectorXd::Constant(n_pts, 1);
    return result;
  }

  // clang-format off
  /** @copydoc lf::fe::ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(1, n_pts);
    result.row(0) = Eigen::RowVectorXd::Constant(n_pts, 0);
    return result;
  }

  /** @brief The only evaluation node is the center of the reference line
   * segment.
   * @copydoc lf::fe::ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::MatrixXd nodes(1, 1);
    nodes << 0.5;
    return nodes;
  }
  /** @brief Only one evaluation node.
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_DISCONTINUOUS_FE_CONSTANT_H
