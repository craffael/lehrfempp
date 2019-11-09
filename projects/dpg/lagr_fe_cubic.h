#ifndef PROJECTS_DPG_LAGR_FE_CUBIC
#define PROJECTS_DPG_LAGR_FE_CUBIC

/**
 * @file
 * @brief Data structures representing cubic Lagrangian finite elments
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <typeinfo>

#include <lf/uscalfe/uscalfe.h>
#include "dpg.h"

namespace projects::dpg {

/**
 * @headerfile projects/dpg/lagr_fe_cubic.h
 * @brief Cubic Lagrangian finite elment on a triangular reference element
 *
 * This is a specialization of lf::uscalfe::ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The reference shape functions are combinations of the barycentric coordinate
 * functions on the reference triangle.
 *
 * The first three shape functions are associated with vertices,
 * the next six with the  edges of the triangle. The last shape function is an
 * interior shape function.
 *
 * @see lf::uscalfe::ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO3Tria final
    : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO3Tria(const FeLagrangeO3Tria&) = default;
  FeLagrangeO3Tria(FeLagrangeO3Tria&&) noexcept = default;
  FeLagrangeO3Tria& operator=(const FeLagrangeO3Tria&) = default;
  FeLagrangeO3Tria& operator=(FeLagrangeO3Tria&&) noexcept = default;
  FeLagrangeO3Tria() = default;
  ~FeLagrangeO3Tria() override = default;

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kTria();
  }

  [[nodiscard]] unsigned Degree() const override { return 3; }

  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 10; }

  /**
   *  One shape function is attached to each node
   * and two to each edge of the triangle. There is one interior
   * shape functions.
   *
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 1) ? 2 : 1;
  }

  /** One shape function is attached to each node
   * and two to each edge of the triangle. There is one interior
   * shape functions.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 1) ? 2 : 1;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        10, refcoords.cols());

    auto x0 = refcoords.row(0).array();
    auto x1 = refcoords.row(1).array();

    // evaluation of the barycentric coordinate functions:
    auto lambda0 = 1 - x0 - x1;
    auto lambda1 = x0;
    auto lambda2 = x1;

    result.row(0) = 4.5 * lambda0 * (lambda0 - 1 / 3.0) * (lambda0 - 2 / 3.0);
    result.row(1) = 4.5 * lambda1 * (lambda1 - 1 / 3.0) * (lambda1 - 2 / 3.0);
    result.row(2) = 4.5 * lambda2 * (lambda2 - 1 / 3.0) * (lambda2 - 2 / 3.0);

    result.row(3) = 13.5 * lambda0 * lambda1 * (lambda0 - 1 / 3.0);
    result.row(4) = 13.5 * lambda0 * lambda1 * (lambda1 - 1 / 3.0);
    result.row(5) = 13.5 * lambda1 * lambda2 * (lambda1 - 1 / 3.0);
    result.row(6) = 13.5 * lambda1 * lambda2 * (lambda2 - 1 / 3.0);
    result.row(7) = 13.5 * lambda2 * lambda0 * (lambda2 - 1 / 3.0);
    result.row(8) = 13.5 * lambda2 * lambda0 * (lambda0 - 1 / 3.0);

    result.row(9) = 27.0 * lambda0 * lambda1 * lambda2;

    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        10, 2 * refcoords.cols());

    // reshape into a 20xn_pts matrix.
    Eigen::Map<Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
               Eigen::AutoAlign>
        temp(&result(0, 0), 20, n_pts);

    auto x0 = refcoords.row(0).array();
    auto x1 = refcoords.row(1).array();

    // evaulate barycentric coordinate functions:
    auto l0 = 1 - x0 - x1;
    auto l1 = x0;
    auto l2 = x1;

    // vertices
    temp.row(0) = -4.5 * ((l0 - 1 / 3.0) * (l0 - 2 / 3.0) +
                          l0 * (l0 - 2 / 3.0) + l0 * (l0 - 1 / 3.0));
    temp.row(10) = -4.5 * ((l0 - 1 / 3.0) * (l0 - 2 / 3.0) +
                           l0 * (l0 - 2 / 3.0) + l0 * (l0 - 1 / 3.0));
    temp.row(1) = 4.5 * ((l1 - 1 / 3.0) * (l1 - 2 / 3.0) + l1 * (l1 - 2 / 3.0) +
                         l1 * (l1 - 1 / 3.0));
    temp.row(11) = 0.0;
    temp.row(2) = 0.0;
    temp.row(12) = 4.5 * ((l2 - 1 / 3.0) * (l2 - 2 / 3.0) +
                          l2 * (l2 - 2 / 3.0) + l2 * (l2 - 1 / 3.0));

    // edges
    temp.row(3) = 13.5 * (-l1 * (l0 - 1 / 3.0) + l0 * (l0 - 1 / 3.0) - l0 * l1);
    temp.row(13) = -13.5 * (l1 * (l0 - 1 / 3.0) + l0 * l1);
    temp.row(4) = 13.5 * (-l1 * (l1 - 1 / 3.0) + l0 * (l1 - 1 / 3.0) + l0 * l1);
    temp.row(14) = -13.5 * l1 * (l1 - 1 / 3.0);
    temp.row(5) = 13.5 * (l2 * (l1 - 1 / 3.0) + l1 * l2);
    temp.row(15) = 13.5 * (l1 * (l1 - 1 / 3.0));
    temp.row(6) = 13.5 * (l2 * (l2 - 1 / 3.0));
    temp.row(16) = 13.5 * (l1 * (l2 - 1 / 3.0) + l1 * l2);
    temp.row(7) = -13.5 * l2 * (l2 - 1 / 3.0);
    temp.row(17) = 13.5 * (l0 * (l2 - 1 / 3.0) - l2 * (l2 - 1 / 3.0) + l0 * l2);
    temp.row(8) = -13.5 * (l2 * (l0 - 1 / 3.0) + l2 * l0);
    temp.row(18) = 13.5 * (l0 * (l0 - 1 / 3.0) - l2 * (l0 - 1 / 3.0) - l2 * l0);

    // center:
    temp.row(9) = 27.0 * (-l1 * l2 + l0 * l2);
    temp.row(19) = 27.0 * (-l1 * l2 + l0 * l1);

    return result;
  }

  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::MatrixXd nodes(2, 10);
    Eigen::MatrixXd vertices(2, 3);
    Eigen::MatrixXd edges(2, 6);
    Eigen::MatrixXd center(2, 1);

    // clang-format off
    vertices << 0.0, 1.0, 0.0,
                0.0, 0.0, 1.0;
    edges << 1.0, 2.0, 2.0, 1.0, 0.0, 0.0,
             0.0, 0.0, 1.0, 2.0, 2.0, 1.0;
    center << 1.0 / 3.0,
              1.0 / 3.0;
    // clang-format on
    edges *= 1.0 / 3.0;
    nodes << vertices, edges, center;

    return nodes;
  }
  /** @brief Ten evaluation nodes
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

/**
 * @headerfile projects/dpg_lagr_fe_quadratic.h
 * @brief Cubic Lagrangian finite element on a line segment
 *
 * This is a specialization of lf::uscalfe::ScalarReferenceFiniteElement
 * Refer to its documentation.
 *
 * The first two shape functins are associated with vertices of the segment.
 * The last two shape functions are interior shape function
 *
 * @see lf::uscalfe::ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO3Segment final
    : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO3Segment(const FeLagrangeO3Segment&) = default;
  FeLagrangeO3Segment(FeLagrangeO3Segment&&) noexcept = default;
  FeLagrangeO3Segment& operator=(const FeLagrangeO3Segment&) = default;
  FeLagrangeO3Segment& operator=(FeLagrangeO3Segment&&) noexcept = default;
  FeLagrangeO3Segment() = default;
  ~FeLagrangeO3Segment() override = default;

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kSegment();
  }

  [[nodiscard]] unsigned Degree() const override { return 3; }

  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 4; }

  /**  One shape function attached to each node. There are two interior
   * shape functions.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return (codim == 0) ? 2 : 1;
  }

  /**  One shape function attached to each node. There are two interior
   * shape functions.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return (codim == 0) ? 2 : 1;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        4, refcoords.cols());

    auto x = refcoords.row(0).array();

    // vertices
    result.row(0) = 4.5 * (1.0 - x) * (1.0 / 3.0 - x) * (2.0 / 3.0 - x);
    result.row(1) = 4.5 * x * (x - 1.0 / 3.0) * (x - 2.0 / 3.0);

    // interior
    result.row(2) = 13.5 * x * (1 - x) * (2.0 / 3.0 - x);
    result.row(3) = 13.5 * x * (1 - x) * (x - 1.0 / 3.0);

    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        4, refcoords.cols());

    auto x = refcoords.row(0).array();

    // vertex
    result.row(0) = -13.5 * x * x + 18.0 * x - 5.5;
    result.row(1) = 13.5 * x * x - 9.0 * x + 1.0;
    // edge
    result.row(2) = 40.5 * x * x - 45 * x + 9;
    result.row(3) = -40.5 * x * x + 36 * x - 4.5;

    return result;
  }

  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::MatrixXd nodes(1, 4);
    nodes << 0.0, 1.0, 1.0 / 3.0, 2.0 / 3.0;
    return nodes;
  }

  /** @brief Four evaluation nodes
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

/**
 * @headerfile projects/dpg/lagr_fe_quadratic.h
 * @brief Cubic Lagrangian finite element on a quadrilateral reference
 * element.
 *
 * This is a specialization of lf::uscalfe::ScalarReferenceFiniteElement.
 * Refer to its documentation
 *
* The shape functions are computed by a tensor product construction:
 * The shape function associated with the interpolation node \f$ \mathbf{p}
 * =(x_j,y_l)\f$ is computed as
 *
 *  \f[ \hat{b}^{\mathbf{p}}(x_1,x_2)=\hat{b}^j(x_1)*\hat{b}^l(x_2) \f]
 *
 * Where \f$ \hat{b}^j\f$ and \f$ \hat{b}^l \f$ are the cubic Lagrangian
 * shape functions on the reference line segment associated to the interpolation
 * nodes \f$ x_j \f$ and \f$ y_l \f$.
 *
 *
 * The first four shape functions are associated with the vertices,
 * the next eight shape functions are associated
 * with edges of the quadrilateral. The last four shape functions are interior
 * shape function.
 *
 * @see lf::uscalfe::ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO3Quad final
    : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO3Quad(const FeLagrangeO3Quad&) = default;
  FeLagrangeO3Quad(FeLagrangeO3Quad&&) noexcept = default;
  FeLagrangeO3Quad& operator=(const FeLagrangeO3Quad&) = default;
  FeLagrangeO3Quad& operator=(FeLagrangeO3Quad&&) noexcept = default;
  FeLagrangeO3Quad() = default;
  ~FeLagrangeO3Quad() override = default;

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kQuad();
  }

  [[nodiscard]] unsigned Degree() const override { return 3; }

  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 16; }

  /**   One shape function is attached to each node and two to each edge
   * of the quadrilateral. Four interior shape function.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    switch (codim) {
      case 0:
        return 4;
      case 1:
        return 2;
      case 2:
        return 1;
      default:
        LF_ASSERT_MSG(false, "Illegal codim" << codim);
    }
    return 0;
  }
  /**  One shape function is attached to each node and two to each edge of
   * the quadrilateral. Four interior shape function.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /* subidx*/) const override {
    return NumRefShapeFunctions(codim);
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        16, refcoords.cols());

    // evaluate "marginal" shape functions (b^j and b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x0_eval =
        (rsf_segment_.EvalReferenceShapeFunctions(refcoords.row(0))).array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x1_eval =
        (rsf_segment_.EvalReferenceShapeFunctions(refcoords.row(1))).array();

    // Evaluate basis functions using the tensor product structure of the basis
    // functions
    for (int i = 0; i < 16; ++i) {
      result.row(i) = (segment_x0_eval.row(ksegment_to_quad_mapping_(i, 0)) *
                       segment_x1_eval.row(ksegment_to_quad_mapping_(i, 1)))
                          .matrix();
    }

    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        16, 2 * refcoords.cols());

    // reshape into a 18xn_pts matrix.
    Eigen::Map<Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
               Eigen::AutoAlign>
        temp(&result(0, 0), 32, refcoords.cols());

    // evaluate "segment" shape functions (b^j and b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x0_eval =
        (rsf_segment_.EvalReferenceShapeFunctions(refcoords.row(0))).array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x1_eval =
        (rsf_segment_.EvalReferenceShapeFunctions(refcoords.row(1))).array();

    // evaluate derivatives of "segment" shape functions (b^j and
    // b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x0_grad =
        (rsf_segment_.GradientsReferenceShapeFunctions(refcoords.row(0)))
            .array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x1_grad =
        (rsf_segment_.GradientsReferenceShapeFunctions(refcoords.row(1)))
            .array();

    // evaluate gradients using the product rule and
    // the tensor product structure of the basis functions.
    for (int i = 0; i < 16; ++i) {
      temp.row(i) = (segment_x0_grad.row(ksegment_to_quad_mapping_(i, 0)) *
                     segment_x1_eval.row(ksegment_to_quad_mapping_(i, 1)))
                        .matrix();
      temp.row(i + 16) =
          (segment_x1_grad.row(ksegment_to_quad_mapping_(i, 1)) *
           segment_x0_eval.row(ksegment_to_quad_mapping_(i, 0)))
              .matrix();
    }

    return result;
  }

  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::MatrixXd nodes(2, 16);

    Eigen::MatrixXd vertices(2, 4);
    Eigen::MatrixXd midpoints(2, 8);
    Eigen::MatrixXd interior(2, 4);

    // clang-format off
    vertices << 0.0, 1.0, 1.0, 0.0,
               0.0, 0.0, 1.0, 1.0;
    midpoints << 1.0, 2.0, 3.0, 3.0, 2.0, 1.0, 0.0, 0.0,
                 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 2.0, 1.0;
    interior << 1.0, 2.0, 2.0, 1.0,
                1.0, 1.0, 2.0, 2.0;
    // clang-format on
    midpoints *= 1.0 / 3.0;
    interior *= 1.0 / 3.0;

    nodes << vertices, midpoints, interior;

    return nodes;
  }

  /** @brief Sixteen evaluation nodes
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }

 private:
  /** description of the marginal shape functions on the reference line segment
   * \f$  [0,1] \f$  */
  const static FeLagrangeO3Segment<SCALAR> krsf_segment_;

  /**
   * The shape functions are computed by a tensor product construction:
   * The shape function associated with the interpolation node \f$ \mathbf{p}
   * =(x_j,y_l)\f$ is computed as
   *
   *  \f[ \hat{b}^{\mathbf{p}}(x_1,x_2)=\hat{b}^j(x_1)*\hat{b}^l(x_2) \f]
   *
   * Where \f$ \hat{b}^j\f$ and \f$ \hat{b}^l \f$ are the cubic Lagrangian
   * shape functions on the reference line segment associated to the
   * interpolation nodes \f$ x_j \f$ and \f$ y_l \f$.
   *
   * This matrix of size NumRefShapeFunctions x 2 contains in each row
   * the indices of the "segment" shape functions, that define
   * the  reference shape function via the above formula.
   * if \f$p\f$ was the \f$k\f$-th interpolation node, the \f$k\f$-th row would
   * contain the indices \f$j\f$ and \f$l\f$.
   */
  const static Eigen::MatrixXi ksegment_to_quad_mapping_;
};

template <typename SCALAR>
const FeLagrangeO3Segment<SCALAR> FeLagrangeO3Quad<SCALAR>::krsf_segment_ =
    FeLagrangeO3Segment<SCALAR>();

template <typename SCALAR>
const Eigen::MatrixXi FeLagrangeO3Quad<SCALAR>::ksegment_to_quad_mapping_ =
    // clang-format off
  (Eigen::MatrixXi(16,2) << 0, 0,
                            1, 0,
                            1, 1,
                            0, 1,
                            2, 0,
                            3, 0,
                            1, 2,
                            1, 3,
                            3, 1,
                            2, 1,
                            0, 3,
                            0, 2,
                            2, 2,
                            3, 2,
                            3, 3,
                            2, 3).finished();
// clang-format on

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_LAGR_FE_CUBIC
