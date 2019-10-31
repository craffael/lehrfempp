#ifndef PROJECTS_DPG_LAGR_FE_QUADRATIC
#define PROJECTS_DPG_LAGR_FE_QUADRATIC

/**
 * @file
 * @brief Data structures representing quadratic Lagrangian finite elments
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <typeinfo>

#include <lf/uscalfe/uscalfe.h>

#include "dpg.h"
namespace projects::dpg {

/**
 * @headerfile projects/dpg/lagr_fe_quadratic.h
 * @brief Quadratic Lagrangian finite element on a triangular reference element
 *
 * This is a specialization of lf::uscalfe::ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The reference shape functions are combinations of the barycentric coordinate
 * functions on the reference triangle.
 *
 * The first three shape functions are associated to the vertices,
 * the other three shape functions to the edges of the triangle. There are
 * no interior shape functions.
 */
template <typename SCALAR>
class FeLagrangeO2Tria final
    : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO2Tria(const FeLagrangeO2Tria&) = default;
  FeLagrangeO2Tria(FeLagrangeO2Tria&&) noexcept = default;
  FeLagrangeO2Tria& operator=(const FeLagrangeO2Tria&) = default;
  FeLagrangeO2Tria& operator=(FeLagrangeO2Tria&&) noexcept = default;
  FeLagrangeO2Tria() = default;
  ~FeLagrangeO2Tria() override = default;

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kTria();
  }

  [[nodiscard]] unsigned Degree() const override { return 2; }

  /**
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 6; }

  /**   One shape function is attached to each node
   * and one to each edge of the triangle. There are no interior
   * shape functions.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 0 : 1;
  }

  /**  One shape function is attached to each node
   * and one to each edge of the triangle. There are no interior
   * shape functions.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 0 : 1;
  }

  // clang-format off
  /** @copydoc lf::uscalfe::ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(6, n_pts);

    Eigen::ArrayXd x0 = refcoords.row(0).array();
    Eigen::ArrayXd x1 = refcoords.row(1).array();

    result.row(0) = (2.0 * (1 - x0 - x1) * (0.5 - x0 - x1)).matrix();
    result.row(1) = (2.0 * x0 * (x0 - 0.5)).matrix();
    result.row(2) = (2.0 * x1 * (x1 - 0.5)).matrix();
    result.row(3) = (4.0 * (1 - x0 - x1) * x0).matrix();
    result.row(4) = (4.0 * x0 * x1).matrix();
    result.row(5) = (4.0 * (1 - x0 - x1) * x1).matrix();

    return result;
  }

  // clang-format off
  /** @copydoc lf::uscalfe::ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(6, 2 * n_pts);

    // reshape into a 12xn_pts matrix.
    Eigen::Map<Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
               Eigen::AutoAlign>
        temp(&result(0, 0), 12, n_pts);

    Eigen::ArrayXd x0 = refcoords.row(0).array();
    Eigen::ArrayXd x1 = refcoords.row(1).array();

    // d/dx_0
    temp.row(0) = -3.0 + 4.0 * x0 + 4.0 * x1;
    temp.row(1) = 4.0 * x0 - 1.0;
    temp.row(2) = 0.0;
    temp.row(3) = 4 - 0 - 8.0 * x0 - 4.0 * x1;
    temp.row(4) = 4.0 * x1;
    temp.row(5) = -4.0 * x1;

    // d/dx_1
    temp.row(6) = -3.0 + 4.0 * x0 + 4.0 * x1;
    temp.row(7) = 0.0;
    temp.row(8) = 4.0 * x1 - 1.0;
    temp.row(9) = -4 * x0;
    temp.row(10) = 4.0 * x0;
    temp.row(11) = 4.0 - 8.0 * x1 - 4.0 * x0;

    return result;
  }

  /** @brief Evaluation nodes are the vertices of the triangle and
   * the midpoints of the edges.
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    // clang-format off
    Eigen::MatrixXd nodes(2, 6);
    nodes << 0.0, 1.0, 0.0, 0.5, 0.5, 0.0,
             0.0, 0.0, 1.0, 0.0, 0.5, 0.5;
    // clang-format on
    return nodes;
  }

  /** @brief Six evaluation nodes
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

/**
 * @headerfile projects/dpg_lagr_fe_quadratic.h
 * @brief Quadratic Lagrangian finite element on a line segment
 *
 * This is a specialization of lf::uscalfe::ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The first two shape functions are associated with the vertices of the
 * segment. The last shape function is an interior shape function
 */
template <typename SCALAR>
class FeLagrangeO2Segment
    : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO2Segment(const FeLagrangeO2Segment&) = default;
  FeLagrangeO2Segment(FeLagrangeO2Segment&&) noexcept = default;
  FeLagrangeO2Segment& operator=(const FeLagrangeO2Segment&) = default;
  FeLagrangeO2Segment& operator=(FeLagrangeO2Segment&&) noexcept = default;
  FeLagrangeO2Segment() = default;
  ~FeLagrangeO2Segment() override = default;

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kSegment();
  }

  [[nodiscard]] unsigned Degree() const override { return 2; }

  /** @copydoc lf::uscalfe::ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 3; }

  /**  One shape function attached to each node and one interior shape
   * function.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return 1;
  }

  /**  One shape function attached to each node and one interior shape
   * function.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return 1;
  }

  // clang-format off
  /** @copydoc lf::uscalfe::ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        3, refcoords.cols());

    Eigen::ArrayXd x = refcoords.row(0).array();

    // endpoints
    result.row(0) = 2.0 * (1.0 - x) * (0.5 - x);
    result.row(1) = 2.0 * x * (x - 0.5);

    // midpoint
    result.row(2) = 4.0 * (1.0 - x) * x;

    return result;
  }

  // clang-format off
  /** @copydoc lf::uscalfe::ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        3, refcoords.cols());

    Eigen::ArrayXd x = refcoords.row(0).array();
    // endpoints
    result.row(0) = (4.0 * x - 3.0).matrix();
    result.row(1) = (4.0 * x - 1.0).matrix();

    // midpoint:
    result.row(2) = (4.0 - 8.0 * x).matrix();
    return result;
  }

  /** @brief Evaluation nodes are the endpoints of the segment and its midpoint
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::MatrixXd nodes(1, 3);
    nodes << 0.0, 1.0, 0.5;
    return nodes;
  }

  /** @brief Three evaluation nodes
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

/**
 * @headerfile projects/dpg/lagr_fe_quadratic.h
 * @brief Quadratic Lagrangian finite element on a quadrilateral reference
 * element.
 *
 * This is a specialization of lf::uscalfe::ScalarReferenceFiniteElement.
 * Refer to its documentation
 *
 * The shape functions are computed via a tensor product construction:
 * The shape function associated with interpolation node \f$ p =(x_j,y_l)\f$ is
 * computed as
 *
 *  \f[ b^p(x_1,x_2)=b^j(x_1)*b^l(x_2) \f]
 *
 * Where \f$ b^j\f$ and \f$ b^l \f$ are the quadratic Lagrangian shape functions
 * on the reference line segment associated to the points \f$ x_j \f$ and \f$
 * y_l \f$.
 *
 * The first four shape functions are associated with the vertices,
 * the next four shape functions are associated
 * with edges and the last shape function is an interior
 * shape function.
 */
template <typename SCALAR>
class FeLagrangeO2Quad final
    : public lf::uscalfe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO2Quad(const FeLagrangeO2Quad&) = default;
  FeLagrangeO2Quad(FeLagrangeO2Quad&&) noexcept = default;
  FeLagrangeO2Quad& operator=(const FeLagrangeO2Quad&) = default;
  FeLagrangeO2Quad& operator=(FeLagrangeO2Quad&&) noexcept = default;
  FeLagrangeO2Quad() = default;
  ~FeLagrangeO2Quad() override = default;

  [[nodiscard]] lf::base::RefEl RefEl() const override {
    return lf::base::RefEl::kQuad();
  }

  [[nodiscard]] unsigned Degree() const override { return 2; }

  /** @copydoc lf::uscalfe::ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 9; }

  /**   One shape function is attached to each node and edge of
   * the quadrilateral. One interior shape function.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return 1;
  }

  /**  One shape function is attached to each node
   * and each edge of the quadrilateral. One interior shape function.
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return 1;
  }

  // clang-format off
  /** @copydoc lf::uscalfe::ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        9, refcoords.cols());

    // evaluate "marginal" shape functions (b^j and b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> marginal_x0_eval =
        (rsf_segment_.EvalReferenceShapeFunctions(refcoords.row(0))).array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> marginal_x1_eval =
        (rsf_segment_.EvalReferenceShapeFunctions(refcoords.row(1))).array();

    // Evaluate basis functions using the tensor product structure
    for (int i = 0; i < 9; ++i) {
      result.row(i) = (marginal_x0_eval.row(marginal_to_shape_mapping_(i, 0)) *
                       marginal_x1_eval.row(marginal_to_shape_mapping_(i, 1)))
                          .matrix();
    }
    return result;
  }

  // clang-format off
  /** @copydoc lf::uscalfe::ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  // clang-format on
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        9, 2 * refcoords.cols());

    // reshape into a 18xn_pts matrix.
    Eigen::Map<Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
               Eigen::AutoAlign>
        temp(&result(0, 0), 18, refcoords.cols());

    // evaluate "marginal" shape functions (b^j and b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> marginal_x0_eval =
        (rsf_segment_.EvalReferenceShapeFunctions(refcoords.row(0))).array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> marginal_x1_eval =
        (rsf_segment_.EvalReferenceShapeFunctions(refcoords.row(1))).array();

    // evaluate derivatives of "marginal" shape functions (b^j and
    // b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> marginal_x0_grad =
        (rsf_segment_.GradientsReferenceShapeFunctions(refcoords.row(0)))
            .array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> marginal_x1_grad =
        (rsf_segment_.GradientsReferenceShapeFunctions(refcoords.row(1)))
            .array();

    // evaluate gradients using the product rule and
    // the tensor product structure of the basis functions.
    for (int i = 0; i < 9; ++i) {
      // d/dx
      temp.row(i) = (marginal_x0_grad.row(marginal_to_shape_mapping_(i, 0)) *
                     marginal_x1_eval.row(marginal_to_shape_mapping_(i, 1)))
                        .matrix();
      // d/dy
      temp.row(i + 9) =
          (marginal_x1_grad.row(marginal_to_shape_mapping_(i, 1)) *
           marginal_x0_eval.row(marginal_to_shape_mapping_(i, 0)))
              .matrix();
    }
    return result;
  }

  /** @brief Evaluation nodes are the vertices of the quadrilateral,
   * the midpoints of the edges and the center of the  quadrilateral.
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::MatrixXd nodes(2, 9);

    Eigen::MatrixXd vertices(2, 4);
    Eigen::MatrixXd midpoints(2, 4);
    Eigen::MatrixXd center(2, 1);

    // clang-format off
    vertices << 0.0, 1.0, 1.0, 0.0,
                0.0, 0.0, 1.0, 1.0;
    midpoints << 0.5, 1.0, 0.5, 0.0,
                 0.0, 0.5, 1.0, 0.5;
    center << 0.5,
              0.5;
    // clang-format on
    nodes << vertices, midpoints, center;
    return nodes;
  }

  /** @brief Nine evaluation nodes
   * @copydoc lf::uscalfe::ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }

 private:
  /** description of the marginal shape functions on the reference line segment
   * \f$ [0,1] \f$  */
  const static FeLagrangeO2Segment<SCALAR> rsf_segment_;

  /**
   * The shape functions are computed via a tensor product construction:
   * The shape function associated with interpolation node \f$ p =(x_j,y_l)\f$
   * is computed as
   *
   *  \f[ b^p(x_1,x_2)=b^j(x_1)*b^l(x_2) \f]
   *
   * Where \f$ b^j\f$ and \f$ b^l \f$ are the quadratic Lagrangian shape
   * functions on the reference line segment associated to the points \f$ x_j
   * \f$ and \f$ y_l \f$.
   *
   * This matrix of size NumRefShapeFunctions x 2 contains in each row
   * the indices of the "marginal" shape functions, that define
   * the correspoding reference shape function via the above formula.
   * if p was the k-th interpolation node, the k-th entries would contain the
   * indices j and l.
   */
  const static Eigen::MatrixXi marginal_to_shape_mapping_;
};

template <typename SCALAR>
const FeLagrangeO2Segment<SCALAR> FeLagrangeO2Quad<SCALAR>::rsf_segment_ =
    FeLagrangeO2Segment<SCALAR>();

template <typename SCALAR>
const Eigen::MatrixXi FeLagrangeO2Quad<SCALAR>::marginal_to_shape_mapping_ =
    // clang-format off
  (Eigen::MatrixXi(9, 2) <<  0, 0,
                             1, 0,
                             1, 1,
                             0, 1,
                             2, 0,
                             1, 2,
                             2, 1,
                             0, 2,
                             2, 2).finished();
// clang-format on

}  // namespace projects::dpg

#endif  // PROJECTS_DPG_LAGR_FE_QUADRATIC
