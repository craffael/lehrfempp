#ifndef LF_FE
#define LF_FE
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementati0on of parametric Lagrangian finite elements
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include <lf/assemble/dofhandler.h>
#include <lf/fe/scalar_reference_finite_element.h>
#include <iostream>
#include <typeinfo>

namespace lf::uscalfe {
/** Type for indices into global matrices/vectors */
using gdof_idx_t = lf::assemble::gdof_idx_t;
/** Type for indices referring to entity matrices/vectors */
using ldof_idx_t = lf::assemble::ldof_idx_t;
/** Type for vector length/matrix sizes */
using size_type = lf::assemble::size_type;
/** Type for (co-)dimensions */
using dim_t = lf::assemble::dim_t;
/** Type for global index of entities */
using glb_idx_t = lf::assemble::glb_idx_t;
/** Type for indexing sub-entities */
using sub_idx_t = lf::base::sub_idx_t;

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Linear Lagrange finite element on triangular reference element
 *
 * This is a specialization of ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The reference shape functions are
   \f[
   \widehat{b}^1(\widehat{\mathbf{x}}) =
 1-\widehat{x}_1-\widehat{x}_2\quad,\quad \widehat{b}^2(\widehat{\mathbf{x}}) =
 \widehat{x}_1\quad,\quad \widehat{b}^3(\widehat{\mathbf{x}}) = \widehat{x}_2
 \;. \f]
 * The basis function \f$\widehat{b}^i\f$ is associated with vertex \f$i\f$
 * of the triangle.
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO1Tria final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO1Tria(const FeLagrangeO1Tria&) = default;
  FeLagrangeO1Tria(FeLagrangeO1Tria&&) noexcept = default;
  FeLagrangeO1Tria& operator=(const FeLagrangeO1Tria&) = default;
  FeLagrangeO1Tria& operator=(FeLagrangeO1Tria&&) noexcept = default;
  FeLagrangeO1Tria() = default;
  ~FeLagrangeO1Tria() override = default;

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kTria();
  }

  [[nodiscard]] unsigned Degree() const override { return 1; }

  /** @brief The local shape functions: barycentric coordinate functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 3; }

  /** @brief Exactly one shape function attached to each node, none to other
   * sub-entities
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim " << codim);
    return (codim == 2) ? 1 : 0;
  }
  /** @brief Exactly one shape function attached to each node, none to other
   * sub-entities
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim " << codim);
    return (codim == 2) ? 1 : 0;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        3, refcoords.cols());
    result.row(0) = Eigen::RowVectorXd::Ones(refcoords.cols()) -
                    refcoords.row(0) - refcoords.row(1);
    result.block(1, 0, 2, refcoords.cols()) = refcoords;
    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const size_type n_pts(refcoords.cols());

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        3, 2 * refcoords.cols());
    result.row(0) = Eigen::RowVectorXd::Constant(2 * n_pts, -1);
    result.row(1) = Eigen::RowVector2d(1., 0.).replicate(1, n_pts);
    result.row(2) = Eigen::RowVector2d(0., 1.).replicate(1, n_pts);
    return result;
  }

  /** @brief Evalutation nodes are just the vertices of the triangle
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return RefEl().NodeCoords();
  }

  /** @brief Three evaluation nodes
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return RefEl().NumNodes();
  }
};

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Linear Lagrange finite element on the quadrilateral reference element
 *
 * This is a specialization of ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The reference shape functions are
 * @f[
       \hat{b}^1(x_1,x_2) = (1-x_1)(1-x_2)\;,\quad
       \hat{b}^2(x_1,x_2) = x_1(1-x_2)\;,\quad
       \hat{b}^3(x_1,x_2) = x_1\cdot x_2\;,\quad
       \hat{b}^4(x_1,x_2) = (1-x_1)x_2\;.
 * @f]
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO1Quad final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO1Quad(const FeLagrangeO1Quad&) = default;
  FeLagrangeO1Quad(FeLagrangeO1Quad&&) noexcept = default;
  FeLagrangeO1Quad& operator=(const FeLagrangeO1Quad&) = default;
  FeLagrangeO1Quad& operator=(FeLagrangeO1Quad&&) noexcept = default;
  FeLagrangeO1Quad() = default;
  ~FeLagrangeO1Quad() override = default;

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kQuad();
  }

  [[nodiscard]] unsigned Degree() const override { return 1; }

  /** @brief The local shape functions: four bilinear basis functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 4; }

  /** @brief Exactly one shape function attached to each node, none to other
   * sub-entities
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    return (codim == 2) ? 1 : 0;
  }
  /** @brief Exactly one shape function attached to each node, none to other
   * sub-entities
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    return (codim == 2) ? 1 : 0;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        4, refcoords.cols());

    result.row(0) =
        ((1 - refcoords.row(0).array()) * (1 - refcoords.row(1).array()))
            .matrix();
    result.row(1) =
        ((refcoords.row(0).array()) * (1 - refcoords.row(1).array())).matrix();
    result.row(2) =
        ((refcoords.row(0).array()) * (refcoords.row(1).array())).matrix();
    result.row(3) =
        ((1 - refcoords.row(0).array()) * (refcoords.row(1).array())).matrix();

    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    size_type n_pts(refcoords.cols());

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(4, 2 * n_pts);

    // reshape the result matrix into a 8xn_pts matrix
    Eigen::Map<Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
               Eigen::AutoAlign>
        temp(&result(0, 0), 8, n_pts);
    temp.row(0) = refcoords.row(1).array() - 1.0;
    temp.row(1) = 1.0 - refcoords.row(1).array();
    temp.row(2) = refcoords.row(1).array();
    temp.row(3) = -refcoords.row(1).array();
    temp.row(4) = refcoords.row(0).array() - 1.0;
    temp.row(5) = -refcoords.row(0).array();
    temp.row(6) = refcoords.row(0).array();
    temp.row(7) = 1.0 - refcoords.row(0).array();

    return result;
  }

  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return RefEl().NodeCoords();
  }

  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return RefEl().NumNodes();
  }
};

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Linear Lagrange finite element on a line segment
 *
 * This is a specialization of ScalarReferenceFiniteElement for an entity
 * of dimension 1! The reference element is the unit interval \f$[0,1]\f$.
 *
 * The reference shape functions are
   \f[
   \widehat{b}^1(\widehat{{x}}) = 1-\widehat{x}\quad,\quad
   \widehat{b}^2(\widehat{{x}}) = \widehat{x}\;.
   \f]
 * The basis function \f$\widehat{b}^i\f$ is associated with vertex \f$i\f$
 * of the edge.
 *
 *  @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO1Segment final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO1Segment(const FeLagrangeO1Segment&) = default;
  FeLagrangeO1Segment(FeLagrangeO1Segment&&) noexcept = default;
  FeLagrangeO1Segment& operator=(const FeLagrangeO1Segment&) = default;
  FeLagrangeO1Segment& operator=(FeLagrangeO1Segment&&) noexcept = default;
  FeLagrangeO1Segment() = default;
  ~FeLagrangeO1Segment() override = default;

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kSegment();
  }

  [[nodiscard]] unsigned Degree() const override { return 1; }

  /** @brief The local shape functions: barycentric coordinate functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 2; }

  /** @brief All shape functions attached to nodes
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return (codim == 1) ? 1 : 0;
  }
  /** @brief All shape functions attached to nodes
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return (codim == 1) ? 1 : 0;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        2, refcoords.cols());
    result.row(0) = Eigen::RowVectorXd::Constant(1, n_pts, 1.0) - refcoords;
    result.row(1) = refcoords;

    return result;
  }

  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(
        2, refcoords.cols());
    result.row(0) = Eigen::MatrixXd::Constant(1, n_pts, -1.0);
    result.row(1) = Eigen::MatrixXd::Constant(1, n_pts, 1.0);

    return result;
  }

  /** @brief Evalutation nodes are just the vertices of the segment
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return RefEl().NodeCoords();
  }

  /** @brief Two evaluation nodes
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return RefEl().NumNodes();
  }
};

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Linear Lagrange finite element on a point
 *
 * This is a specialization of ScalarReferenceFiniteElement for an entity
 * of dimension 0, which is exactly one scalar value. It is an ingredient
 * of all Lagrange type finite element spaces (any degree).
 */
template <class SCALAR>
class FeLagrangePoint : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  /**
   * @brief Create a new FeLagrangePoint by specifying the degree of the shape
   * functions.
   * @param degree The degree of the shape function.
   */
  explicit FeLagrangePoint(unsigned degree) : degree_(degree) {}

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
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 0, "refcoords has too many rows.");
    return Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>::Ones(1, refcoords.cols());
  }
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& /*refcoords*/) const override {
    LF_VERIFY_MSG(false, "gradients not defined in points of mesh.");
  }
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    return Eigen::MatrixXd(0, 1);
  }
  [[nodiscard]] size_type NumEvaluationNodes() const override { return 1; }

 private:
  unsigned degree_;
};

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Quadratic Lagrangian finite element on a triangular reference element
 *
 * This is a specialization of @ref ScalarReferenceFiniteElement.
 *
 * The reference shape functions are combinations of the barycentric coordinate
 * functions on the reference triangle, see @lref{ex:quadLFE}.
 *
 * The first three shape functions are associated to the vertices,
 * the other three shape functions to the edges of the triangle. There are
 * no interior shape functions, see @lref{eq:quadlsf}:
 *
 * \f{eqnarray}{
    \hat{b}^1({\bf x}) &=& 2(1-x_1-x_2)(\frac12 -x_1 -x_2)\;,\\
    \hat{b}^2({\bf x}) &=& 2x_1(x_1-\frac12)\;,\\
    \hat{b}^3({\bf x}) &=& 2x_2(x_2-\frac12)\;,\\
    \hat{b}^4({\bf x}) &=& 4(1-x_1-x_2)x_1\;,\\
    \hat{b}^5({\bf x}) &=& 4x_1x_2\;,\\
    \hat{b}^6({\bf x}) &=& 4(1-x_1-x_2)x_2\;.
 * \f}
 *
 * The numbering convention for the local shape functions is explained in
 @lref{ex:quadnodes}.
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO2Tria final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO2Tria(const FeLagrangeO2Tria&) = default;
  FeLagrangeO2Tria(FeLagrangeO2Tria&&) noexcept = default;
  FeLagrangeO2Tria& operator=(const FeLagrangeO2Tria&) = default;
  FeLagrangeO2Tria& operator=(FeLagrangeO2Tria&&) noexcept = default;
  FeLagrangeO2Tria() = default;
  ~FeLagrangeO2Tria() override = default;

  /** @copydoc ScalarReferenceFiniteElement::RefEl()
   */
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kTria();
  }

  /**  @brief Quadratic Lagrangian finite elements rely on polynomials of degree
   * 2
   */
  [[nodiscard]] unsigned Degree() const override { return 2; }

  /** @brief Six local shape functions are associated with a triangular cell in
   *         the case of quadratic Lagrangian finite elements.
   *
   * @sa ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 6; }

  /**  @brief One shape function attached to each node
   * and one to each edge of the triangle. There are no interior
   * shape functions.
   *
   * @sa ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 0 : 1;
  }

  /**  @brief One shape function attached to each node
   * and one to each edge of the triangle. There are no interior
   * shape functions.
   *
   * @sa ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 0 : 1;
  }

  /** @brief Point evaluation of reference shape functions
   *
   * The formulas of the reference shape functions are given in the class
   * documentation of @ref FeLagrangeO2Tria.
   *
   * @sa ScalarReferenceFiniteElement::EvalReferenceShapeFunctions()
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    // Number of evaluation points
    const size_type n_pts(refcoords.cols());
    // Result returned in 6 x P matrix
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(6, n_pts);
    // Convert to array type for componentwise operations
    auto x0 = refcoords.row(0).array();
    auto x1 = refcoords.row(1).array();
    // Compound evaluation of formulas for reference shape functions
    result.row(0) = (2.0 * (1 - x0 - x1) * (0.5 - x0 - x1)).matrix();
    result.row(1) = (2.0 * x0 * (x0 - 0.5)).matrix();
    result.row(2) = (2.0 * x1 * (x1 - 0.5)).matrix();
    result.row(3) = (4.0 * (1 - x0 - x1) * x0).matrix();
    result.row(4) = (4.0 * x0 * x1).matrix();
    result.row(5) = (4.0 * (1 - x0 - x1) * x1).matrix();
    return result;
  }

  /** @brief Point evaluations of gradients of reference shape functions
   *
   *  Refer to @lref{ex:qfemelmat} for related explanations
   *
   * @sa ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions()
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    // Number of evaluation points
    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(6, 2 * n_pts);
    // reshape into a 12xn_pts matrix.
    Eigen::Map<Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
               Eigen::AutoAlign>
        temp(result.data(), 12, n_pts);
    // Convert to array type for componentwise operations
    auto x0 = refcoords.row(0).array();
    auto x1 = refcoords.row(1).array();
    // d/dx_0 for all reference shape functions
    temp.row(0) = -3.0 + 4.0 * x0 + 4.0 * x1;
    temp.row(1) = 4.0 * x0 - 1.0;
    temp.row(2) = 0.0;
    temp.row(3) = 4 - 0 - 8.0 * x0 - 4.0 * x1;
    temp.row(4) = 4.0 * x1;
    temp.row(5) = -4.0 * x1;
    // d/dx_1 for the reference shape functions
    temp.row(6) = -3.0 + 4.0 * x0 + 4.0 * x1;
    temp.row(7) = 0.0;
    temp.row(8) = 4.0 * x1 - 1.0;
    temp.row(9) = -4 * x0;
    temp.row(10) = 4.0 * x0;
    temp.row(11) = 4.0 - 8.0 * x1 - 4.0 * x0;
    return result;
  }

  /** @brief Evaluation nodes are the three vertices of the triangle and
   * the three midpoints of its edges.
   *
   * The reference shape functions as implemented by this class are a _cardinal
   * basis_ of the space of quadratic two-variate polynomials  with respect to
   * these evaluation nodes.
   *
   * The numbering of evluation nodes is fixed by the numbering of the local
   * shape functions.
   *
   * @sa ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    // clang-format off
    Eigen::Matrix<double, 2,6> nodes;
    // Reference coordinates of evaluation nodes in the columns of a 2 x 6 - matrix. 
    nodes << 0.0, 1.0, 0.0,  0.5, 0.5, 0.0,
             0.0, 0.0, 1.0,  0.0, 0.5, 0.5;
    // clang-format on
    return nodes;
  }

  /** @brief Six evaluation nodes, the same number as local shape functions
   *
   * @sa ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Quadratic Lagrangian finite element on a line segment
 *
 * This is a specialization of ScalarReferenceFiniteElement for quadratic
 * Lagrangian finite elements on a segment.
 *
 * The first two shape functions are associated with the vertices of the
 * segment. The last one is an interior shape function. This complies with the
 * ordering convention of LehreFEM++ for local shape functions, see
 * @lref{par:betlordlsf}.
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO2Segment
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO2Segment(const FeLagrangeO2Segment&) = default;
  FeLagrangeO2Segment(FeLagrangeO2Segment&&) noexcept = default;
  FeLagrangeO2Segment& operator=(const FeLagrangeO2Segment&) = default;
  FeLagrangeO2Segment& operator=(FeLagrangeO2Segment&&) noexcept = default;
  FeLagrangeO2Segment() = default;
  ~FeLagrangeO2Segment() override = default;

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kSegment();
  }

  /** @brief Quadratic Lagrangian finite element rely on polynomials of degree 2
   */
  [[nodiscard]] unsigned Degree() const override { return 2; }

  /** @brief Three local shape functions are associated with an edge
   *
   * @sa ScalarReferenceFiniteElement::NumRefShapeFunctions()
   * @sa FeLagrangeO2Segment::EvalReferenceShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 3; }

  /** @brief One shape function attached to each node and one interior shape
   * function.
   *
   * Hence (sub-)entities of any co-dimension possess exactly one local shape
   * function.
   *
   * @sa ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return 1;
  }

  /** @brief One shape function attached to each node and one interior shape
   * function.
   *
   * @sa ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return 1;
  }

  /** @brief Evaluation of shape functions on reference segment (unit interval)
   *
   * The local shape functions are associated with a segment in the case of
   * quadratic Lagrangian finite elements. In reference coordinates those are:
   * \f[
         \hat{b}^1(x) = 2(1-x)(\frac12-x)\;,\quad
         \hat{b}^2(x) = 2x(x-\frac12)\;,\quad
         \hat{b}^3(x) = 4(1-x)x\;.
   * \f]
   * This is a _cardinal basis_ of the space of quadratic univariate polynomials
   with respect to the points \f$x=0,\frac12,1\f$.
   *
   * @sa ScalarReferenceFiniteElement::EvalReferenceShapeFunctions()
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());
    // Matrix for returning values of local shape functions
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(3, n_pts);
    // Convert into an array type for componentwise operations
    auto x = refcoords.row(0).array();
    // local shape functions belonging to the endpoints
    result.row(0) = 2.0 * (1.0 - x) * (0.5 - x);
    result.row(1) = 2.0 * x * (x - 0.5);
    // local shape function sitting at midpoint (belonging to segment)
    result.row(2) = 4.0 * (1.0 - x) * x;
    return result;
  }

  /** @brief Evaluation of derivatives of local shape function on reference
   * interval
   *
   * For the passed points return the derivative of the three local shape
   * functions on the reference interval in the rows of a matrix.
   *
   * @sa FeLagrangeO2Segment::EvalReferenceShapeFunctions()
   * @sa ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions()
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());
    // Matrix for returning the derivatives of the local shape functions
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(3, n_pts);
    // Convert to array for componentwise operations
    auto x = refcoords.row(0).array();
    // LSF at endpoints
    result.row(0) = (4.0 * x - 3.0).matrix();
    result.row(1) = (4.0 * x - 1.0).matrix();
    // LSF at midpoint:
    result.row(2) = (4.0 - 8.0 * x).matrix();
    return result;
  }

  /** @brief Evaluation nodes are the endpoints of the segment and its midpoint
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::Matrix<double, 1, 3> nodes;
    // Reference coordinates of the three evaluation nodes
    nodes << 0.0, 1.0, 0.5;
    return nodes;
  }

  /** @brief Three evaluation nodes, same number as local shape functions
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Quadratic Lagrangian finite element on a quadrilateral reference
 * element.
 *
 * This is a specialization of ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The shape functions are computed by a tensor product construction:
 * The shape function associated with the interpolation node \f$ \mathbf{p}
 * =(x_j,y_l)\f$ is computed as
 *
 *  \f[ \hat{b}^{\mathbf{p}}(x_1,x_2)=\hat{b}^j(x_1)*\hat{b}^l(x_2) \f]
 *
 * Where \f$ \hat{b}^j\f$ and \f$ \hat{b}^l \f$ are the quadratic Lagrangian
 * shape functions on the reference line segment associated to the interpolation
 * nodes \f$ x_j \f$ and \f$ y_l \f$.
 *
 * The first four shape functions are associated with the vertices,
 * the next four with the edges of the quadrilateral. The last shape function is
 * an interior shape function.
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO2Quad final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO2Quad(const FeLagrangeO2Quad&) = default;
  FeLagrangeO2Quad(FeLagrangeO2Quad&&) noexcept = default;
  FeLagrangeO2Quad& operator=(const FeLagrangeO2Quad&) = default;
  FeLagrangeO2Quad& operator=(FeLagrangeO2Quad&&) noexcept = default;
  FeLagrangeO2Quad() = default;
  ~FeLagrangeO2Quad() override = default;

  /** @copydoc ScalarReferenceFiniteElement::RefEl() */
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kQuad();
  }

  /** @brief Quadratic Lagrangian finite elements sport polynomials of degree 2
   */
  [[nodiscard]] unsigned Degree() const override { return 2; }

  /** @brief Nine local shape functions are associated with a quadrilateral cell
   *
   * Refer to @lref{ex:QuadTPLFE}.
   *
   * @sa ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 9; }

  /** @brief One shape function is attached to each node and edge of
   * the quadrilateral. One interior shape function.
   *
   * @sa ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t codim)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return 1;
  }

  /** @brief One shape function is attached to each node
   * and each edge of the quadrilateral. One interior shape function.
   *
   * @sa ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,
   * sub_idx_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return 1;
  }

  /** @brief Point evaluation of reference local shape functions in unit square
   *
   * The nine local shape functions for quadratic Lagrangian finite elements on
   * the unit square can be obtained by forming tensor product of the three
   * local shape functions on the unit interval. This entails taking into
   * account the numbering of the local shape functions: LSF 0-3 belong to the
   * four corners, LSF 4-7 belong to the edges, LSF 8 belongs to the cell. This
   * complies with the numbering convention for local shape functions as
   * explained in @lref{par:betlordlsf}, @lref{ex:quadnodes}.
   *
   * Here the local shape functions are chosen as _cardinal basis_ of the space
   * of tensor product polynomials of degree 2 with respect to the following
   * nine evaluation nodes: the four corners of the unit square, the four
   * midpoints of its edges, and its center of gravity, see @lref{tplfeldof}.
   *
   * @see ScalarReferenceFiniteElement::EvalReferenceShapeFunctions()
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    // Number of evaluation points
    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(9, n_pts);

    // evaluate "segment" shape functions (b^j and b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x0_eval =
        (krsf_segment_.EvalReferenceShapeFunctions(refcoords.row(0))).array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x1_eval =
        (krsf_segment_.EvalReferenceShapeFunctions(refcoords.row(1))).array();

    // Evaluate basis functions using the tensor product structure
    for (int i = 0; i < 9; ++i) {
      result.row(i) = (segment_x0_eval.row(ksegment_to_quad_mapping_(i, 0)) *
                       segment_x1_eval.row(ksegment_to_quad_mapping_(i, 1)))
                          .matrix();
    }
    return result;
  }

  /** @brief Point evaluation of gradient of reference shape functions for
   * quadratic Lagrangian finite elements on quadrilateral cells
   *
   * Gradients are returned packed in the rows of a 9 x 2P matrix, where
   * P is the number of evaluation points passed as a 2 x P matrix.
   *
   * @sa FeLagrangeO2Quad::EvalReferenceShapeFunctions()
   * @sa ScalarReferenceFiniteElement::EvalReferenceShapeFunctions()
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(9, 2 * n_pts);

    // reshape into a 18xn_pts matrix.
    Eigen::Map<Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
               Eigen::AutoAlign>
        temp(&result(0, 0), 18, n_pts);

    // evaluate "segment" shape functions (b^j and b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x0_eval =
        (krsf_segment_.EvalReferenceShapeFunctions(refcoords.row(0))).array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x1_eval =
        (krsf_segment_.EvalReferenceShapeFunctions(refcoords.row(1))).array();

    // evaluate derivatives of "segment" shape functions (b^j and
    // b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x0_grad =
        (krsf_segment_.GradientsReferenceShapeFunctions(refcoords.row(0)))
            .array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x1_grad =
        (krsf_segment_.GradientsReferenceShapeFunctions(refcoords.row(1)))
            .array();

    // evaluate gradients  using the product rule and
    // the tensor product structure of the basis functions.
    for (int i = 0; i < 9; ++i) {
      // d/dx
      temp.row(i) = (segment_x0_grad.row(ksegment_to_quad_mapping_(i, 0)) *
                     segment_x1_eval.row(ksegment_to_quad_mapping_(i, 1)))
                        .matrix();
      // d/dy
      temp.row(i + 9) = (segment_x1_grad.row(ksegment_to_quad_mapping_(i, 1)) *
                         segment_x0_eval.row(ksegment_to_quad_mapping_(i, 0)))
                            .matrix();
    }
    return result;
  }

  /** @brief Evaluation nodes are the vertices,
   * the midpoints of the edges and the center of the  quadrilateral.
   *
   * Location and numbering of evaluation nodes:
   * @image html lnlfe2.png
   *
   * @return a 2x9 matrix containing the refrence coordinates of the evaluation
   * nodes in its columns
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::Matrix<double, 2, 9> nodes;

    // clang-format off
    nodes << 0.0, 1.0, 1.0, 0.0,  0.5, 1.0, 0.5, 0.0,  0.5,
               0.0, 0.0, 1.0, 1.0,  0.0, 0.5, 1.0, 0.5,  0.5;
    // clang-format on
    return nodes;
  }

  /** @brief Nine evaluation nodes, same as the number of local shape functions
   *
   * @sa FeLagrangeO2Quad::EvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }

 private:
  /** description of the shape functions on the reference line segment
   * \f$ [0,1] \f$  that are used in the tensor product construction of the
   * shape functions*/
  const static FeLagrangeO2Segment<SCALAR> krsf_segment_;

  /**
   * The shape functions are computed by a tensor product construction:
   * The shape function associated with the interpolation node \f$ \mathbf{p}
   * =(x_j,y_l)\f$ is computed as
   *
   *  \f[ \hat{b}^{\mathbf{p}}(x_1,x_2)=\hat{b}^j(x_1)*\hat{b}^l(x_2) \f]
   *
   * Where \f$ \hat{b}^j\f$ and \f$ \hat{b}^l \f$ are the quadratic Lagrangian
   * shape functions on the reference line segment associated to the
   * interpolation nodes \f$ x_j \f$ and \f$ y_l \f$.
   *
   * This matrix of size NumRefShapeFunctions x 2 contains in each row
   * the indices of the "segment" shape functions, that define
   * the  reference shape function via the above formula.
   * if \f$p\f$ was the \f$k\f$-th interpolation node, the \f$k\f$-th row would
   * contain the indices \f$j\f$ and \f$l\f$.
   */
  const static Eigen::Matrix<int, 9, 2> ksegment_to_quad_mapping_;
};

template <typename SCALAR>
const FeLagrangeO2Segment<SCALAR> FeLagrangeO2Quad<SCALAR>::krsf_segment_ =
    FeLagrangeO2Segment<SCALAR>();

template <typename SCALAR>
const Eigen::Matrix<int, 9, 2>
    FeLagrangeO2Quad<SCALAR>::ksegment_to_quad_mapping_ =
        // clang-format off
  (Eigen::Matrix<int,9,2>() <<  0, 0,
                                1, 0,
                                1, 1,
                                0, 1,
                                2, 0,
                                1, 2,
                                2, 1,
                                0, 2,
                                2, 2).finished();
// clang-format on

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Cubic Lagrangian finite elment on a triangular reference element
 *
 * This is a specialization of ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The 10 reference shape functions are combinations of the barycentric
 * coordinate functions on the reference triangle.
 *
 * The first three shape functions are associated with vertices,
 * the next six with edges of the triangle. The last shape function is an
 * interior shape function, see @lref{ex:inlagfe} for the location of
 * interpolation nodes.
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO3Tria final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO3Tria(const FeLagrangeO3Tria&) = default;
  FeLagrangeO3Tria(FeLagrangeO3Tria&&) noexcept = default;
  FeLagrangeO3Tria& operator=(const FeLagrangeO3Tria&) = default;
  FeLagrangeO3Tria& operator=(FeLagrangeO3Tria&&) noexcept = default;
  FeLagrangeO3Tria() = default;
  ~FeLagrangeO3Tria() override = default;

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kTria();
  }

  /**  @brief Quadratic Lagrangian finite elements rely on polynomials of degree
   * 3
   */
  [[nodiscard]] unsigned Degree() const override { return 3; }

  /** @brief Ten local shape functions are associated with a triangular cell in
   *         the case of cubic Lagrangian finite elements.
   *
   * @sa ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 10; }

  /**
   * @brief One shape function attached to each node
   * and two to each edge of the triangle. There is one interior
   * shape function.
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 1) ? 2 : 1;
  }

  /**
   * @brief One shape function attached to each node
   * and two to each edge of the triangle. There is one interior
   * shape function.
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 1) ? 2 : 1;
  }

  /** @brief Point evaluation of reference shape functions
   *
   * The reference shape function are cardinal basis functions of the space of
   * cubic polynomials with respect to the ten interpolation nodes supplied by
   * the lattive drawn in @lref{lfedof3}.
   *
   * @sa ScalarReferenceFiniteElement::EvalReferenceShapeFunctions()
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    // number of points for which evaluation is desired
    const size_type n_pts(refcoords.cols());
    // Matrix for returning the values of the reference shape functions at the
    // specified nodes. Each row corresponds to a reference shape function
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(10, n_pts);

    // evaluation of the barycentric coordinate functions
    Eigen::Array<double, 1, Eigen::Dynamic> lambda0 =
        1 - refcoords.row(0).array() - refcoords.row(1).array();
    Eigen::Array<double, 1, Eigen::Dynamic> lambda1 = refcoords.row(0).array();
    Eigen::Array<double, 1, Eigen::Dynamic> lambda2 = refcoords.row(1).array();

    // evaluation of the shape functions
    // The LSF associated with vertices
    result.row(0) = 4.5 * lambda0 * (lambda0 - 1 / 3.0) * (lambda0 - 2 / 3.0);
    result.row(1) = 4.5 * lambda1 * (lambda1 - 1 / 3.0) * (lambda1 - 2 / 3.0);
    result.row(2) = 4.5 * lambda2 * (lambda2 - 1 / 3.0) * (lambda2 - 2 / 3.0);
    // Six LSF are attached to edges
    result.row(3) = 13.5 * lambda0 * lambda1 * (lambda0 - 1 / 3.0);
    result.row(4) = 13.5 * lambda0 * lambda1 * (lambda1 - 1 / 3.0);
    result.row(5) = 13.5 * lambda1 * lambda2 * (lambda1 - 1 / 3.0);
    result.row(6) = 13.5 * lambda1 * lambda2 * (lambda2 - 1 / 3.0);
    result.row(7) = 13.5 * lambda2 * lambda0 * (lambda2 - 1 / 3.0);
    result.row(8) = 13.5 * lambda2 * lambda0 * (lambda0 - 1 / 3.0);
    // LSF associated with the cell: cubic bubble function
    result.row(9) = 27.0 * lambda0 * lambda1 * lambda2;

    return result;
  }

  /** @brief Point evaluations of gradients of reference shape functions
   *
   *  Compute values of gradients of reference shape functions at specified
   * nodes
   *
   * @sa ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions()
   */
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

    const size_type n_pts(refcoords.cols());
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(10, 2 * n_pts);

    // reshape into a 20xn_pts matrix.
    Eigen::Map<Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
               Eigen::AutoAlign>
        temp(&result(0, 0), 20, n_pts);

    // evaulate barycentric coordinate functions:
    Eigen::Array<double, 1, Eigen::Dynamic> l0 =
        1 - refcoords.row(0).array() - refcoords.row(1).array();
    auto l1 = refcoords.row(0).array();
    auto l2 = refcoords.row(1).array();

    // vertex-associated LSF: numbers 0,1,2
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
    // edge-associated basis functions
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
    // cell-associated LSF
    temp.row(9) = 27.0 * (-l1 * l2 + l0 * l2);
    temp.row(19) = 27.0 * (-l1 * l2 + l0 * l1);

    return result;
  }

  /** @brief Evaluation nodes are the three vertices of the triangle, two
   * equispaced interio nodes for each edge, a the barycentre
   *
   * The reference shape functions as implemented by this class are a _cardinal
   * basis_ of the space of cubic two-variate polynomials  with respect to
   * these evaluation nodes.
   *
   * The numbering of evluation nodes is fixed by the numbering of the local
   * shape functions.
   *
   * @sa ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::Matrix<double, 2, 10> nodes;
    Eigen::Matrix<double, 2, 3> vertices;
    Eigen::Matrix<double, 2, 6> edges;
    Eigen::Matrix<double, 2, 1> center;

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
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Cubic Lagrangian finite element on a line segment
 *
 * This is a specialization of ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The first two shape functins are associated with vertices of the segment.
 * The last two are interior shape functions.
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO3Segment final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO3Segment(const FeLagrangeO3Segment&) = default;
  FeLagrangeO3Segment(FeLagrangeO3Segment&&) noexcept = default;
  FeLagrangeO3Segment& operator=(const FeLagrangeO3Segment&) = default;
  FeLagrangeO3Segment& operator=(FeLagrangeO3Segment&&) noexcept = default;
  FeLagrangeO3Segment() = default;
  ~FeLagrangeO3Segment() override = default;

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kSegment();
  }

  [[nodiscard]] unsigned Degree() const override { return 3; }

  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 4; }

  /** @brief One shape function attached to each node and two interior
   * shape functions.
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return (codim == 0) ? 2 : 1;
  }

  /** @brief One shape function attached to each node and two interior
   * shape functions.
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,
   * sub_idx_t)
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

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(4, n_pts);

    auto x = refcoords.row(0).array();

    // endpoints
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

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(4, n_pts);

    auto x = refcoords.row(0).array();

    // endpoints
    result.row(0) = -13.5 * x * x + 18.0 * x - 5.5;
    result.row(1) = 13.5 * x * x - 9.0 * x + 1.0;
    // interior
    result.row(2) = 40.5 * x * x - 45 * x + 9;
    result.row(3) = -40.5 * x * x + 36 * x - 4.5;

    return result;
  }

  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::Matrix<double, 1, 4> nodes;
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
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Cubic Lagrangian finite element on a quadrilateral reference
 * element.
 *
 * This is a specialization of ScalarReferenceFiniteElement.
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
 * nodes \f$ x_j \f$ and \f$ y_l \f$, see @lref{tplfedof}.
 *
 * The first four shape functions are associated with the vertices,
 * the next eight = 4*2 with edges of the quadrilateral. The last four shape
 * functions are interior shape functions.
 *
 * @see ScalarReferenceFiniteElement
 */
template <typename SCALAR>
class FeLagrangeO3Quad final
    : public lf::fe::ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO3Quad(const FeLagrangeO3Quad&) = default;
  FeLagrangeO3Quad(FeLagrangeO3Quad&&) noexcept = default;
  FeLagrangeO3Quad& operator=(const FeLagrangeO3Quad&) = default;
  FeLagrangeO3Quad& operator=(FeLagrangeO3Quad&&) noexcept = default;
  FeLagrangeO3Quad() = default;
  ~FeLagrangeO3Quad() override = default;

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kQuad();
  }

  [[nodiscard]] unsigned Degree() const override { return 3; }

  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 16; }

  /** @brief One shape function attached to each node and two to each edge
   * of the quadrilateral. Four interior shape functions.
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
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

  /** @brief One shape function attached to each node and two to each edge
   * of the quadrilateral. Four interior shape functions.
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t)
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

    // evaluate "segment" shape functions (b^j and b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x0_eval =
        (krsf_segment_.EvalReferenceShapeFunctions(refcoords.row(0))).array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x1_eval =
        (krsf_segment_.EvalReferenceShapeFunctions(refcoords.row(1))).array();

    // Evaluate basis functions using the tensor product structure
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
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(16, 2 * n_pts);

    // reshape into a 32xn_pts matrix.
    Eigen::Map<Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic>,
               Eigen::AutoAlign>
        temp(&result(0, 0), 32, n_pts);

    // evaluate "segment" shape functions (b^j and b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x0_eval =
        (krsf_segment_.EvalReferenceShapeFunctions(refcoords.row(0))).array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x1_eval =
        (krsf_segment_.EvalReferenceShapeFunctions(refcoords.row(1))).array();

    // evaluate derivatives of "segment" shape functions (b^j and
    // b^l)
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x0_grad =
        (krsf_segment_.GradientsReferenceShapeFunctions(refcoords.row(0)))
            .array();
    Eigen::Array<SCALAR, Eigen::Dynamic, Eigen::Dynamic> segment_x1_grad =
        (krsf_segment_.GradientsReferenceShapeFunctions(refcoords.row(1)))
            .array();

    // evaluate gradients  using the product rule and
    // the tensor product structure of the basis functions.
    for (int i = 0; i < 16; ++i) {
      // d/dx
      temp.row(i) = (segment_x0_grad.row(ksegment_to_quad_mapping_(i, 0)) *
                     segment_x1_eval.row(ksegment_to_quad_mapping_(i, 1)))
                        .matrix();
      // d/dy
      temp.row(i + 16) = (segment_x1_grad.row(ksegment_to_quad_mapping_(i, 1)) *
                          segment_x0_eval.row(ksegment_to_quad_mapping_(i, 0)))
                             .matrix();
    }

    return result;
  }

  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::Matrix<double, 2, 16> nodes;

    Eigen::Matrix<double, 2, 4> vertices;
    Eigen::Matrix<double, 2, 8> midpoints;
    Eigen::Matrix<double, 2, 4> interior;

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
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }

 private:
  /** description of the shape functions on the reference line segment
   * \f$ [0,1] \f$  that are used in the tensor product construction of the
   * shape functions*/
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
  const static Eigen::Matrix<int, 16, 2> ksegment_to_quad_mapping_;
};

template <typename SCALAR>
const FeLagrangeO3Segment<SCALAR> FeLagrangeO3Quad<SCALAR>::krsf_segment_ =
    FeLagrangeO3Segment<SCALAR>();

template <typename SCALAR>
const Eigen::Matrix<int, 16, 2>
    FeLagrangeO3Quad<SCALAR>::ksegment_to_quad_mapping_ =
        // clang-format off
  (Eigen::Matrix<int,16,2>() << 0, 0,
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

}  // namespace lf::uscalfe

#endif
