#ifndef LF_FE
#define LF_FE
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Data structures representing simple finite elements
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include <lf/assemble/dofhandler.h>
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
 * @brief Interface class for parametric scalar valued finite elements
 *
 * @tparam SCALAR underlying scalar type, usually either `double` or
 * `complex<double>`
 *
 * A scalar parametric finite element is defined through a set of
 * reference shape functions (RSFs) on a particular reference entity.
 *
 * Each reference shape function is associated with a unique sub-entity
 * of the reference entity.
 *
 * Specializations of this class support the evaluation of RSFs in
 * arbitrary points in the reference element and the computation of their
 * gradients. The also provide local components for the defintion of nodal
 * interpolants.
 *
 * ### Numbering _convention_ for reference shape functions
 *
 * The numbering of reference shape functions is done according to the
 * following convention:
 *
 * -# RSFs  belonging to _sub-entities_ of a larger (relative)
 * co-dimension are arranged before dofs associated with _sub-entities_ of a
 * smaller co-dimension (first nodes, then edges, finally cells).
 * -# RSFs belonging to the same sub-entity are numbered contiguously.
 * -# RSFs for _sub-entities_ of the same co-dimension are taken into account in
 * the order given by the _local indexing_ of the sub-entities.
 *
 * The numbering scheme must the same as that adopted for the definition
 * of local-to-global maps, see the \ref lf::assemble::DofHandler class.
 *
 * ### Example for numbering of reference shape functions
 *
 * For triangular cubic Lagrangian finite elements there is a single reference
 * shape function associated with each vertex, two reference shape functions
 * belonging to every edge and a single interior reference shape function. In
 * this case the RSFs are numbered a follows
 * - RSF 0 -> vertex 0
 * - RSF 1 -> vertex 1
 * - RSF 2 -> vertex 2
 * - RSF 3,4 -> edge 0 (3 closer to endpoint 0, 4 closer to endpoint 1)
 * - RSF 5,6 -> edge 1 (5 closer to endpoint 0, 6 closer to endpoint 1)
 * - RSF 7,8 -> edge 2 (7 closer to endpoint 0, 8 closer to endpoint 1)
 * - RSF 9  -> triangle
 */
template <typename SCALAR>
class ScalarReferenceFiniteElement {
 protected:
  ScalarReferenceFiniteElement() = default;
  ScalarReferenceFiniteElement(const ScalarReferenceFiniteElement&) = default;
  // NOLINTNEXTLINE
  ScalarReferenceFiniteElement(ScalarReferenceFiniteElement&&) noexcept =
      default;
  // NOLINTNEXTLINE
  ScalarReferenceFiniteElement& operator=(const ScalarReferenceFiniteElement&) =
      default;
  // NOLINTNEXTLINE
  ScalarReferenceFiniteElement& operator=(
      ScalarReferenceFiniteElement&&) noexcept = default;

 public:
  /** @brief The underlying scalar type */
  using Scalar = SCALAR;

  /** Virtual destructor */
  virtual ~ScalarReferenceFiniteElement() = default;

  /**
   * @brief Tells the type of reference cell underlying the parametric finite
   * element
   */
  [[nodiscard]] virtual base::RefEl RefEl() const = 0;

  /**
   * @brief Request the polynomial degree of the finite element space
   *
   * @sa ScalarReferenceFiniteElement(lf::base::RefEl ref_el, unsigned int
   * degree)
   */
  [[nodiscard]] virtual unsigned int Degree() const = 0;

  /**
   * @brief Returns the spatial dimension of the reference cell
   */
  [[nodiscard]] dim_t Dimension() const { return RefEl().Dimension(); }

  /**
   * @brief Total number of reference shape functions associated with the
   * reference cell
   *
   * @note the _total_ number of shape functions is the sum of the number
   * of interior shape functions for all sub-entities and the entity itself.
   */
  [[nodiscard]] virtual size_type NumRefShapeFunctions() const {
    size_type cnt = 0;
    for (dim_t codim = 0; codim <= Dimension(); ++codim) {
      for (sub_idx_t subidx = 0; subidx < RefEl().NumSubEntities(codim);
           ++subidx) {
        cnt += NumRefShapeFunctions(codim, subidx);
      }
    }
    return cnt;
  }

  /**
   * @brief The number of _interior_ reference shape functions for sub-entities
   *        of a particular co-dimension
   *
   * @param codim do-dimension of the subentity
   * @return number of _interior_ reference shape function belonging to the
   *         specified sub-entity.
   *
   * @note this method will throw an exception when different numbers of
   *       reference shape functions are assigned to different sub-entities
   *       of the same co-dimension
   */
  [[nodiscard]] virtual size_type NumRefShapeFunctions(dim_t /*codim*/) const {
    LF_VERIFY_MSG(false, "Illegal call for non-uniform sub-entity dof numbers");
    return 0;
  }

  /**
   * @brief The number of _interior_ reference shape functions for every
   * sub-entity
   *
   * @param codim do-dimension of the subentity
   * @param subidx local index of the sub-entity
   * @return number of _interior_ reference shape function belonging to the
   *         specified sub-entity.
   */
  [[nodiscard]] virtual size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t subidx) const = 0;

  /**
   * @brief Evaluation of _all_ reference shape functions in a number of points
   *
   * @param refcoords coordinates of N points in the reference cell passed as
   * columns of a matrix of size dim x N, where dim is the dimension of the
   * reference element, that is =0 for points, =1 for edges, =2 for triangles
   * and quadrilaterals
   *
   * @return An Eigen Matrix of size `NumRefShapeFunctions() x refcoords.cols()`
   * which contains the shape functions evaluated at every quadrature point.
   *
   * Concerning the numbering of local shape functions, please consult
   * the documentation of lf::assemble::DofHandler or the documentation of the
   * class.
   *
   * @note shape functions are assumed to be real-valued.
   *
   * ### Example: Triangular Linear Lagrangian finite elements
   *
   * There are three reference shape functions
   * \f$\hat{b}^1,\hat{b}^2,\hat{b}^3\f$ associated with the vertices of the
   * reference triangle. Let us assume that the `refcoords` argument is a 2x2
   * matrix \f$[\mathbf{x}_1\;\mathbf{x}_2]\f$, which corresponds to passing
   * the coordinates of two points in the reference triangle. Then this method
   * will return a `3x2` matrix:
   * \f[
        \begin{pmatrix}\hat{b}^1(\mathbf{x}_1) & \hat{b}^1(\mathbf{x}_2) \\
          \hat{b}^2(\mathbf{x}_1) & \hat{b}^2(\mathbf{x}_2) \\
          \hat{b}^3(\mathbf{x}_1)\ & \hat{b}^3(\mathbf{x}_2) \end{pmatrix}
     \f]
   *
   */
  [[nodiscard]] virtual Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const = 0;

  /**
   * @brief Computation of the gradients of _all_ reference shape functions in a
   * number of points
   *
   * @param refcoords coordinates of N points in the reference cell passed as
   *                  columns of a matrix of size dim x N.
   * @return An Eigen Matrix of size `NumRefShapeFunctions() x (dim *
   * refcoords.cols())` where `dim` is the dimension of the reference
   * finite element. The gradients are returned in chunks of rows of this
   * matrix.
   *
   * Concerning the numbering of local shape functions, please consult
   * the documentation of lf::assemble::DofHandler.
   *
   * ### Example: Triangular Linear Lagrangian finite elements
   *
   * There are three reference shape functions
   * \f$\hat{b}^1,\hat{b}^2,\hat{b}^3\f$ associated with the vertices of the
   * reference triangle. Let us assume that the `refcoords` argument is a 2x2
   * matrix \f$[\mathbf{x}_1\;\mathbf{x}_2]\f$, which corresponds to passing
   * the coordinates of two points in the reference triangle. Then this method
   * will return a `3x4` matrix:
   * \f[
  \begin{pmatrix}
  \grad^T\hat{b}^1(\mathbf{x}_1) & \grad^T\hat{b}^1(\mathbf{x}_2) \\
  \grad^T\hat{b}^2(\mathbf{x}_1) & \grad^T\hat{b}^2(\mathbf{x}_2) \\
  \grad^T\hat{b}^3(\mathbf{x}_1) & \grad^T\hat{b}^3(\mathbf{x}_2)
  \end{pmatrix} \f]
   */
  [[nodiscard]] virtual Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  GradientsReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const = 0;

  /**
   * @brief Returns reference coordinates of "evaluation nodes" for evaluation
   *        of parametric degrees of freedom, nodal interpolation in the
   * simplest case.
   *
   * @return A d x N matrix, where d is the dimension of the reference cell,
   * and N the number of interpolation nodes. The columns of this matrix contain
   * their reference coordinates.
   *
   * Every parametric scalar finite element implicitly defines a local
   * interpolation operator by duality with the reference shape functions.
   * This interpolation operator can be realized through evaluations at certain
   * evaluation nodes, which are provided by this method.
   *
   * ### Unisolvence
   *
   * The evaluation points must satisfy the following requirement: If the values
   * of a function belonging to the span of the reference shape functions are
   * known in the evaluation nodes, then this function is uniquely determined.
   * This entails that the number of evaluation nodes must be at least as big as
   * the number of reference shape functions.
   *
   * @note It is not required that any vector a values at evaluation nodes can
   * be produced by a suitable linear combination of reference shape functions.
   *       For instance, this will not be possible, if there are more evaluation
   * points than reference shape functions. If both sets have the same size,
   * however, the interpolation problem has a solution for any vector of values
   * at the evluation points.
   *
   * ### Example: Principal lattice
   * For triangular Lagrangian finite elements of degree p the evaluation nodes,
   * usually called "interpolation nodes" in this context, can be chosen as
   * \f$\left(\frac{j}{p},\frac{k}{p}\right),\; 0\leq j,k \leq p, j+k\leq p\f$.
   *
   * ### Moment-based degrees of freedom
   * For some finite element spaces the interpolation functional may be defined
   * based on integrals over edges. In this case the evaluation nodes will be
   * quadrature nodes for the approximate evaluation of these integrals.
   *
   * The quadrature rule must be exact for the polynomials contained in the
   * local finite element spaces.
   */
  [[nodiscard]] virtual Eigen::MatrixXd EvaluationNodes() const = 0;

  /**
   * @brief Tell the number of evaluation (interpolation) nodes
   */
  [[nodiscard]] virtual size_type NumEvaluationNodes() const = 0;

  /**
   * @brief Computes the linear combination of reference shape functions
   *        matching function values at evaluation nodes.
   *
   * @param nodvals row vector of function values at evaluation nodes
   * The length of this vector must agree with NumEvaluationNodes().
   *
   * @return The coefficients of the local "nodal interpolant" with respect to
   * the reference shape functions. This is a row vector of length
   * NumRefShapeFunctions().
   *
   * If the evaluation nodes are interpolation nodes, that is, if the set of
   * reference shape functions forms a cardinal basis with respect to these
   * nodes, then we have NumEvaluationNodes() == NumRefShapeFunctions() and
   * the linear mapping realized by NodalValuesToDofs() is the identity mapping.
   *
   * @note default implementation is the identity mapping
   *
   * ### Requirement: reproduction of local finite element functions
   *
   * If the vector of values at the evaluation nodes agrees with a vector
   * of function values of a linear combination of reference shape functions
   * at the evaluation nodes, then this method must return the very coefficients
   * of the linear combination.
   */
  [[nodiscard]] virtual Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>
  NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>& nodvals) const {
    LF_ASSERT_MSG(nodvals.cols() == NumEvaluationNodes(),
                  "nodvals = " << nodvals << " <-> " << NumEvaluationNodes());
    return nodvals;
  }

  /**
   * @brief output function
   *
   * @param o stream to which output is to be sent
   * @return reference to the same stream
   *
   * Default implementation just prints degree, and the numbers of evaluation
   * nodes and local shape functions
   */
  virtual std::ostream& print(std::ostream& o) const {
    o << typeid(*this).name() << ", degree = " << Degree()
      << ", n_rsf = " << NumRefShapeFunctions()
      << ", n_evln = " << NumEvaluationNodes();
    if ((ctrl_ & kout_evln) != 0) {
      o << "\n evl nodes = " << EvaluationNodes();
    }
    return o;
  }

 public:
  /** @brief output control variable */
  static unsigned int ctrl_;
  static const unsigned int kout_evln = 1;
};

// Definition of output control variable
template <typename SCALAR>
unsigned int ScalarReferenceFiniteElement<SCALAR>::ctrl_ = 0;

/** @brief Stream output operator: just calls the
 * ScalarReferenceFiniteElement::print() method
 * @relates ScalarReferenceFiniteElement
 */
template <typename SCALAR>
std::ostream& operator<<(std::ostream& o,
                         const ScalarReferenceFiniteElement<SCALAR>& fe_desc) {
  return fe_desc.print(o);
}

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
 */
template <typename SCALAR>
class FeLagrangeO1Tria final : public ScalarReferenceFiniteElement<SCALAR> {
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

  /** @copydoc ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
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

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
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
 * The reference shape functions are
 * @f[
       \hat{b}^1(x_1,x_2) = (1-x_1)(1-x_2)\;,\quad
       \hat{b}^2(x_1,x_2) = x_1(1-x_2)\;,\quad
       \hat{b}^3(x_1,x_2) = x_1\cdot x_2\;,\quad
       \hat{b}^4(x_1,x_2) = (1-x_1)x_2\;.
 * @f]
 * This is a specialization of ScalarReferenceFiniteElement.
 * Refer to its documentation.
 */
template <typename SCALAR>
class FeLagrangeO1Quad final : public ScalarReferenceFiniteElement<SCALAR> {
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

  /** @copydoc ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
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

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
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
 */
template <typename SCALAR>
class FeLagrangeO1Segment final : public ScalarReferenceFiniteElement<SCALAR> {
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

  /** @copydoc ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
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

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
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
class FeLagrangePoint : public ScalarReferenceFiniteElement<SCALAR> {
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
 * This is a specialization of ScalarReferenceFiniteElement.
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
class FeLagrangeO2Tria final : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO2Tria(const FeLagrangeO2Tria&) = default;
  FeLagrangeO2Tria(FeLagrangeO2Tria&&) noexcept = default;
  FeLagrangeO2Tria& operator=(const FeLagrangeO2Tria&) = default;
  FeLagrangeO2Tria& operator=(FeLagrangeO2Tria&&) noexcept = default;
  FeLagrangeO2Tria() = default;
  ~FeLagrangeO2Tria() override = default;

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kTria();
  }

  [[nodiscard]] unsigned Degree() const override { return 2; }

  /**
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 6; }

  /**  @brief One shape function attached to each node
   * and one to each edge of the triangle. There are no interior
   * shape functions.
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 0 : 1;
  }

  /**  @brief One shape function attached to each node
   * and one to each edge of the triangle. There are no interior
   * shape functions.
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return (codim == 0) ? 0 : 1;
  }

  /** @copydoc ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
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

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
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
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    // clang-format off
    Eigen::Matrix<double, 2,6> nodes;
    nodes << 0.0, 1.0, 0.0, 0.5, 0.5, 0.0,
             0.0, 0.0, 1.0, 0.0, 0.5, 0.5;
    // clang-format on
    return nodes;
  }

  /** @brief Six evaluation nodes
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  [[nodiscard]] size_type NumEvaluationNodes() const override {
    return NumRefShapeFunctions();
  }
};

/**
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Quadratic Lagrangian finite element on a line segment
 *
 * This is a specialization of ScalarReferenceFiniteElement.
 * Refer to its documentation.
 *
 * The first two shape functions are associated with the vertices of the
 * segment. The last one is an interior shape function.
 */
template <typename SCALAR>
class FeLagrangeO2Segment : public ScalarReferenceFiniteElement<SCALAR> {
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

  [[nodiscard]] unsigned Degree() const override { return 2; }

  /** @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions() */
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 3; }

  /** @brief One shape function attached to each node and one interior shape
   * function.
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return 1;
  }

  /** @brief One shape function attached to each node and one interior shape
   * function.
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return 1;
  }

  /** @copydoc ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());

    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> result(3, n_pts);

    Eigen::ArrayXd x = refcoords.row(0).array();

    // endpoints
    result.row(0) = 2.0 * (1.0 - x) * (0.5 - x);
    result.row(1) = 2.0 * x * (x - 0.5);

    // midpoint
    result.row(2) = 4.0 * (1.0 - x) * x;

    return result;
  }

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
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
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::Matrix<double, 1, 3> nodes;
    nodes << 0.0, 1.0, 0.5;
    return nodes;
  }

  /** @brief Three evaluation nodes
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
 * Refer to its documentation
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
 */
template <typename SCALAR>
class FeLagrangeO2Quad final : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  FeLagrangeO2Quad(const FeLagrangeO2Quad&) = default;
  FeLagrangeO2Quad(FeLagrangeO2Quad&&) noexcept = default;
  FeLagrangeO2Quad& operator=(const FeLagrangeO2Quad&) = default;
  FeLagrangeO2Quad& operator=(FeLagrangeO2Quad&&) noexcept = default;
  FeLagrangeO2Quad() = default;
  ~FeLagrangeO2Quad() override = default;

  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kQuad();
  }

  [[nodiscard]] unsigned Degree() const override { return 2; }

  /** @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()*/
  [[nodiscard]] size_type NumRefShapeFunctions() const override { return 9; }

  /** @brief One shape function is attached to each node and edge of
   * the quadrilateral. One interior shape function.
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return 1;
  }

  /** @brief One shape function is attached to each node
   * and each edge of the quadrilateral. One interior shape function.
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,
   * sub_idx_t)
   */
  [[nodiscard]] size_type NumRefShapeFunctions(
      dim_t codim, sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim" << codim);
    return 1;
  }

  /** @copydoc ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  [[nodiscard]] Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");

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
                       segment_x0_eval.row(ksegment_to_quad_mapping_(i, 1)))
                          .matrix();
    }
    return result;
  }

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
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
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  [[nodiscard]] Eigen::MatrixXd EvaluationNodes() const override {
    Eigen::Matrix<double, 2, 9> nodes;

    Eigen::Matrix<double, 2, 4> vertices;
    Eigen::Matrix<double, 2, 4> midpoints;
    Eigen::Matrix<double, 2, 1> center;

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
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
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
  const static Eigen::MatrixXi ksegment_to_quad_mapping_;
};

template <typename SCALAR>
const FeLagrangeO2Segment<SCALAR> FeLagrangeO2Quad<SCALAR>::krsf_segment_ =
    FeLagrangeO2Segment<SCALAR>();

template <typename SCALAR>
const Eigen::MatrixXi FeLagrangeO2Quad<SCALAR>::ksegment_to_quad_mapping_ =
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

}  // namespace lf::uscalfe

#endif
