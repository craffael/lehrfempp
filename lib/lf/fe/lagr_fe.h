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

namespace lf::fe {
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
 * @brief Interface class for parametric scalar valued finite elements
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
 * smaller co-dimension.
 * -# RSFss belonging to the same sub-entity are numbered contiguously.
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

  /** @brief The scalar type of the shape function */
  using Scalar = SCALAR;

  /**
   * @brief Constructor setting topological type and degree
   *
   * @param ref_el reference cell on which the finite element is defined
   * @param degree polynomial degree of the finite element; meant to provide
   *              a hint about the appropriate choice of quadrature rules
   *
   * The degree of a scalar valued finite element will usually agree with the
   * degree of the largest full polynomial space contained in its local
   * space.
   */
  explicit ScalarReferenceFiniteElement(lf::base::RefEl ref_el,
                                        unsigned int degree)
      : ref_el_(std::move(ref_el)), degree_(degree) {}

  virtual ~ScalarReferenceFiniteElement() = default;

  /**
   * @brief Tells the type of reference cell underlying the parametric finite
   * element
   */
  lf::base::RefEl RefEl() const { return ref_el_; }
  /**
   * @brief Request the polynomial degree of the finite element space
   *
   * @sa ScalarReferenceFiniteElement(lf::base::RefEl ref_el, unsigned int
   * degree)
   */
  unsigned int degree() const { return degree_; }

  /**
   * @brief Returns the spatial dimension of the reference cell
   */
  dim_t Dimension() const { return ref_el_.Dimension(); }

  /**
   * @brief Total number of reference shape functions associated with the
   * reference cell
   *
   * @note the _total_ number of shape functions is the sum of the number
   * of interior shape functions for all sub-entities and the entity itself.
   */
  virtual size_type NumRefShapeFunctions() const {
    size_type cnt = 0;
    for (dim_t codim = 0; codim <= Dimension(); ++codim) {
      for (sub_idx_t subidx = 0; subidx < ref_el_.NumSubEntities(codim);
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
  virtual size_type NumRefShapeFunctions(dim_t codim) const {
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
  virtual size_type NumRefShapeFunctions(dim_t codim,
                                         sub_idx_t subidx) const = 0;

  /**
   * @brief Evaluation of _all_ reference shape functions in a number of points
   *
   * @param refcoords coordinates of N points in the reference cell passed as
   * columns of a matrix of size dim x N, where dim is the dimension of the
   * reference element, that is =0 for points, =1 for edges, =2 for triangles
   * and quadrilaterals
   *
   * @return array of *row vectors*. The i-th array element contains the values
   *         of the i-th shape function in the passed points.
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
   * will return a vector of length 3 of row vectors of size 2:
   * \f[
        \left\{[\hat{b}^1(\mathbf{x}_1)\; \hat{b}^1(\mathbf{x}_2)],
          [\hat{b}^2(\mathbf{x}_1)\; \hat{b}^2(\mathbf{x}_2)],
          [\hat{b}^3(\mathbf{x}_1)\; \hat{b}^3(\mathbf{x}_2)]\right\}\;.
     \f]
   *
   */
  virtual std::vector<Eigen::RowVectorXd> EvalReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const = 0;

  /**
   * @brief Computation of the gradients of _all_ reference shape functions in a
   * number of points
   *
   * @param refcoords coordinates of N points in the reference cell passed as
   *                  columns of a matrix of size dim x N.
   * @return array of _matrices_. The i-th array element contains the values
   *         of the gradients of the i-th shape function at the points passed
   *         in `refcoords` in its _columns_.
   *
   * Concerning the numbering of local shape functions, please consult
   * the documentation of lf::assemble::DofHandler.
   */
  virtual std::vector<Eigen::MatrixXd> GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const = 0;

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
  virtual Eigen::MatrixXd EvaluationNodes() const = 0;

  /**
   * @brief Tell the number of evaluation (interpolation) nodes
   */
  virtual size_type NumEvaluationNodes() const = 0;

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
  virtual Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
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
    o << typeid(*this).name() << ", degree = " << degree()
      << ", n_rsf = " << NumRefShapeFunctions()
      << ", n_evln = " << NumEvaluationNodes();
    if ((ctrl_ & kout_evln) != 0) {
      o << "\n evl nodes = " << EvaluationNodes();
    }
    return o;
  }

 protected:
  /** type of underlying reference cell */
  const lf::base::RefEl ref_el_;
  /** polynomial degree */
  const unsigned int degree_;

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
 */
template <typename SCALAR>
std::ostream& operator<<(std::ostream& o,
                         const ScalarReferenceFiniteElement<SCALAR>& fe_desc) {
  return fe_desc.print(o);
}

/**
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
class TriaLinearLagrangeFE final : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  TriaLinearLagrangeFE(const TriaLinearLagrangeFE&) = default;
  TriaLinearLagrangeFE(TriaLinearLagrangeFE&&) noexcept = default;
  TriaLinearLagrangeFE& operator=(const TriaLinearLagrangeFE&) = default;
  TriaLinearLagrangeFE& operator=(TriaLinearLagrangeFE&&) noexcept = default;
  TriaLinearLagrangeFE()
      : ScalarReferenceFiniteElement<SCALAR>(lf::base::RefEl::kTria(), 1) {}
  virtual ~TriaLinearLagrangeFE() = default;

  /** @brief The local shape functions: barycentric coordinate functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  size_type NumRefShapeFunctions() const override { return 3; }

  /** @brief Exactly one shape function attached to each node, none to other
   * sub-entities
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim " << codim);
    return (codim == 2) ? 1 : 0;
  }
  /** @brief Exactly one shape function attached to each node, none to other
   * sub-entities
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t)
   */
  size_type NumRefShapeFunctions(dim_t codim,
                                 sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 2, "Illegal codim " << codim);
    return (codim == 2) ? 1 : 0;
  }

  /** @copydoc ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  std::vector<Eigen::RowVectorXd> EvalReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const size_type n_pts(refcoords.cols());
    std::vector<Eigen::RowVectorXd> ret{
        (Eigen::RowVectorXd::Constant(1, n_pts, 1.0) - refcoords.row(0) -
         refcoords.row(1)),
        refcoords.row(0), refcoords.row(1)};
    return ret;
  }

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  std::vector<Eigen::MatrixXd> GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const size_type n_pts(refcoords.cols());

    std::vector<Eigen::MatrixXd> ret{
        (Eigen::Vector2d::Constant(-1.0).replicate(1, n_pts)),
        ((Eigen::Vector2d() << 1.0, 0.0).finished()).replicate(1, n_pts),
        ((Eigen::Vector2d() << 0.0, 1.0).finished()).replicate(1, n_pts)};
    return ret;
  }

  /** @brief Evalutation nodes are just the vertices of the triangle
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  Eigen::MatrixXd EvaluationNodes() const override {
    return ScalarReferenceFiniteElement<double>::RefEl().NodeCoords();
  }

  /** @brief Three evaluation nodes
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  size_type NumEvaluationNodes() const override {
    return ScalarReferenceFiniteElement<double>::RefEl().NumNodes();
  }
};

/**
 * @brief Linear Lagrange finite element on the quadrilateral reference element
 *
 * The reference shape functions are
 * @f[
       \hat{b}^1(x_1,x_2) = (1-x_1)(1-x_2)\;,\quad
       \hat{b}^1(x_1,x_2) = x_1)(1-x_2)\;,\quad
       \hat{b}^1(x_1,x_2) = x_1\cdot x_2\;,\quad
       \hat{b}^1(x_1,x_2) = (1-x_1)x_2\;.
 * @f]
 * This is a specialization of ScalarReferenceFiniteElement.
 * Refer to its documentation.
 */
template <typename SCALAR>
class QuadLinearLagrangeFE final : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  QuadLinearLagrangeFE(const QuadLinearLagrangeFE&) = default;
  QuadLinearLagrangeFE(QuadLinearLagrangeFE&&) noexcept = default;
  QuadLinearLagrangeFE& operator=(const QuadLinearLagrangeFE&) = default;
  QuadLinearLagrangeFE& operator=(QuadLinearLagrangeFE&&) noexcept = default;
  QuadLinearLagrangeFE()
      : ScalarReferenceFiniteElement<SCALAR>(lf::base::RefEl::kQuad(), 1) {}
  virtual ~QuadLinearLagrangeFE() = default;

  /** @brief The local shape functions: four bilinear basis functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  size_type NumRefShapeFunctions() const override { return 4; }

  /** @brief Exactly one shape function attached to each node, none to other
   * sub-entities
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  size_type NumRefShapeFunctions(dim_t codim) const override {
    return (codim == 2) ? 1 : 0;
  }
  /** @brief Exactly one shape function attached to each node, none to other
   * sub-entities
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  size_type NumRefShapeFunctions(dim_t codim,
                                 sub_idx_t /*subidx*/) const override {
    return (codim == 2) ? 1 : 0;
  }

  /** @copydoc ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  std::vector<Eigen::RowVectorXd> EvalReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const std::vector<Eigen::RowVectorXd> ret{
        ((1 - refcoords.row(0).array()) * (1 - refcoords.row(1).array()))
            .matrix(),
        ((refcoords.row(0).array()) * (1 - refcoords.row(1).array())).matrix(),
        ((refcoords.row(0).array()) * (refcoords.row(1).array())).matrix(),
        ((1 - refcoords.row(0).array()) * (refcoords.row(1).array())).matrix()};
    return ret;
  }

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  std::vector<Eigen::MatrixXd> GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const size_type n_pts(refcoords.cols());

    std::vector<Eigen::MatrixXd> ret{
        (Eigen::MatrixXd(2, n_pts) << (refcoords.row(1).array() - 1.0).matrix(),
         (refcoords.row(0).array() - 1.0).matrix())
            .finished(),
        (Eigen::MatrixXd(2, n_pts) << (1.0 - refcoords.row(1).array()).matrix(),
         (-refcoords.row(0).array()).matrix())
            .finished(),
        (Eigen::MatrixXd(2, n_pts) << (refcoords.row(1).array()).matrix(),
         (refcoords.row(0).array()).matrix())
            .finished(),
        (Eigen::MatrixXd(2, n_pts) << (-refcoords.row(1).array()).matrix(),
         (1.0 - refcoords.row(0).array()).matrix())
            .finished()};
    return ret;
  }

  Eigen::MatrixXd EvaluationNodes() const override {
    return ScalarReferenceFiniteElement<SCALAR>::RefEl().NodeCoords();
  }

  size_type NumEvaluationNodes() const override {
    return ScalarReferenceFiniteElement<SCALAR>::RefEl().NumNodes();
  }
};

/**
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
class SegmentLinearLagrangeFE final
    : public ScalarReferenceFiniteElement<SCALAR> {
 public:
  SegmentLinearLagrangeFE(const SegmentLinearLagrangeFE&) = default;
  SegmentLinearLagrangeFE(SegmentLinearLagrangeFE&&) noexcept = default;
  SegmentLinearLagrangeFE& operator=(const SegmentLinearLagrangeFE&) = default;
  SegmentLinearLagrangeFE& operator=(SegmentLinearLagrangeFE&&) noexcept =
      default;
  SegmentLinearLagrangeFE()
      : ScalarReferenceFiniteElement<SCALAR>(lf::base::RefEl::kSegment(), 1) {}
  virtual ~SegmentLinearLagrangeFE() = default;

  /** @brief The local shape functions: barycentric coordinate functions
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions()
   */
  size_type NumRefShapeFunctions() const override { return 2; }

  /** @brief All shape functions attached to nodes
   * @copydoc ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t)
   */
  size_type NumRefShapeFunctions(dim_t codim) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return (codim == 1) ? 1 : 0;
  }
  /** @brief All shape functions attached to nodes
   * @copydoc
   * ScalarReferenceFiniteElement::NumRefShapeFunctions(dim_t,sub_idx_t)
   */
  size_type NumRefShapeFunctions(dim_t codim,
                                 sub_idx_t /*subidx*/) const override {
    LF_ASSERT_MSG(codim <= 1, "Illegal codim " << codim);
    return (codim == 1) ? 1 : 0;
  }

  /** @copydoc ScalarReferenceFiniteElement::EvalReferenceShapeFunctions() */
  std::vector<Eigen::RowVectorXd> EvalReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());
    std::vector<Eigen::RowVectorXd> ret{
        (Eigen::RowVectorXd::Constant(1, n_pts, 1.0) - refcoords.row(0)),
        refcoords.row(0)};
    return ret;
  }

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  std::vector<Eigen::MatrixXd> GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());

    std::vector<Eigen::MatrixXd> ret{
        (Eigen::MatrixXd::Constant(1, n_pts, -1.0)),
        (Eigen::MatrixXd::Constant(1, n_pts, 1.0))};
    return ret;
  }

  /** @brief Evalutation nodes are just the vertices of the triangle
   * @copydoc ScalarReferenceFiniteElement::EvaluationNodes()
   */
  Eigen::MatrixXd EvaluationNodes() const override {
    return ScalarReferenceFiniteElement<SCALAR>::RefEl().NodeCoords();
  }

  /** @brief Three evaluation nodes
   * @copydoc ScalarReferenceFiniteElement::NumEvaluationNodes()
   */
  size_type NumEvaluationNodes() const override {
    return ScalarReferenceFiniteElement<SCALAR>::RefEl().NumNodes();
  }
};

}  // namespace lf::fe

#endif
