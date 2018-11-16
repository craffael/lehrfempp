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
 * of local-to-global maps, see lf::asssemble::DofHandler.
 *
 * Specializations of this class support the evaluation of RSFs in arbitrary
 * points in the reference element and the computation of their gradients.
 * The also provide local components for the defintion of nodal interpolants.
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
  /**
   * @brief Constructor setting topological type and order
   *
   * @param ref_el reference cell on which the finite element is defined
   * @param order "polynomial order" of the finite element; meant to provide
   *              a hint about the appropriate choice of quadrature rules
   *
   * The order of a scalar valued finite element will usually agree with the
   * degree (+1) of the largest full polynomial space contained in its local
   * space.
   */
  explicit ScalarReferenceFiniteElement(lf::base::RefEl ref_el,
                                        unsigned int order)
      : ref_el_(std::move(ref_el)), order_(order) {}

  virtual ~ScalarReferenceFiniteElement() = default;

  /**
   * @brief Tells the type of reference cell underlying the parametric finite
   * element
   */
  lf::base::RefEl RefEl() const { return ref_el_; }
  /**
   * @brief Request the "polynomial order" of the finite element space
   *
   * @sa ScalarReferenceFiniteElement(lf::base::RefEl ref_el, unsigned int
   * order)
   */
  unsigned int order() const { return order_; }

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
    LF_VERIFY_MSG(false,"Illegal call for non-uniform sub-entity dof numbers");
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
   * columns of a matrix of size dim x N.
   * @return array of row vectors. The i-th array element contains the values
   *         of the i-th shape function in the passed points.
   *
   * Concerning the numbering of local shape functions, please consult
   * the documentation of lf::assemble::DofHandler.
   */
  virtual std::vector<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const = 0;

  /**
   * @brief Computation of the gradients of _all_ reference shape functions in a
   * number of points
   *
   * @param refcoords coordinates of N points in the reference cell passed as
   *                  columns of a matrix of size dim x N.
   * @return array of matrices. The i-th array element contains the values
   *         of the gradients of the i-th shape function at the points passed
   *         in `refcoords` in its _columns_.
   *
   * Concerning the numbering of local shape functions, please consult
   * the documentation of lf::assemble::DofHandler.
   */
  virtual std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>>
  GradientsReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const = 0;

  /**
   * @brief Returns positions of "reference" points for nodal interpolation
   *
   * @return A d x N matrix, where d is the dimension of the reference cell,
   * and N the number of interpolation nodes. The columns of this matrix contain
   * their reference coordinates.
   *
   * Every parametric scalar finite element implicitly defines a local
   * interpolation operator by duality with the reference shape functions.
   * This interpolation operator can be realized through evaluations at certain
   * evaluation nodes, which are provided by this method.
   */
  virtual Eigen::MatrixXd EvaluationNodes() const = 0;

  /**
   * @brief Tell the number of evaluation nodes
   */
  virtual size_type NumEvaluationNodes() const = 0;

  /**
   * @brief Computes the local nodal interpolant from function values at
   * interpolation nodes.
   *
   * @param nodvals row vector of function values at interpolation nodes
   * The length of this vector must agree with NumEvaluationNodes().
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
   */
  virtual Eigen::Matrix<SCALAR, 1, Eigen::Dynamic> NodalValuesToDofs(
      const Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>& nodvals) const {
    LF_ASSERT_MSG(nodvals.cols() == NumEvaluationNodes(),
                  "nodvals = " << nodvals << " <-> " << NumEvaluationNodes());
    return nodvals;
  }

 protected:
  /** type of underlying reference cell */
  const lf::base::RefEl ref_el_;
  /** polynomial order */
  const unsigned int order_;
};

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
  std::vector<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const size_type n_pts(refcoords.cols());
    std::vector<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>> ret{
        (Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>::Constant(1, n_pts, 1.0) -
         refcoords.row(0) - refcoords.row(1)),
        refcoords.row(0), refcoords.row(1)};
    return ret;
  }

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const size_type n_pts(refcoords.cols());

    std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>> ret{
        (Eigen::Vector2d::Constant(-1.0).replicate(1, n_pts)),
        ((Eigen::Vector2d() << 1.0, 0.0).finished()).replicate(1, n_pts),
        ((Eigen::Vector2d() << 0.0, 1.0).finished()).replicate(1, n_pts)};
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

/**
 * @brief Linear Lagrange finite element on the quadrilateral reference element
 *
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
  std::vector<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const std::vector<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>> ret{
        ((1 - refcoords.row(0).array()) * (1 - refcoords.row(1).array()))
            .matrix(),
        ((refcoords.row(0).array()) * (1 - refcoords.row(1).array())).matrix(),
        ((refcoords.row(0).array()) * (refcoords.row(1).array())).matrix(),
        ((1 - refcoords.row(0).array()) * (refcoords.row(1).array())).matrix()};
    return ret;
  }

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 2,
                  "Reference coordinates must be 2-vectors");
    const size_type n_pts(refcoords.cols());

    std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>> ret{
        (Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>(2, n_pts)
             << (refcoords.row(1).array() - 1.0).matrix(),
         (refcoords.row(0).array() - 1.0).matrix())
            .finished(),
        (Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>(2, n_pts)
             << (1.0 - refcoords.row(1).array()).matrix(),
         (-refcoords.row(0).array()).matrix())
            .finished(),
        (Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>(2, n_pts)
             << (refcoords.row(1).array()).matrix(),
         (refcoords.row(0).array()).matrix())
            .finished(),
        (Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>(2, n_pts)
             << (-refcoords.row(1).array()).matrix(),
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
  std::vector<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>>
  EvalReferenceShapeFunctions(const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());
    std::vector<Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>> ret{
        (Eigen::Matrix<SCALAR, 1, Eigen::Dynamic>::Constant(1, n_pts, 1.0) -
         refcoords.row(0)),
        refcoords.row(0)};
    return ret;
  }

  /** @copydoc ScalarReferenceFiniteElement::GradientsReferenceShapeFunctions*/
  std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>>
  GradientsReferenceShapeFunctions(
      const Eigen::MatrixXd& refcoords) const override {
    LF_ASSERT_MSG(refcoords.rows() == 1,
                  "Reference coordinates must be 1-vectors");
    const size_type n_pts(refcoords.cols());

    std::vector<Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>> ret{
        (Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>::Constant(
            1, n_pts, -1.0)),
        (Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>::Constant(
            1, n_pts, 1.0))};
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
