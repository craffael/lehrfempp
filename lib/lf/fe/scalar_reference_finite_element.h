#ifndef LF_FE_SCALAR_REFERENCE_FINITE_ELEMENT_H_
#define LF_FE_SCALAR_REFERENCE_FINITE_ELEMENT_H_

/**
 * @file
 * @brief Data structure representing simple finite elements
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>

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
 * @headerfile lf/fe/fe.h
 * @brief Interface class for parametric scalar valued finite elements
 *
 * @tparam SCALAR underlying scalar type, usually either `double` or
 * `complex<double>`
 *
 * A scalar parametric finite element is defined through a set of
 * reference shape functions (RSFs) on a particular reference entity, cf.
 * @lref{def:parfe}.
 *
 * Each reference shape function is associated with a unique sub-entity
 * of the reference entity according to @lref{eq:lbaff}.
 *
 * Specializations of this class support the evaluation of RSFs in
 * arbitrary points in the reference element and the computation of their
 * gradients. The also provide local components for the defintion of nodal
 * interpolants.
 *
 * This class is discussed in detail in @lref{par:lfppparfe}.
 *
 * ### Numbering _convention_ for reference shape functions
 *
 * The numbering of reference shape functions is done according to the
 * following convention, see also @lref{par:betlordlsf}
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
 *
 * The rules governing this numbering in LehrFEM++ are explained above and in
 * @lref{par:betlordlsf}.
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
   * @brief Request the maximal polynomial degree of the basis functions in this
   * finite element.
   *
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
   * @param codim co-dimension of the subentity
   * @return number of _interior_ reference shape function belonging to the
   *         specified sub-entity.
   *
   * @note this method will throw an exception when different numbers of
   *       reference shape functions are assigned to different sub-entities
   *       of the same co-dimension
   */
  // NOLINTNEXTLINE(misc-unused-parameters)
  [[nodiscard]] virtual size_type NumRefShapeFunctions(dim_t codim) const {
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
};

/**
 * @brief Print information about a ScalarReferenceFiniteElement to the given
 * stream.
 *
 * @param o stream to which output is to be sent
 * @param srfe The ScalarReferenceFiniteElement that should be printed.
 * @param ctrl controls the level of detail that is printed (see below)
 *
 *
 * #### Level of output:
 * - ctrl = 0: Only the type of the FiniteElementSpace, degree and number of
 * ShapeFunctions/Evaluation nodes are printed
 * - ctrl > 0: Also the coordinates of the evaluation nodes are
 * printed.
 *
 * @relates ScalarReferenceFiniteElement
 */
template <class SCALAR>
void PrintInfo(std::ostream& o,
               const ScalarReferenceFiniteElement<SCALAR>& srfe,
               unsigned int ctrl = 0) {
  o << typeid(srfe).name() << ", degree = " << srfe.Degree()
    << ", n_rsf = " << srfe.NumRefShapeFunctions()
    << ", n_evln = " << srfe.NumEvaluationNodes();
  if (ctrl > 0) {
    o << "\n evl nodes = " << srfe.EvaluationNodes();
  }
}

/** @brief Stream output operator: just calls the
 * ScalarReferenceFiniteElement::print() method
 * @relates ScalarReferenceFiniteElement
 */
template <typename SCALAR>
std::ostream& operator<<(std::ostream& o,
                         const ScalarReferenceFiniteElement<SCALAR>& fe_desc) {
  PrintInfo(o, fe_desc, 0);
  return o;
}

}  // end namespace lf::fe

/// \cond
/**
 * @brief Make lf::fe::ScalarReferenceFiniteElement formattable by fmt
 * (https://fmt.dev/latest/api.html#ostream-api)
 */
template <class SCALAR>
struct fmt::formatter<lf::fe::ScalarReferenceFiniteElement<SCALAR>>
    : ostream_formatter {};
/// \endcond

#endif  // LF_FE_SCALAR_REFERENCE_FINITE_ELEMENT_H_
