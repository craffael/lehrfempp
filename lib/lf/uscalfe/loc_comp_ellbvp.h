/**
 * @file
 * @brief Classes taking care of local computations for scalar 2nd-order
 * elliptic BVPs
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */
#include <map>

#ifndef LF_LOCCOMPELLBVP
#define LF_LOCCOMPELLBVP
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>

#include <iostream>

#include "precomputed_scalar_reference_finite_element.h"
#include "uscalfe.h"

namespace lf::uscalfe {
/** @brief Auxiliary data structure for passing collections of quadrature rules
 *
 * This type can be used to pass several quadrature rules to a function, when
 * different quadrature rules for different types of entities are required.
 */
using quad_rule_collection_t = std::map<lf::base::RefEl, lf::quad::QuadRule>;

/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Class for local quadrature based computations for Lagrangian finite
 * elements and second-order scalar elliptic BVPs.
 *
 * @tparam SCALAR scalar type of the UniformScalarFESpace. Must be a field
 *                     type such as `double` or `std::complex<double>`
 * @tparam DIFF_COEFF a \ref mesh_function "MeshFunction" that defines the
 *                    diffusion coefficient \f$ \mathbf{\alpha} \f$.
 *                    It should be either scalar- or matrix-valued.
 * @tparam REACTION_COEFF a \ref mesh_function "MeshFunction" that defines the
 *                        reaction coefficient \f$ \mathbf{\gamma} \f$.
 *                    It should be either scalar- or matrix-valued.
  *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * The element matrix is corresponds to the (local) bilinear form
 * @f[
    (u,v) \mapsto\int\limits_{K}\boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,u
          \cdot\mathbf{grad}\,v +
 \gamma(\mathbf{x})u\,\overline{v}\,\mathrm{d}\mathbf{x}
 \;,
 * @f]
 * with _diffusion coefficient_ @f$\mathbf{\alpha}@f$ and reaction coefficient
 * @f$\gamma@f$, see also @lref{ex:rdemp}
 *
 * ## Template parameter requirement
 *
 * - SCALAR must be a type like `double`
 * - DIFF_COEFF must provide an evaluation operator
 * `operator (const Entity &,ref_coord_t)` that returns either a scalar
 * or a matrix type that is compatible with Eigen's matrices. Usually it will
 * be an Eigen::Matrix either of variable of fixed size.
 *
 * @note The constructors of this class want an object of type @ref
 * UniformScalarFESpace, which holds a pointer to a mesh. However, for local
 * builder classes global information about the mesh is irrelevant, and,
 * therefore this object is used only to obtain information about the local
 * shape functions.
 * A revised implementation should directly pass this information to the
 * constructor.
 *
 *
 * ## Logger
 * This class logs additional information to
 * \ref ReactionDiffusionElementMatrixProviderLogger().
 * See \ref loggers for more information.
 */
template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
class ReactionDiffusionElementMatrixProvider {
  static_assert(mesh::utils::isMeshFunction<DIFF_COEFF>);
  static_assert(mesh::utils::isMeshFunction<REACTION_COEFF>);

 public:
  /**
   * @brief type of returned element matrix
   */
  using Scalar =
      typename decltype(mesh::utils::MeshFunctionReturnType<DIFF_COEFF>() *
                            Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>() +
                        mesh::utils::MeshFunctionReturnType<REACTION_COEFF>() *
                            Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>())::Scalar;
  using ElemMat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

  /** @brief standard constructors */
  /** @{ */
  ReactionDiffusionElementMatrixProvider(
      const ReactionDiffusionElementMatrixProvider &) = delete;
  ReactionDiffusionElementMatrixProvider(
      ReactionDiffusionElementMatrixProvider &&) noexcept = default;
  ReactionDiffusionElementMatrixProvider &operator=(
      const ReactionDiffusionElementMatrixProvider &) = delete;
  ReactionDiffusionElementMatrixProvider &operator=(
      ReactionDiffusionElementMatrixProvider &&) = delete;
  /** @} */

  /**
   * @brief Constructor: cell-independent precomputations
   *
   * @param fe_space collection of specifications for scalar-valued parametric
   * reference elements
   * @param alpha mesh function for the (possibly matrix-valued) diffusion
   * coefficient
   * @param gamma mesh function providing scalar-valued diffusion coefficient
   *
   * @see LocCompLagrFEPreprocessor::LocCompLagrFEPreprocessor()
   *
   * This constructor uses local quadature rules with double the degree of
   * exactness as the polynomial degree of the finite element space.
   */
  ReactionDiffusionElementMatrixProvider(
      std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space,
      DIFF_COEFF alpha, REACTION_COEFF gamma);

  /** @brief Constructor: cell-independent precomputations and custom quadrature
   * rule
   * @param fe_space collection of specifications for scalar-valued parametric
   * reference elements
   * @param alpha mesh function for the (possibly matrix-valued) diffusion
   * coefficient
   * @param gamma mesh function providing scalar-valued diffusion coefficient
   * @param qr_collection collection of quadrature rules. A quadrature rule is
   *  required for every cell type for which the finite element space provides
   *  local shape functions. If a quadrature rule is not specified for a cell
   *  type and the Eval() method is called for such a cell, exception will be
   * thrown.
   *
   * @see LocCompLagrFEPreprocessor::LocCompLagrFEPreprocessor()
   */
  ReactionDiffusionElementMatrixProvider(
      std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space,
      DIFF_COEFF alpha, REACTION_COEFF gamma,
      quad_rule_collection_t qr_collection);

  /**
   * @brief All cells are considered active in the default implementation
   *
   * This method is meant to be overloaded if assembly should be restricted to a
   * subset of cells.
   */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /**
   * @brief main routine for the computation of element matrices
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a small dense, containing the element matrix.
   *
   * Actual computation of the element matrix based on numerical quadrature and
   * mapping techniques. The order of the quadrature rule is tied to the
   * polynomial degree of the underlying Lagrangian finite element spaces: for
   * polynomial degree p a quadrature rule is chosen that is exact for
   * polynomials o degree 2p.
   *
   * @throw base::LfException in case the finite element specification is
   * missing for the type of the cell or if there is no quadrature rule
   * specified for the given cell type.
   */
  ElemMat Eval(const lf::mesh::Entity &cell);

  /** Virtual destructor */
  virtual ~ReactionDiffusionElementMatrixProvider() = default;

 private:
  /** @name functors providing coefficient functions
   * @{ */
  /** Diffusion coefficient */
  DIFF_COEFF alpha_;
  /** Reaction coefficient */
  REACTION_COEFF gamma_;
  /** @} */

  // fe_precomp_[i] contains precomputed reference finite element for ref_el i.
  std::array<PrecomputedScalarReferenceFiniteElement<SCALAR>, 5> fe_precomp_;
};

/**
 * @brief logger for ReactionDiffusionElementMatrixProvider
 */
std::shared_ptr<spdlog::logger> &ReactionDiffusionElementMatrixProviderLogger();

template <class PTR, class DIFF_COEFF, class REACTION_COEFF>
ReactionDiffusionElementMatrixProvider(PTR fe_space, DIFF_COEFF alpha,
                                       REACTION_COEFF gamma)
    -> ReactionDiffusionElementMatrixProvider<
        typename PTR::element_type::Scalar, DIFF_COEFF, REACTION_COEFF>;

template <class PTR, class DIFF_COEFF, class REACTION_COEFF>
ReactionDiffusionElementMatrixProvider(
    PTR fe_space, DIFF_COEFF alpha, REACTION_COEFF gamma,
    std::map<lf::base::RefEl, lf::quad::QuadRule>)
    -> ReactionDiffusionElementMatrixProvider<
        typename PTR::element_type::Scalar, DIFF_COEFF, REACTION_COEFF>;

// First constructor (internal construction of quadrature rules
template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
ReactionDiffusionElementMatrixProvider<SCALAR, DIFF_COEFF, REACTION_COEFF>::
    ReactionDiffusionElementMatrixProvider(
        std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space,
        DIFF_COEFF alpha, REACTION_COEFF gamma)
    : alpha_(std::move(alpha)), gamma_(std::move(gamma)), fe_precomp_() {
  for (auto ref_el : {base::RefEl::kTria(), base::RefEl::kQuad()}) {
    auto fe = fe_space->ShapeFunctionLayout(ref_el);
    // Check whether shape functions for that entity type are available.
    // Note that the corresponding PrecomputedScalarReferenceFiniteElement local
    // object is not initialized if the associated description of local shape
    // functions is missing.
    if (fe != nullptr) {
      // Precompute cell-independent quantities based on quadrature rules
      // with twice the degree of exactness compared to the degree of the
      // finite element space.
      fe_precomp_[ref_el.Id()] = PrecomputedScalarReferenceFiniteElement(
          fe, quad::make_QuadRule(ref_el, 2 * fe->Degree()));
    }
  }
}

// Second constructor (quadrature rules passed as arguments)
template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
ReactionDiffusionElementMatrixProvider<SCALAR, DIFF_COEFF, REACTION_COEFF>::
    ReactionDiffusionElementMatrixProvider(
        std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space,
        DIFF_COEFF alpha, REACTION_COEFF gamma,
        quad_rule_collection_t qr_collection)
    : alpha_(std::move(alpha)), gamma_(std::move(gamma)), fe_precomp_() {
  for (auto ref_el : {base::RefEl::kTria(), base::RefEl::kQuad()}) {
    // Obtain pointer to an object describing local shape functions
    auto fe = fe_space->ShapeFunctionLayout(ref_el);
    // Check whether shape functions for that entity type are available
    // Note that the corresponding PrecomputedScalarReferenceFiniteElement local
    // object is not initialized if the associated description of local shape
    // functions is missing.
    if (fe != nullptr) {
      // Obtain quadrature rule from user-supplied collection.
      auto qr_coll_ptr = qr_collection.find(ref_el);
      if (qr_coll_ptr != qr_collection.end()) {
        // A quadrature rule for the current entity type is available
        lf::quad::QuadRule qr = qr_coll_ptr->second;
        LF_ASSERT_MSG(qr.RefEl() == ref_el,
                      "qr.RefEl() = " << qr.RefEl() << " <-> " << ref_el);
        // Precomputations of cell-independent quantities
        fe_precomp_[ref_el.Id()] =
            PrecomputedScalarReferenceFiniteElement(fe, qr);
      }
    }
  }
}

// Main method for the computation of the element matrix
template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
typename lf::uscalfe::ReactionDiffusionElementMatrixProvider<
    SCALAR, DIFF_COEFF, REACTION_COEFF>::ElemMat
ReactionDiffusionElementMatrixProvider<
    SCALAR, DIFF_COEFF, REACTION_COEFF>::Eval(const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  // Obtain precomputed information about values of local shape functions
  // and their gradients at quadrature points.
  PrecomputedScalarReferenceFiniteElement<SCALAR> &pfe =
      fe_precomp_[ref_el.Id()];
  if (!pfe.isInitialized()) {
    // Accident: cell is of a type not covered by finite element
    // specifications or there is no quadrature rule available for this
    // reference element type
    std::stringstream temp;
    temp << "No local shape function information or no quadrature rule for "
            "reference element type "
         << ref_el;
    throw base::LfException(temp.str());
  }

  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  SPDLOG_LOGGER_TRACE(ReactionDiffusionElementMatrixProviderLogger(),
                      "{}, shape = \n{}", ref_el,
                      geo_ptr->Global(ref_el.NodeCoords()));

  // Physical dimension of the cell
  const dim_t world_dim = geo_ptr->DimGlobal();
  // Gram determinant at quadrature points
  const Eigen::VectorXd determinants(
      geo_ptr->IntegrationElement(pfe.Qr().Points()));
  LF_ASSERT_MSG(
      determinants.size() == pfe.Qr().NumPoints(),
      "Mismatch " << determinants.size() << " <-> " << pfe.Qr().NumPoints());
  // Fetch the transformation matrices for the gradients
  const Eigen::MatrixXd JinvT(
      geo_ptr->JacobianInverseGramian(pfe.Qr().Points()));
  LF_ASSERT_MSG(
      JinvT.cols() == 2 * pfe.Qr().NumPoints(),
      "Mismatch " << JinvT.cols() << " <-> " << 2 * pfe.Qr().NumPoints());
  LF_ASSERT_MSG(JinvT.rows() == world_dim,
                "Mismatch " << JinvT.rows() << " <-> " << world_dim);

  // compute values of coefficients alpha, gamma at quadrature points:
  auto alphaval = alpha_(cell, pfe.Qr().Points());
  auto gammaval = gamma_(cell, pfe.Qr().Points());

  // Element matrix
  ElemMat mat(pfe.NumRefShapeFunctions(), pfe.NumRefShapeFunctions());
  mat.setZero();

  // Loop over quadrature points
  for (base::size_type k = 0; k < pfe.Qr().NumPoints(); ++k) {
    const double w = pfe.Qr().Weights()[k] * determinants[k];
    // Transformed gradients
    const auto trf_grad(
        JinvT.block(0, 2 * static_cast<Eigen::Index>(k), world_dim, 2) *
        pfe.PrecompGradientsReferenceShapeFunctions()
            .block(0, 2 * k, mat.rows(), 2)
            .transpose());
    // Transformed gradients multiplied with coefficient
    const auto alpha_trf_grad(alphaval[k] * trf_grad);
    mat += w * (trf_grad.adjoint() * alpha_trf_grad +
                (gammaval[k] * pfe.PrecompReferenceShapeFunctions().col(k)) *
                    (pfe.PrecompReferenceShapeFunctions().col(k).adjoint()));
  }
  return mat;
}

/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Quadrature-based computation of local mass matrix for an edge
 *
 * @tparam SCALAR underlying scalar type, usually double or complex<double>
 * @tparam COEFF \ref mesh_function "MeshFunction" that defines the
 * scalar valued coefficient \f$ \gamma \f$
 * @tparam EDGESELECTOR predicate defining which edges are included
 *
 * This helper class corresponds to the the element matrix
 * for the bilinear form
 * @f[
 *     (u,v) \mapsto \int\limits_e
 * \gamma(x)u(x)\overline{v(x)}\,\mathrm{d}S(x)\;,
 * @f]
 * where @f$e@f$ is an edge of the mesh, and @f$\gamma@f$ a scalar-valued
 * coefficient function.
 *
 * #### Logger
 * This class logs additional information to
 * \ref MassEdgeMatrixProviderLogger().
 * See \ref loggers for more information.
 *
 */
template <typename SCALAR, typename COEFF, typename EDGESELECTOR>
class MassEdgeMatrixProvider {
 public:
  /// Scalar type of the element matrix
  using scalar_t =
      decltype(static_cast<SCALAR>(0) * static_cast<mesh::utils::MeshFunctionReturnType<COEFF>>(0));
  using ElemMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

  /** @name standard constructors
   * @{ */
  MassEdgeMatrixProvider(const MassEdgeMatrixProvider &) = delete;
  MassEdgeMatrixProvider(MassEdgeMatrixProvider &&) noexcept = default;
  MassEdgeMatrixProvider &operator=(const MassEdgeMatrixProvider &) = delete;
  MassEdgeMatrixProvider &operator=(MassEdgeMatrixProvider &&) = delete;
  /** @} */
  /**
   * @brief Constructor performing cell-independent initializations and choosing
   * a suitable 1D quadrature rule
   *
   * @param fe_space Describes the shapefunctions
   * @param gamma coefficient function through functor object
   * @param edge_selector predicate object selecting active to be covered in
   * the assembly
   *
   * This constructor chooses a local quadature rule with double the degree of
   * exactness as the polynomial degree of the finite element space.
   */
  MassEdgeMatrixProvider(
      std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space, COEFF gamma,
      EDGESELECTOR edge_selector = base::PredicateTrue{})
      : gamma_(std::move(gamma)),
        edge_sel_(std::move(edge_selector)),
        fe_precomp_() {
    auto fe = fe_space->ShapeFunctionLayout(base::RefEl::kSegment());
    LF_ASSERT_MSG(fe != nullptr, "No shape functions specified for edges");
    // Precompute entity-independent quantities based on a LehrFEM++ built-in
    // quadrature rule
    fe_precomp_ = PrecomputedScalarReferenceFiniteElement(
        fe, quad::make_QuadRule(base::RefEl::kSegment(), 2 * fe->Degree()));
  }
  /**
   * @brief Constructor performing cell-independent initializations
   *
   * @param fe_space Describes the shapefunctions
   * @param gamma coefficient function through functor object
   * @param quadrule quadrature rule for EDGE entities
   * @param edge_selector predicate object selecting active to be covered in
   * the assembly
   *
   * This constructor takes a user-supplied quadrature rule.
   */
  MassEdgeMatrixProvider(
      std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space, COEFF gamma,
      lf::quad::QuadRule quadrule,
      EDGESELECTOR edge_selector = base::PredicateTrue{})
      : gamma_(std::move(gamma)),
        edge_sel_(std::move(edge_selector)),
        fe_precomp_() {
    auto fe = fe_space->ShapeFunctionLayout(base::RefEl::kSegment());
    LF_ASSERT_MSG(fe != nullptr, "No shape functions specified for edges");
    LF_ASSERT_MSG(quadrule.RefEl() == base::RefEl::kSegment(),
                  "Quadrature rule not meant for EDGE entities!");
    // Precompute entity-independent quantities
    fe_precomp_ = PrecomputedScalarReferenceFiniteElement(fe, quadrule);
  }

  /**
   * @brief If true, then an edge is taken into account during assembly
   *
   * The information about "active" edges is supplied through the
   * `edge_selector` argument of the constructor.
   */
  bool isActive(const lf::mesh::Entity &edge) {
    LF_ASSERT_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
                  "Wrong type for an edge");
    return edge_sel_(edge);
  }

  /**
   * @brief actual computation of edge mass matrix
   *
   * @param edge reference to the edge for
   *        which the mass matrix is needed
   * @return a small dense matrix, containing the element matrix.
   *
   * Actual computation of the local edge mass based on numerical quadrature
   and
   * mapping techniques. The order of the quadrature rule is tied to the
   * polynomial degree of the underlying Lagrangian finite element spaces:
   for
   * polynomial degree p a quadrature rule is chosen that is exact for
   * polynomials o degree 2p.
   */
  ElemMat Eval(const lf::mesh::Entity &edge);

  virtual ~MassEdgeMatrixProvider() = default;

 private:
  COEFF gamma_;            // functor for coefficient
  EDGESELECTOR edge_sel_;  // Defines the active edges
  // Precomputed quantities at quadrature points
  PrecomputedScalarReferenceFiniteElement<SCALAR> fe_precomp_;
};

/**
 * @brief logger for MassEdgeMatrixProvider
 */
std::shared_ptr<spdlog::logger> &MassEdgeMatrixProviderLogger();

// deduction guide:
template <class PTR, class COEFF, class EDGESELECTOR = base::PredicateTrue>
MassEdgeMatrixProvider(PTR, COEFF coeff,
                       EDGESELECTOR edge_predicate = base::PredicateTrue{})
    -> MassEdgeMatrixProvider<typename PTR::element_type::Scalar, COEFF,
                              EDGESELECTOR>;

// Eval() method
// TODO(craffael) remove const once
// https://
// developercommunity.visualstudio.com/content/problem/180948/vs2017-155-c-cv-qualifiers-lost-on-type-alias-used.html
// is resolved
template <class SCALAR, class COEFF, class EDGESELECTOR>
typename MassEdgeMatrixProvider<SCALAR, COEFF, EDGESELECTOR>::ElemMat
MassEdgeMatrixProvider<SCALAR, COEFF, EDGESELECTOR>::Eval(
    const lf::mesh::Entity &edge) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{edge.RefEl()};
  LF_ASSERT_MSG(ref_el == lf::base::RefEl::kSegment(),
                "Edge must be of segment type");
  // Query the shape of the edge
  const lf::geometry::Geometry *geo_ptr = edge.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");

  // Obtain the metric factors for the quadrature points
  const Eigen::VectorXd determinants(
      geo_ptr->IntegrationElement(fe_precomp_.Qr().Points()));
  LF_ASSERT_MSG(determinants.size() == fe_precomp_.Qr().NumPoints(),
                "Mismatch " << determinants.size() << " <-> "
                            << fe_precomp_.Qr().NumPoints());

  // Element matrix
  ElemMat mat(fe_precomp_.NumRefShapeFunctions(),
              fe_precomp_.NumRefShapeFunctions());
  mat.setZero();

  auto gammaval = gamma_(edge, fe_precomp_.Qr().Points());

  // Loop over quadrature points
  for (long k = 0; k < determinants.size(); ++k) {
    // Build local matrix by summing rank-1 contributions
    // from quadrature points.
    const auto w =
        (fe_precomp_.Qr().Weights()[k] * determinants[k]) * gammaval[k];
    mat += ((fe_precomp_.PrecompReferenceShapeFunctions().col(k)) *
            (fe_precomp_.PrecompReferenceShapeFunctions().col(k).adjoint())) *
           w;
  }
  return mat;
}

/**
 * @ingroup entity_vector_provider
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Local computation of general element (load) vector for scalar
 finite
 * elements; volume contributions only
 *
 * @tparam SCALAR Scalar type of the Finite Element Space.
 * @tparam MESH_FUNCTION \ref mesh_function "MeshFunction" which defines the
 source
 * function \f$ f \f$
 *
 * The underlying local linear form is
 * @f[
      v \mapsto \int_K
 f(\mathbf{x})\,\overline{v(\mathbf{x})}\,\mathrm{d}\mathbf{x}\;,
 * @f]
 * where \f$f\f$ is supposed to be a locally continuous source function.
 *
 * Computation is based on a quadrature rules supplied by the LehrFEM++
 * lf::quad::QuadRule module.
 *
 * This class complies with the requirements for the template parameter
 * `ELEM_VEC_COMP` of the function AssembleVectorLocally().
 *
 * #### Logger
 * This class logs additional information to
 * \ref ScalarLoadElementVectorProvider().
 * See \ref loggers for more information.
 *
 */
template <typename SCALAR, typename MESH_FUNCTION>
class ScalarLoadElementVectorProvider {
  static_assert(mesh::utils::isMeshFunction<MESH_FUNCTION>);

 public:
  /// Scalar type of the element matrix
  using scalar_t =
      decltype(static_cast<SCALAR>(0) *
               static_cast<mesh::utils::MeshFunctionReturnType<MESH_FUNCTION>>(0));
  using ElemVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;

  /** @name standard constructors
   *@{*/
  ScalarLoadElementVectorProvider(const ScalarLoadElementVectorProvider &) =
      delete;
  ScalarLoadElementVectorProvider(ScalarLoadElementVectorProvider &&) noexcept =
      default;
  ScalarLoadElementVectorProvider &operator=(
      const ScalarLoadElementVectorProvider &) = delete;
  ScalarLoadElementVectorProvider &operator=(
      ScalarLoadElementVectorProvider &&) = delete;
  /**@}*/

  /** @brief Constructor, performs precomputations
   *
   * @param fe_space specification of local shape functions
   * @param f functor object for source function
   *
   * Uses quadrature rule of double the degree of exactness compared to the
   * degree of the finite element space.
   */
  ScalarLoadElementVectorProvider(
      std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space,
      MESH_FUNCTION f);
  /** @brief Constructor, performs precomputations based on user-supplied
   * quadrature rules.
   *
   * @param fe_space specification of local shape functions
   * @param f functor object for source function
   * @param qr_collection collection of quadrature rule.
   *
   */
  ScalarLoadElementVectorProvider(
      std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space,
      MESH_FUNCTION f, quad_rule_collection_t qr_collection);
  /** @brief Default implement: all cells are active */
  virtual bool isActive(const lf::mesh::Entity & /*cell*/) { return true; }
  /*
   * @brief Main method for computing the element vector
   *
   * @param cell current cell for which the element vector is desired
   * @return local load vector as column vector
   *
   */
  ElemVec Eval(const lf::mesh::Entity &cell);

  virtual ~ScalarLoadElementVectorProvider() = default;

 private:
  /** @brief An object providing the source function */
  MESH_FUNCTION f_;

  std::array<PrecomputedScalarReferenceFiniteElement<SCALAR>, 5> fe_precomp_;
};

/**
 * @brief logger used by ScalarLoadElementVectorProvider
 */
std::shared_ptr<spdlog::logger> &ScalarLoadElementVectorProviderLogger();

// Deduction guide
template <class PTR, class MESH_FUNCTION>
ScalarLoadElementVectorProvider(PTR fe_space, MESH_FUNCTION mf)
    -> ScalarLoadElementVectorProvider<typename PTR::element_type::Scalar,
                                       MESH_FUNCTION>;

// Constructors
template <typename SCALAR, typename MESH_FUNCTION>
ScalarLoadElementVectorProvider<SCALAR, MESH_FUNCTION>::
    ScalarLoadElementVectorProvider(
        std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space,
        MESH_FUNCTION f)
    : f_(std::move(f)) {
  for (auto ref_el : {base::RefEl::kTria(), base::RefEl::kQuad()}) {
    auto fe = fe_space->ShapeFunctionLayout(ref_el);
    // Check whether shape functions for that entity type are available
    if (fe != nullptr) {
      // Precompute cell-independent quantities based on quadrature rules
      // with twice the degree of exactness compared to the degree of the
      // finite element space.
      fe_precomp_[ref_el.Id()] =
          PrecomputedScalarReferenceFiniteElement<SCALAR>(
              fe, quad::make_QuadRule(ref_el, 2 * fe->Degree()));
    }
  }
}

template <typename SCALAR, typename MESH_FUNCTION>
ScalarLoadElementVectorProvider<SCALAR, MESH_FUNCTION>::
    ScalarLoadElementVectorProvider(
        std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space,
        MESH_FUNCTION f, quad_rule_collection_t qr_collection)
    : f_(std::move(f)) {
  for (auto ref_el : {base::RefEl::kTria(), base::RefEl::kQuad()}) {
    auto fe = fe_space->ShapeFunctionLayout(ref_el);
    // Check whether shape functions for that entity type are available
    if (fe != nullptr) {
      // Obtain quadrature rule from user-supplied collection.
      auto qr_coll_ptr = qr_collection.find(ref_el);
      if (qr_coll_ptr != qr_collection.end()) {
        // A quadrature rule for the current entity type is available
        lf::quad::QuadRule qr = qr_coll_ptr->second;
        LF_ASSERT_MSG(qr.RefEl() == ref_el,
                      "qr.RefEl() = " << qr.RefEl() << " <-> " << ref_el);
        // Precompute cell-independent quantities using the user-supplied
        // quadrature rules
        fe_precomp_[ref_el.Id()] =
            PrecomputedScalarReferenceFiniteElement<SCALAR>(fe, qr);
      } else {
        // Quadrature rule is missing for an entity type for which
        // local shape functions are available
        LF_ASSERT_MSG(false, "Quadrature rule missing for " << ref_el);
      }
    }
  }
}

// TODO(craffael) remove const once
// http://developercommunity.visualstudio.com/content/problem/180948/vs2017-155-c-cv-qualifiers-lost-on-type-alias-used.html
// is resolved
template <typename SCALAR, typename MESH_FUNCTION>
typename ScalarLoadElementVectorProvider<SCALAR, MESH_FUNCTION>::ElemVec
ScalarLoadElementVectorProvider<SCALAR, MESH_FUNCTION>::Eval(
    const lf::mesh::Entity &cell) {
  // Type for source function
  using source_fn_t = mesh::utils::MeshFunctionReturnType<MESH_FUNCTION>;
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  // Obtain precomputed information about values of local shape functions
  // and their gradients at quadrature points.
  auto &pfe = fe_precomp_[ref_el.Id()];
  // Accident: cell is of a type not coverence by finite element specifications
  LF_ASSERT_MSG(
      pfe.isInitialized(),
      "No local shape function information for entity type " << ref_el);

  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  SPDLOG_LOGGER_TRACE(ScalarLoadElementVectorProviderLogger(),
                      "{}, shape = \n{}", ref_el,
                      geo_ptr->Global(ref_el.NodeCoords()));

  // Obtain the metric factors for the quadrature points
  const Eigen::VectorXd determinants(
      geo_ptr->IntegrationElement(pfe.Qr().Points()));
  LF_ASSERT_MSG(
      determinants.size() == pfe.Qr().NumPoints(),
      "Mismatch " << determinants.size() << " <-> " << pfe.Qr().NumPoints());
  SPDLOG_LOGGER_TRACE(ScalarLoadElementVectorProviderLogger(),
                      "LOCVEC({}): Metric factors :\n{}", ref_el,
                      determinants.transpose());

  // Element vector
  ElemVec vec(pfe.NumRefShapeFunctions());
  vec.setZero();

  auto fval = f_(cell, pfe.Qr().Points());

  // Loop over quadrature points
  for (long k = 0; k < determinants.size(); ++k) {
    SPDLOG_LOGGER_TRACE(ScalarLoadElementVectorProviderLogger(),
                        "LOCVEC: [{}] -> [weight = {}]",
                        pfe.Qr().Points().transpose(), pfe.Qr().Weights()[k]);
    // Contribution of current quadrature point
    vec += (pfe.Qr().Weights()[k] * determinants[k] * fval[k]) *
           pfe.PrecompReferenceShapeFunctions().col(k).conjugate();
  }

  SPDLOG_LOGGER_TRACE(ScalarLoadElementVectorProviderLogger(), "LOCVEC = \n{}",
                      vec.transpose());

  return vec;
}

/**
 * @ingroup entity_vector_provider
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Local edge contributions to element vector
 *
 * @tparam SCALAR underlying scalar type, usually double or complex<double>
 * @tparam FUNCTOR `SCALAR` valued \ref mesh_function "MeshFunction" which
 * defines the function \f$ g \f$
 * @tparam EDGESELECTOR selector type for active edges
 *
 * The underlying local linear form for an edge @f$e@f$ is
 * @f[
    v \mapsto \int_e g(\mathbf{x})\,
 \overline{v(\mathbf{x})}\,\mathrm{d}S\mathbf{x}\;,
 * @f]
 * where \f$g\f$ is supposed to be a locally continuous source function.
 *
 * Computations are either based on a quadrature rules supplied by the LehrFEM++
 * lf::quad::QuadRule module or on a user-supplied quadrature rule.
 *
 * This class complies with the requirements for the template parameter
 * `ELEM_VEC_COMP` of the function AssembleVectorLocally().
 *
 * ### Logger
 * This class logs additional information to
 * \ref ScalarLoadEdgeVectorProvider().
 * See \ref loggers for more information.
 *
 * ### Type requirements
 *
 * - The EDGESELECTOR type must provide
 * ~~~
 bool operator(const lf::mesh::Entity &edge) const
 * ~~~
 * which returns true, if the edge is to be included in assembly.
 */
template <class SCALAR, class FUNCTOR, class EDGESELECTOR = base::PredicateTrue>
class ScalarLoadEdgeVectorProvider {
 public:
  static_assert(mesh::utils::isMeshFunction<FUNCTOR>,
                "FUNCTOR does not fulfill the concept of a mesh function.");
  using Scalar =
      decltype(static_cast<SCALAR>(0) * static_cast<mesh::utils::MeshFunctionReturnType<FUNCTOR>>(0));
  using ElemVec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

  /** @name standard constructors
   *@{*/
  ScalarLoadEdgeVectorProvider(const ScalarLoadEdgeVectorProvider &) = delete;
  ScalarLoadEdgeVectorProvider(ScalarLoadEdgeVectorProvider &&) noexcept =
      default;
  ScalarLoadEdgeVectorProvider &operator=(
      const ScalarLoadEdgeVectorProvider &) = delete;
  ScalarLoadEdgeVectorProvider &operator=(ScalarLoadEdgeVectorProvider &&) =
      delete;
  /**@}*/

  /** @brief Constructor, performs precomputations
   *
   * @param fe_space UniformScalarFESpace which describes the shape functions.
   * @param g functor object providing edge data
   * @param edge_sel selector predicate for active edges.
   *
   * This constructor selects one of LehrFEM++'s built-in quadrature rules
   * with a degree of exactness twice as big as the polynomial degree of the
   * finite element space.
   */
  ScalarLoadEdgeVectorProvider(
      std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space, FUNCTOR g,
      EDGESELECTOR edge_sel = base::PredicateTrue{})
      : g_(std::move(g)), edge_sel_(std::move(edge_sel)), pfe_() {
    auto fe = fe_space->ShapeFunctionLayout(base::RefEl::kSegment());
    LF_ASSERT_MSG(fe != nullptr, "No shape functions specified for edges");
    // Precompute entity-independent quantities based on a LehrFEM++ built-in
    // quadrature rule
    pfe_ = PrecomputedScalarReferenceFiniteElement(
        fe, quad::make_QuadRule(base::RefEl::kSegment(), 2 * fe->Degree()));
  }

  /** @brief Constructor, performs precomputations
   *
   * @param fe_space UniformScalarFESpace which describes the shape functions.
   * @param g functor object providing edge data
   * @param quadrule user-supplied quadrature rule object
   * @param edge_sel selector predicate for active edges.
   */
  ScalarLoadEdgeVectorProvider(
      std::shared_ptr<const UniformScalarFESpace<SCALAR>> fe_space, FUNCTOR g,
      lf::quad::QuadRule quadrule,
      EDGESELECTOR edge_sel = base::PredicateTrue{})
      : g_(std::move(g)), edge_sel_(std::move(edge_sel)), pfe_() {
    auto fe = fe_space->ShapeFunctionLayout(base::RefEl::kSegment());
    LF_ASSERT_MSG(fe != nullptr, "No shape functions specified for edges");
    LF_ASSERT_MSG(quadrule.RefEl() == base::RefEl::kSegment(),
                  "Quadrature rule not meant for EDGE entities!");
    // Precompute entity-independent quantities
    pfe_ = PrecomputedScalarReferenceFiniteElement(fe, quadrule);
  }

  /** @brief Default implement: all edges are active */
  virtual bool isActive(const lf::mesh::Entity &cell) {
    return edge_sel_(cell);
  }
  /*
   * @brief Main method for computing the element vector
   *
   * @param cell current cell for which the element vector is desired
   * @return local load vector as column vector
   *
   */
  ElemVec Eval(const lf::mesh::Entity &edge);

  virtual ~ScalarLoadEdgeVectorProvider() = default;

 private:
  FUNCTOR g_;              // source function
  EDGESELECTOR edge_sel_;  // selects edges
  PrecomputedScalarReferenceFiniteElement<SCALAR> pfe_;
};

/**
 * @brief logger for ScalarLoadEdgeVectorProvider class template.
 */
std::shared_ptr<spdlog::logger> &ScalarLoadEdgeVectorProviderLogger();

// deduction guide
template <class PTR, class FUNCTOR, class EDGESELECTOR = base::PredicateTrue>
ScalarLoadEdgeVectorProvider(PTR, FUNCTOR, EDGESELECTOR = base::PredicateTrue{})
    -> ScalarLoadEdgeVectorProvider<typename PTR::element_type::Scalar, FUNCTOR,
                                    EDGESELECTOR>;

// Eval() method
// TODO(craffael) remove const once
// https://developercommunity.visualstudio.com/content/problem/180948/vs2017-155-c-cv-qualifiers-lost-on-type-alias-used.html
// is resolved
template <class SCALAR, class FUNCTOR, class EDGESELECTOR>
typename ScalarLoadEdgeVectorProvider<SCALAR, FUNCTOR, EDGESELECTOR>::ElemVec
ScalarLoadEdgeVectorProvider<SCALAR, FUNCTOR, EDGESELECTOR>::Eval(
    const lf::mesh::Entity &edge) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{edge.RefEl()};
  LF_ASSERT_MSG(ref_el == lf::base::RefEl::kSegment(),
                "Edge must be of segment type");
  // Query the shape of the edge
  const lf::geometry::Geometry *geo_ptr = edge.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");

  // Quadrature points on physical edge
  const Eigen::MatrixXd mapped_qpts(geo_ptr->Global(pfe_.Qr().Points()));
  LF_ASSERT_MSG(
      mapped_qpts.cols() == pfe_.Qr().NumPoints(),
      "Mismatch " << mapped_qpts.cols() << " <-> " << pfe_.Qr().NumPoints());

  // Obtain the metric factors for the quadrature points
  const Eigen::VectorXd determinants(
      geo_ptr->IntegrationElement(pfe_.Qr().Points()));
  LF_ASSERT_MSG(
      determinants.size() == pfe_.Qr().NumPoints(),
      "Mismatch " << determinants.size() << " <-> " << pfe_.Qr().NumPoints());

  // Element vector
  ElemVec vec(pfe_.NumRefShapeFunctions());
  vec.setZero();

  auto g_vals = g_(edge, pfe_.Qr().Points());

  // Loop over quadrature points
  for (base::size_type k = 0; k < pfe_.Qr().NumPoints(); ++k) {
    // Add contribution of quadrature point to local vector
    const auto w = (pfe_.Qr().Weights()[k] * determinants[k]) * g_vals[k];
    vec += pfe_.PrecompReferenceShapeFunctions().col(k).conjugate() * w;
  }
  return vec;
}

}  // namespace lf::uscalfe

#endif
