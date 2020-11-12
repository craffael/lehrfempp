/**
 * @file
 * @brief Classes taking care of local computations for scalar 2nd-order
 * elliptic BVPs
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */
#include <map>

#ifndef LF_FE_LOCCOMPELLBVP
#define LF_FE_LOCCOMPELLBVP
/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <iostream>
#include "scalar_fe_space.h"
#include "scalar_reference_finite_element.h"

namespace lf::fe {
/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/uscalfe/uscalfe.h
 * @brief Class for local quadrature based computations for Lagrangian finite
 * elements and second-order scalar elliptic BVPs.
 *
 * @tparam SCALAR type for the entries of the element matrices. Must be a field
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
          \cdot\mathbf{grad}\,v + \gamma(\mathbf{x})u\,v\,\mathrm{d}\mathbf{x}
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
 */
template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
class ReactionDiffusionElementMatrixProvider {
  static_assert(mesh::utils::isMeshFunction<DIFF_COEFF>);
  static_assert(mesh::utils::isMeshFunction<REACTION_COEFF>);

 public:
  /**
   * @brief type of returned element matrix
   */
  using ElemMat = Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic>;

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
      std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space, DIFF_COEFF alpha,
      REACTION_COEFF gamma);

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
   * Throws an assertion in case the finite element specification is missing for
   * the type of the cell.
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

  /** FE Space */
  std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space_;
};

/**
 * @brief logger for ReactionDiffusionElementMatrixProvider
 */
extern std::shared_ptr<spdlog::logger>
    reaction_diffusion_element_matrix_provider_logger;

template <class PTR, class DIFF_COEFF, class REACTION_COEFF>
ReactionDiffusionElementMatrixProvider(PTR fe_space, DIFF_COEFF alpha,
                                       REACTION_COEFF gamma)
    ->ReactionDiffusionElementMatrixProvider<typename PTR::element_type::Scalar,
                                             DIFF_COEFF, REACTION_COEFF>;

// First constructor (internal construction of quadrature rules
template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
ReactionDiffusionElementMatrixProvider<SCALAR, DIFF_COEFF, REACTION_COEFF>::
    ReactionDiffusionElementMatrixProvider(
        std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space, DIFF_COEFF alpha,
        REACTION_COEFF gamma)
    : alpha_(std::move(alpha)),
      gamma_(std::move(gamma)),
      fe_space_(std::move(fe_space)) {}

// TODO(craffael) remove const once
// https://developercommunity.visualstudio.com/content/problem/180948/vs2017-155-c-cv-qualifiers-lost-on-type-alias-used.html
// is resolved
template <typename SCALAR, typename DIFF_COEFF, typename REACTION_COEFF>
typename lf::fe::ReactionDiffusionElementMatrixProvider<SCALAR, DIFF_COEFF,
                                                        REACTION_COEFF>::ElemMat
ReactionDiffusionElementMatrixProvider<
    SCALAR, DIFF_COEFF, REACTION_COEFF>::Eval(const lf::mesh::Entity &cell) {
  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  SPDLOG_LOGGER_TRACE(reaction_diffusion_element_matrix_provider_logger,
                      "{}, shape = \n{}", cell.RefEl(),
                      geo_ptr->Global(cell.RefEl().NodeCoords()));
  // Physical dimension of the cell
  const dim_t world_dim = geo_ptr->DimGlobal();

  // Get a quadrature rule of sufficiently high degree on the element
  const auto sfl = fe_space_->ShapeFunctionLayout(cell);
  const lf::quad::QuadRule qr =
      lf::quad::make_QuadRule(cell.RefEl(), 2 * sfl->Degree());

  const Eigen::VectorXd determinants(geo_ptr->IntegrationElement(qr.Points()));
  LF_ASSERT_MSG(
      determinants.size() == qr.NumPoints(),
      "Mismatch " << determinants.size() << " <-> " << qr.NumPoints());
  // Fetch the transformation matrices for the gradients
  const Eigen::MatrixXd JinvT(geo_ptr->JacobianInverseGramian(qr.Points()));
  LF_ASSERT_MSG(JinvT.cols() == 2 * qr.NumPoints(),
                "Mismatch " << JinvT.cols() << " <-> " << 2 * qr.NumPoints());
  LF_ASSERT_MSG(JinvT.rows() == world_dim,
                "Mismatch " << JinvT.rows() << " <-> " << world_dim);

  // compute values of alpha, gamma at quadrature points:
  auto alphaval = alpha_(cell, qr.Points());
  auto gammaval = gamma_(cell, qr.Points());

  // Element matrix
  ElemMat mat(sfl->NumRefShapeFunctions(), sfl->NumRefShapeFunctions());
  mat.setZero();

  // Compute the reference shape functions
  const auto rsf = sfl->EvalReferenceShapeFunctions(qr.Points());
  const auto grsf = sfl->GradientsReferenceShapeFunctions(qr.Points());

  // Loop over quadrature points
  for (base::size_type k = 0; k < qr.NumPoints(); ++k) {
    const double w = qr.Weights()[k] * determinants[k];
    // Transformed gradients
    const auto trf_grad(JinvT.block(0, 2 * k, world_dim, 2) *
                        grsf.block(0, 2 * k, mat.rows(), 2).transpose());
    // Transformed gradients multiplied with coefficient
    const auto alpha_trf_grad(alphaval[k] * trf_grad);
    mat += w * (alpha_trf_grad.transpose() * trf_grad +
                (gammaval[k] * rsf.col(k)) * (rsf.col(k).transpose()));
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
 *     (u,v) \mapsto \int\limits_e \gamma(x)u(x)v(x)\,\mathrm{d}S(x)\;,
 * @f]
 * where @f$e@f$ is an edge of the mesh, and @f$\gamma@f$ a scalar-valued
 * coefficient function.
 *
 */
template <typename SCALAR, typename COEFF, typename EDGESELECTOR>
class MassEdgeMatrixProvider {
 public:
  using scalar_t =
      decltype(SCALAR(0) * mesh::utils::MeshFunctionReturnType<COEFF>(0));
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
  MassEdgeMatrixProvider(std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space,
                         COEFF gamma,
                         EDGESELECTOR edge_selector = base::PredicateTrue{})
      : gamma_(std::move(gamma)),
        edge_sel_(std::move(edge_selector)),
        fe_space_(fe_space) {}

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
  // FE Space
  std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space_;
};

/**
 * @brief logger for MassEdgeMatrixProvider
 */
extern std::shared_ptr<spdlog::logger> mass_edge_matrix_provider_logger;

// deduction guide:
template <class PTR, class COEFF, class EDGESELECTOR = base::PredicateTrue>
MassEdgeMatrixProvider(PTR, COEFF coeff,
                       EDGESELECTOR edge_predicate = base::PredicateTrue{})
    ->MassEdgeMatrixProvider<typename PTR::element_type::Scalar, COEFF,
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
  const lf::geometry::Geometry *geo_ptr = edge.Geometry();
  // Get the shape function layout on the edge
  const auto sfl = fe_space_->ShapeFunctionLayout(edge);
  // Compute a quadrature rule on the given entity
  const lf::quad::QuadRule qr =
      lf::quad::make_QuadRule(edge.RefEl(), 2 * sfl->Degree());
  // Obtain the metric factors for the quadrature points
  const Eigen::VectorXd determinants(geo_ptr->IntegrationElement(qr.Points()));
  LF_ASSERT_MSG(
      determinants.size() == qr.NumPoints(),
      "Mismatch " << determinants.size() << " <-> " << qr.NumPoints());

  // Element matrix
  ElemMat mat(sfl->NumRefShapeFunctions(), sfl->NumRefShapeFunctions());
  mat.setZero();

  auto gammaval = gamma_(edge, qr.Points());

  // Compute the reference shape functions
  const auto rsf = sfl->EvalReferenceShapeFunctions(qr.Points());
  // Loop over quadrature points
  for (long k = 0; k < determinants.size(); ++k) {
    // Build local matrix by summing rank-1 contributions
    // from quadrature points.
    const auto w = (qr.Weights()[k] * determinants[k]) * gammaval[k];
    mat += ((rsf.col(k)) * (rsf.col(k).transpose())) * w;
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
 * @tparam SCALAR underlying scalar type, usually double or complex<double>
 * @tparam MESH_FUNCTION \ref mesh_function "MeshFunction" which defines the
 source
 * function \f$ f \f$
 *
 * The underlying local linear form is
 * @f[
      v \mapsto \int_K f(\mathbf{x})\,v(\mathbf{x}\,\mathrm{d}\mathbf{x}\;,
 * @f]
 * where \f$f\f$ is supposed to be a locally continuous source function.
 *
 * Computation is based on a quadrature rules supplied by the LehrFEM++
 * lf::quad::QuadRule module.
 *
 * This class complies with the requirements for the template parameter
 * `ELEM_VEC_COMP` of the function AssembleVectorLocally().
 */
template <typename SCALAR, typename MESH_FUNCTION>
class ScalarLoadElementVectorProvider {
  static_assert(mesh::utils::isMeshFunction<MESH_FUNCTION>);

 public:
  using ElemVec = Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>;

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
      std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space, MESH_FUNCTION f);
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

  std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space_;
};

/**
 * @brief logger used by ScalarLoadElementVectorProvider
 */
extern std::shared_ptr<spdlog::logger>
    scalar_load_element_vector_provider_logger;

// Deduction guide
template <class PTR, class MESH_FUNCTION>
ScalarLoadElementVectorProvider(PTR fe_space, MESH_FUNCTION mf)
    ->ScalarLoadElementVectorProvider<typename PTR::element_type::Scalar,
                                      MESH_FUNCTION>;

// Constructors
template <typename SCALAR, typename FUNCTOR>
ScalarLoadElementVectorProvider<SCALAR, FUNCTOR>::
    ScalarLoadElementVectorProvider(
        std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space, FUNCTOR f)
    : f_(std::move(f)), fe_space_(std::move(fe_space)) {}

// TODO(craffael) remove const once
// http://developercommunity.visualstudio.com/content/problem/180948/vs2017-155-c-cv-qualifiers-lost-on-type-alias-used.html
// is resolved
template <typename SCALAR, typename MESH_FUNCTION>
typename ScalarLoadElementVectorProvider<SCALAR, MESH_FUNCTION>::ElemVec
ScalarLoadElementVectorProvider<SCALAR, MESH_FUNCTION>::Eval(
    const lf::mesh::Entity &cell) {
  // Type for source function
  using source_fn_t = mesh::utils::MeshFunctionReturnType<MESH_FUNCTION>;

  // Get the shape function layout for the given cell
  const auto sfl = fe_space_->ShapeFunctionLayout(cell);

  // Initialize a quadrature rule of sufficiently high degree
  const lf::quad::QuadRule qr =
      lf::quad::make_QuadRule(cell.RefEl(), 2 * sfl->Degree());

  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  SPDLOG_LOGGER_TRACE(scalar_load_element_vector_provider_logger,
                      "{}, shape = \n{}", cell.RefEl(),
                      geo_ptr->Global(cell.RefEl().NodeCoords()));

  // Obtain the metric factors for the quadrature points
  const Eigen::VectorXd determinants(geo_ptr->IntegrationElement(qr.Points()));
  LF_ASSERT_MSG(
      determinants.size() == qr.NumPoints(),
      "Mismatch " << determinants.size() << " <-> " << qr.NumPoints());
  SPDLOG_LOGGER_TRACE(scalar_load_element_vector_provider_logger,
                      "LOCVEC({}): Metric factors :\n{}", cell.RefEl(),
                      determinants.transpose());
  // Element vector
  ElemVec vec(sfl->NumRefShapeFunctions());
  vec.setZero();

  auto fval = f_(cell, qr.Points());

  // Compute the reference shape functions
  const auto rsf = sfl->EvalReferenceShapeFunctions(qr.Points());

  // Loop over quadrature points
  for (long k = 0; k < determinants.size(); ++k) {
    SPDLOG_LOGGER_TRACE(scalar_load_element_vector_provider_logger,
                        "LOCVEC: [{}] -> [weight = {}]",
                        qr.Points().transpose(), qr.Weights()[k]);
    // Contribution of current quadrature point
    vec += (qr.Weights()[k] * determinants[k] * fval[k]) * rsf.col(k);
  }

  SPDLOG_LOGGER_TRACE(scalar_load_element_vector_provider_logger,
                      "LOCVEC = \n{}", vec.transpose());
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
    v \mapsto \int_e g(\mathbf{x})\,v(\mathbf{x})\,\mathrm{d}S\mathbf{x}\;,
 * @f]
 * where \f$g\f$ is supposed to be a locally continuous source function.
 *
 * Computations are either based on a quadrature rules supplied by the LehrFEM++
 * lf::quad::QuadRule module or on a user-supplied quadrature rule.
 *
 * This class complies with the requirements for the template parameter
 * `ELEM_VEC_COMP` of the function AssembleVectorLocally().
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
  using ElemVec = Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>;

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
   * @param fe_edge_p FE specification on edge
   * @param g functor object providing edge data
   * @param edge_sel selector predicate for active edges.
   *
   * This constructor selects one of LehrFEM++'s built-in quadrature rules
   * with a degree of exactness twice as big as the polynomial degree of the
   * finite element space.
   */
  ScalarLoadEdgeVectorProvider(
      std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space, FUNCTOR g,
      EDGESELECTOR edge_sel = base::PredicateTrue{})
      : g_(std::move(g)), edge_sel_(std::move(edge_sel)), fe_space_(fe_space) {}

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
  std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space_;
};

/**
 * @brief logger for ScalarLoadEdgeVectorProvider class template.
 */
extern std::shared_ptr<spdlog::logger> scalar_load_edge_vector_provider_logger;

// deduction guide
template <class PTR, class FUNCTOR, class EDGESELECTOR = base::PredicateTrue>
ScalarLoadEdgeVectorProvider(PTR, FUNCTOR, EDGESELECTOR = base::PredicateTrue{})
    ->ScalarLoadEdgeVectorProvider<typename PTR::element_type::Scalar, FUNCTOR,
                                   EDGESELECTOR>;

// Eval() method
// TODO(craffael) remove const once
// https://developercommunity.visualstudio.com/content/problem/180948/vs2017-155-c-cv-qualifiers-lost-on-type-alias-used.html
// is resolved
template <class SCALAR, class FUNCTOR, class EDGESELECTOR>
typename ScalarLoadEdgeVectorProvider<SCALAR, FUNCTOR, EDGESELECTOR>::ElemVec
ScalarLoadEdgeVectorProvider<SCALAR, FUNCTOR, EDGESELECTOR>::Eval(
    const lf::mesh::Entity &edge) {
  // Query the shape of the edge
  const lf::geometry::Geometry *geo_ptr = edge.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");

  // Get the shape function layout of the given edge
  const auto sfl = fe_space_->ShapeFunctionLayout(edge);

  // Quadrature points on physical edge
  const lf::quad::QuadRule qr =
      lf::quad::make_QuadRule(edge.RefEl(), 2 * sfl->Degree());
  const Eigen::MatrixXd mapped_qpts(geo_ptr->Global(qr.Points()));
  LF_ASSERT_MSG(mapped_qpts.cols() == qr.NumPoints(),
                "Mismatch " << mapped_qpts.cols() << " <-> " << qr.NumPoints());

  // Obtain the metric factors for the quadrature points
  const Eigen::VectorXd determinants(geo_ptr->IntegrationElement(qr.Points()));
  LF_ASSERT_MSG(
      determinants.size() == qr.NumPoints(),
      "Mismatch " << determinants.size() << " <-> " << qr.NumPoints());

  // Element vector
  ElemVec vec(sfl->NumRefShapeFunctions());
  vec.setZero();

  auto g_vals = g_(edge, qr.Points());

  // Compute the reference shape functions
  const auto rsf = sfl->EvalReferenceShapeFunctions(qr.Points());

  // Loop over quadrature points
  for (base::size_type k = 0; k < qr.NumPoints(); ++k) {
    // Add contribution of quadrature point to local vector
    const auto w = (qr.Weights()[k] * determinants[k]) * g_vals[k];
    vec += rsf.col(k) * w;
  }
  return vec;
}

}  // namespace lf::fe

#endif
