/**
 * @file
 * @brief Classes taking care of local computations for scalar 2nd-order
 * elliptic BVPs
 * @author Ralf Hiptmair
 * @date October 2018
 * @copyright MIT License
 */
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
#include <map>

#include "scalar_fe_space.h"
#include "scalar_reference_finite_element.h"

namespace lf::fe {
/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/fe/fe.h
 * @brief Class for computing element matrices for general scalar-valued finite
 * elements and homogeneous 2nd-order elliptic bilinear forms
 *
 * @tparam SCALAR scalar type of the ScalarFESpace Must be a field
 *                     type such as `double` or `std::complex<double>`
 * @tparam DIFF_COEFF a \ref mesh_function "MeshFunction" that defines the
 *                    diffusion coefficient \f$ \mathbf{\alpha} \f$.
 *                    It should be either scalar- or matrix-valued.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * The element matrix is corresponds to the (local) bilinear form
 * @f[
    (u,v) \mapsto\int\limits_{K}\boldsymbol{\alpha}(\mathbf{x})\mathbf{grad}\,u
          \cdot\mathbf{grad}\,v\,\mathrm{d}\mathbf{x}
 \;,
 * @f]
 * with _diffusion coefficient_ @f$\mathbf{\alpha}@f$, see also @lref{ex:rdemp}
 *
 * ## Template parameter requirement
 *
 * - SCALAR must be a type like `double`
 * - DIFF_COEFF must model the concept of a \ref mesh_function "MeshFunction"
 *   that returns either scalars or matrices.
 *
 *
 */
template <typename SCALAR, typename DIFF_COEFF>
class DiffusionElementMatrixProvider final {
  static_assert(mesh::utils::isMeshFunction<DIFF_COEFF>);

 public:
  /**
   * @brief type of returned element matrix
   */
  using Scalar =
      typename decltype(mesh::utils::MeshFunctionReturnType<DIFF_COEFF>() *
                        Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>())::Scalar;
  using ElemMat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

  /** @brief standard constructors */
  /** @{ */
  DiffusionElementMatrixProvider(const DiffusionElementMatrixProvider &) =
      delete;
  DiffusionElementMatrixProvider(DiffusionElementMatrixProvider &&) noexcept =
      default;
  DiffusionElementMatrixProvider &operator=(
      const DiffusionElementMatrixProvider &) = delete;
  DiffusionElementMatrixProvider &operator=(DiffusionElementMatrixProvider &&) =
      delete;
  /** @} */

  /**
   * @brief Constructor: cell-independent precomputations
   *
   * @param fe_space collection of specifications for scalar-valued parametric
   * reference elements
   * @param alpha mesh function for the (possibly matrix-valued) diffusion
   * coefficient
   *
   * @see LocCompLagrFEPreprocessor::LocCompLagrFEPreprocessor()
   *
   * This constructor uses local quadature rules with double the degree of
   * exactness as the polynomial degree of the finite element space.
   */
  DiffusionElementMatrixProvider(
      std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space, DIFF_COEFF alpha);

  /**
   * @brief All cells are considered active.
   */
  bool isActive(const lf::mesh::Entity & /*cell*/) const { return true; }
  /**
   * @brief main routine for the computation of element matrices
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a small dense, containing the element matrix.
   *
   * Actual computation of the element matrix based on numerical quadrature and
   * mapping techniques. The order of the quadrature rule is tied to the
   * polynomial degree of the underlying finite element spaces: for
   * polynomial degree p a quadrature rule is chosen that is exact for
   * polynomials of degree 2p.
   *
   * Throws an assertion in case the finite element specification is missing for
   * the type of the cell.
   */
  ElemMat Eval(const lf::mesh::Entity &cell) const;

  /** destructor */
  ~DiffusionElementMatrixProvider() = default;

 private:
  /** @name functors providing coefficient functions
   * @{ */
  /** Diffusion coefficient */
  DIFF_COEFF alpha_;
  /** @} */

  /** FE Space */
  std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space_;

  quad::QuadRuleCache qr_cache_;
};

/**
 * @brief logger for DiffusionElementMatrixProvider
 */
extern std::shared_ptr<spdlog::logger> diffusion_element_matrix_provider_logger;

template <class PTR, class DIFF_COEFF>
DiffusionElementMatrixProvider(PTR fe_space, DIFF_COEFF alpha)
    -> DiffusionElementMatrixProvider<typename PTR::element_type::Scalar,
                                      DIFF_COEFF>;

// First constructor (internal construction of quadrature rules
template <typename SCALAR, typename DIFF_COEFF>
DiffusionElementMatrixProvider<SCALAR, DIFF_COEFF>::
    DiffusionElementMatrixProvider(
        std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space, DIFF_COEFF alpha)
    : alpha_(std::move(alpha)), fe_space_(std::move(fe_space)) {}

// TODO(craffael) remove const once
// https://developercommunity.visualstudio.com/content/problem/180948/vs2017-155-c-cv-qualifiers-lost-on-type-alias-used.html
// is resolved
template <typename SCALAR, typename DIFF_COEFF>
typename lf::fe::DiffusionElementMatrixProvider<SCALAR, DIFF_COEFF>::ElemMat
DiffusionElementMatrixProvider<SCALAR, DIFF_COEFF>::Eval(
    const lf::mesh::Entity &cell) const {
  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  SPDLOG_LOGGER_TRACE(diffusion_element_matrix_provider_logger,
                      "{}, shape = \n{}", cell.RefEl(),
                      geo_ptr->Global(cell.RefEl().NodeCoords()));
  // Physical dimension of the cell
  const dim_t world_dim = geo_ptr->DimGlobal();

  // Get a quadrature rule of sufficiently high degree on the element
  const auto sfl = fe_space_->ShapeFunctionLayout(cell);
  const lf::quad::QuadRule qr = qr_cache_.Get(cell.RefEl(), 2 * sfl->Degree());

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
  // Fetch values of diffusion coefficient alpha at quadrature points:
  auto alphaval = alpha_(cell, qr.Points());
  // The requested element matrix is square, size = number of local shape
  // functions
  ElemMat mat(sfl->NumRefShapeFunctions(), sfl->NumRefShapeFunctions());
  mat.setZero();
  // Compute the gradients of the reference shape functions in the quadrature
  // points
  const auto grsf = sfl->GradientsReferenceShapeFunctions(qr.Points());
  // Quadrature formula: loop over quadrature points
  for (base::size_type k = 0; k < qr.NumPoints(); ++k) {
    const double w = qr.Weights()[k] * determinants[k];
    // Transformed gradient in current quadrature node
    const auto trf_grad(JinvT.block(0, 2 * k, world_dim, 2) *
                        grsf.block(0, 2 * k, mat.rows(), 2).transpose());
    // Transformed gradients multiplied with coefficient
    mat += w * trf_grad.adjoint() * (alphaval[k] * trf_grad);
  }
  return mat;
}

/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/fe/fe.h
 * @brief Class for local quadrature based computation of element matrix for
 * Lagrangian finite elements and a weighted \f$L^2\f$ inner product.
 *
 * @tparam SCALAR scalar type of the FiniteElementSpace. Must be a field
 *                     type such as `double` or `std::complex<double>`
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
    (u,v)
 \mapsto\int\limits_{K}\gamma(\mathbf{x})u\,\overline{v}\,\mathrm{d}\mathbf{x}
 \;,
 * @f]
 * with reaction coefficient @f$\gamma@f$, see also @lref{ex:rdemp}
 *
 * ## Template parameter requirement
 *
 * - SCALAR must be a type like `double`
 * - REACTION_COEFF should be compatible with a scalar-valued \ref mesh_function
 *
 */
template <typename SCALAR, typename REACTION_COEFF>
class MassElementMatrixProvider final {
  static_assert(mesh::utils::isMeshFunction<REACTION_COEFF>);

 public:
  /**
   * @brief type of returned element matrix
   */
  using Scalar =
      typename decltype(mesh::utils::MeshFunctionReturnType<REACTION_COEFF>() *
                        Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>())::Scalar;
  using ElemMat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

  /** @brief standard constructors */
  /** @{ */
  MassElementMatrixProvider(const MassElementMatrixProvider &) = delete;
  MassElementMatrixProvider(MassElementMatrixProvider &&) noexcept = default;
  MassElementMatrixProvider &operator=(const MassElementMatrixProvider &) =
      delete;
  MassElementMatrixProvider &operator=(MassElementMatrixProvider &&) = delete;
  /** @} */

  /**
   * @brief Constructor: cell-independent precomputations
   *
   * @param fe_space collection of specifications for scalar-valued parametric
   * reference elements
   * @param gamma mesh function providing scalar-valued diffusion coefficient
   *
   * @see LocCompLagrFEPreprocessor::LocCompLagrFEPreprocessor()
   *
   * This constructor uses local quadature rules with double the degree of
   * exactness as the polynomial degree of the finite element space.
   */
  MassElementMatrixProvider(
      std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space,
      REACTION_COEFF gamma);

  /**
   * @brief All cells are considered active.
   *
   */
  bool isActive(const lf::mesh::Entity & /*cell*/) const { return true; }
  /**
   * @brief main routine for the computation of element matrices
   *
   * @param cell reference to the (triangular or quadrilateral) cell for
   *        which the element matrix should be computed.
   * @return a small dense, containing the element matrix.
   *
   * Actual computation of the element matrix based on numerical quadrature and
   * mapping techniques. The order of the quadrature rule is tied to the
   * polynomial degree of the underlying finite element spaces: for
   * polynomial degree p a quadrature rule is chosen that is exact for
   * polynomials o degree 2p.
   *
   * Throws an assertion in case the finite element specification is missing for
   * the type of the cell.
   */
  ElemMat Eval(const lf::mesh::Entity &cell) const;

  /** destructor */
  ~MassElementMatrixProvider() = default;

 private:
  /** @name functors providing coefficient functions
   * @{ */
  /** Reaction coefficient */
  REACTION_COEFF gamma_;
  /** @} */

  /** FE Space */
  std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space_;

  quad::QuadRuleCache qr_cache_;
};

/**
 * @brief logger for MassElementMatrixProvider
 */
extern std::shared_ptr<spdlog::logger> mass_element_matrix_provider_logger;

template <class PTR, class REACTION_COEFF>
MassElementMatrixProvider(PTR fe_space, REACTION_COEFF gamma)
    -> MassElementMatrixProvider<typename PTR::element_type::Scalar,
                                 REACTION_COEFF>;

// First constructor
template <typename SCALAR, typename REACTION_COEFF>
MassElementMatrixProvider<SCALAR, REACTION_COEFF>::MassElementMatrixProvider(
    std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space, REACTION_COEFF gamma)
    : gamma_(std::move(gamma)), fe_space_(std::move(fe_space)) {}

// TODO(craffael) remove const once
// https://developercommunity.visualstudio.com/content/problem/180948/vs2017-155-c-cv-qualifiers-lost-on-type-alias-used.html
// is resolved
template <typename SCALAR, typename REACTION_COEFF>
typename lf::fe::MassElementMatrixProvider<SCALAR, REACTION_COEFF>::ElemMat
MassElementMatrixProvider<SCALAR, REACTION_COEFF>::Eval(
    const lf::mesh::Entity &cell) const {
  // Query the shape of the cell
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG((geo_ptr->DimLocal() == 2),
                "Only 2D implementation available!");
  SPDLOG_LOGGER_TRACE(mass_element_matrix_provider_logger, "{}, shape = \n{}",
                      cell.RefEl(), geo_ptr->Global(cell.RefEl().NodeCoords()));
  // Physical dimension of the cell
  const dim_t world_dim = geo_ptr->DimGlobal();

  // Get a quadrature rule of sufficiently high degree on the element
  const auto sfl = fe_space_->ShapeFunctionLayout(cell);
  const lf::quad::QuadRule qr = qr_cache_.Get(cell.RefEl(), 2 * sfl->Degree());
  // Metric factors in quadrature nodes
  const Eigen::VectorXd determinants(geo_ptr->IntegrationElement(qr.Points()));
  LF_ASSERT_MSG(
      determinants.size() == qr.NumPoints(),
      "Mismatch " << determinants.size() << " <-> " << qr.NumPoints());
  // Fetch values of reaction coefficient gamma at quadrature points:
  auto gammaval = gamma_(cell, qr.Points());
  // Alocated element matrix
  ElemMat mat(sfl->NumRefShapeFunctions(), sfl->NumRefShapeFunctions());
  mat.setZero();
  // Compute the reference shape functions in the quadrature points
  const auto rsf = sfl->EvalReferenceShapeFunctions(qr.Points());
  // Loop over quadrature points to evaluate quadrature formula
  for (base::size_type k = 0; k < qr.NumPoints(); ++k) {
    const double w = qr.Weights()[k] * determinants[k];
    mat += w * ((gammaval[k] * rsf.col(k)) * (rsf.col(k).adjoint()));
  }
  return mat;
}

/**
 * @ingroup entity_matrix_provider
 * @headerfile lf/fe/fe.h
 * @brief Quadrature-based computation of local mass matrix for an edge
 *
 * @tparam SCALAR underlying scalar type of the ScalarFESpace, usually double or
 * complex<double>
 * @tparam COEFF \ref mesh_function "MeshFunction" that defines the
 * scalar valued coefficient \f$ \gamma \f$
 * @tparam EDGESELECTOR predicate defining which edges are included
 *
 * This \ref entity_matrix_provider "EntityMatrixProvider" class corresponds to
 * the the element matrix for the bilinear form
 * @f[
 *     (u,v) \mapsto \int\limits_e
 * \gamma(x)u(x)\overline{v(x)}\,\mathrm{d}S(x)\;,
 * @f]
 * where @f$e@f$ is an edge of the mesh, and @f$\gamma@f$ a scalar-valued
 * coefficient function.
 *
 */
template <typename SCALAR, typename COEFF, typename EDGESELECTOR>
class MassEdgeMatrixProvider final {
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
   * exactness as the polynomial degree of the finite element.
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
  bool isActive(const lf::mesh::Entity &edge) const {
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
  ElemMat Eval(const lf::mesh::Entity &edge) const;

  ~MassEdgeMatrixProvider() = default;

 private:
  COEFF gamma_;            // functor for coefficient
  EDGESELECTOR edge_sel_;  // Defines the active edges
  // FE Space
  std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space_;

  quad::QuadRuleCache qr_cache_;
};

/**
 * @brief logger for MassEdgeMatrixProvider
 */
extern std::shared_ptr<spdlog::logger> mass_edge_matrix_provider_logger;

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
    const lf::mesh::Entity &edge) const {
  const lf::geometry::Geometry *geo_ptr = edge.Geometry();
  // Get the shape function layout on the edge
  const auto sfl = fe_space_->ShapeFunctionLayout(edge);
  // Compute a quadrature rule on the given entity
  const lf::quad::QuadRule qr = qr_cache_.Get(edge.RefEl(), 2 * sfl->Degree());
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
    mat += ((rsf.col(k)) * (rsf.col(k).adjoint())) * w;
  }
  return mat;
}

/**
 * @ingroup entity_vector_provider
 * @headerfile lf/fe/fe.h
 * @brief Local computation of general element (load) vector for scalar
 finite
 * elements; volume contributions only
 *
 * @tparam SCALAR underlying scalar type of the ScalarFESpace, usually double or
 *                complex<double>
 * @tparam MESH_FUNCTION \ref mesh_function "MeshFunction" which defines the
 *                       source function \f$ f \f$
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
 * `ELEM_VEC_COMP` of the function assemble::AssembleVectorLocally().
 */
template <typename SCALAR, typename MESH_FUNCTION>
class ScalarLoadElementVectorProvider final {
  static_assert(mesh::utils::isMeshFunction<MESH_FUNCTION>);

 public:
  using scalar_t = decltype(
      SCALAR(0) * mesh::utils::MeshFunctionReturnType<MESH_FUNCTION>(0));
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
      std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space, MESH_FUNCTION f);
  /** @brief all cells are active */
  bool isActive(const lf::mesh::Entity & /*cell*/) const { return true; }
  /*
   * @brief Main method for computing the element vector
   *
   * @param cell current cell for which the element vector is desired
   * @return local load vector as column vector
   *
   */
  ElemVec Eval(const lf::mesh::Entity &cell) const;

  ~ScalarLoadElementVectorProvider() = default;

 private:
  /** @brief An object providing the source function */
  MESH_FUNCTION f_;

  std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space_;

  quad::QuadRuleCache qr_cache_;
};

/**
 * @brief logger used by ScalarLoadElementVectorProvider
 */
extern std::shared_ptr<spdlog::logger>
    scalar_load_element_vector_provider_logger;

// Deduction guide
template <class PTR, class MESH_FUNCTION>
ScalarLoadElementVectorProvider(PTR fe_space, MESH_FUNCTION mf)
    -> ScalarLoadElementVectorProvider<typename PTR::element_type::Scalar,
                                       MESH_FUNCTION>;

// Constructors
template <typename SCALAR, typename MESH_FUNCTION>
ScalarLoadElementVectorProvider<SCALAR, MESH_FUNCTION>::
    ScalarLoadElementVectorProvider(
        std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space, MESH_FUNCTION f)
    : f_(std::move(f)), fe_space_(std::move(fe_space)) {}

// TODO(craffael) remove const once
// http://developercommunity.visualstudio.com/content/problem/180948/vs2017-155-c-cv-qualifiers-lost-on-type-alias-used.html
// is resolved
template <typename SCALAR, typename MESH_FUNCTION>
typename ScalarLoadElementVectorProvider<SCALAR, MESH_FUNCTION>::ElemVec
ScalarLoadElementVectorProvider<SCALAR, MESH_FUNCTION>::Eval(
    const lf::mesh::Entity &cell) const {
  // Type for source function
  using source_fn_t = mesh::utils::MeshFunctionReturnType<MESH_FUNCTION>;

  // Get the shape function layout for the given cell
  const auto sfl = fe_space_->ShapeFunctionLayout(cell);

  // Initialize a quadrature rule of sufficiently high degree
  const lf::quad::QuadRule qr = qr_cache_.Get(cell.RefEl(), 2 * sfl->Degree());

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
    vec +=
        (qr.Weights()[k] * determinants[k] * fval[k]) * rsf.col(k).conjugate();
  }

  SPDLOG_LOGGER_TRACE(scalar_load_element_vector_provider_logger,
                      "LOCVEC = \n{}", vec.transpose());
  return vec;
}

/**
 * @ingroup entity_vector_provider
 * @headerfile lf/fe/fe.h
 * @brief Local edge contributions to element vector
 *
 * @tparam SCALAR underlying scalar type of the FESpace, usually double or
 *                complex<double>
 * @tparam FUNCTOR `SCALAR` valued \ref mesh_function "MeshFunction" which
 *                 defines the function \f$ g \f$
 * @tparam EDGESELECTOR selector type for active edges
 *
 * The underlying local linear form for an edge @f$e@f$ is
 * @f[
    v \mapsto \int_e
 g(\mathbf{x})\,\overline{v(\mathbf{x})}\,\mathrm{d}S\mathbf{x}\;,
 * @f]
 * where \f$g\f$ is supposed to be a locally continuous source function.
 *
 * Computations are based on quadrature rules supplied by the LehrFEM++
 * lf::quad::make_QuadRule() method.
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
class ScalarLoadEdgeVectorProvider final {
 public:
  static_assert(mesh::utils::isMeshFunction<FUNCTOR>,
                "FUNCTOR does not fulfill the concept of a mesh function.");
  using Scalar =
      decltype(SCALAR(0) * mesh::utils::MeshFunctionReturnType<FUNCTOR>(0));
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

  /** @brief Constructor
   *
   * @param fe_space ScalarFESpace that supplied the finite elements.
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

  /** @brief all edges are active */
  bool isActive(const lf::mesh::Entity &cell) const { return edge_sel_(cell); }
  /*
   * @brief Main method for computing the element vector
   *
   * @param cell current cell for which the element vector is desired
   * @return local load vector as column vector
   *
   */
  ElemVec Eval(const lf::mesh::Entity &edge) const;

  ~ScalarLoadEdgeVectorProvider() = default;

 private:
  FUNCTOR g_;              // source function
  EDGESELECTOR edge_sel_;  // selects edges
  std::shared_ptr<const ScalarFESpace<SCALAR>> fe_space_;

  quad::QuadRuleCache qr_cache_;
};

/**
 * @brief logger for ScalarLoadEdgeVectorProvider class template.
 */
extern std::shared_ptr<spdlog::logger> scalar_load_edge_vector_provider_logger;

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
    const lf::mesh::Entity &edge) const {
  // Query the shape of the edge
  const lf::geometry::Geometry *geo_ptr = edge.Geometry();
  LF_ASSERT_MSG(geo_ptr != nullptr, "Invalid geometry!");
  LF_ASSERT_MSG(geo_ptr->DimLocal() == 1, "The passed entity is not an edge!");

  // Get the shape function layout of the given edge
  const auto sfl = fe_space_->ShapeFunctionLayout(edge);

  // Quadrature points on physical edge
  const lf::quad::QuadRule qr = qr_cache_.Get(edge.RefEl(), 2 * sfl->Degree());

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
    vec += rsf.col(k).conjugate() * w;
  }
  return vec;
}

}  // namespace lf::fe

#endif
