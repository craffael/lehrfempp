/**
 * @file
 * @brief Defines a compile time constant that can be used to check if a given
 *        type fulfills the MeshFunction concept + a way to determine the
 *        type of objects returned by a mesh function.
 * @author Raffael Casagrande
 * @date   2018-12-15 03:49:01
 * @copyright MIT License
 */

#ifndef INCG7e4ffaa81e244723acbfbbaea68e03b1
#define INCG7e4ffaa81e244723acbfbbaea68e03b1

#include <lf/mesh/mesh.h>

#include <concepts>
#include <type_traits>

namespace lf::mesh::utils {
namespace internal {

// TODO(craffael): Replace this with std::invoke_result_t once issue in msvc is
// fixed: https://github.com/microsoft/STL/issues/1288
template <class T>
using MeshFunctionReturnType_t = decltype(std::declval<T>()(
    std::declval<const Entity>(), std::declval<const Eigen::MatrixXd>()));

template <typename>
struct VectorElementType {
  using type = void;
};

template <typename T, typename A>
struct VectorElementType<std::vector<T, A>> {
  using type = T;
};

template <class T>
using VectorElement_t = typename VectorElementType<T>::type;

template <class T, class RETURN_TYPE>
  requires(!std::is_same_v<VectorElement_t<MeshFunctionReturnType_t<T>>, void>)
constexpr bool IsMeshFunctionCallable(int /*unused*/) {
  if constexpr (std::is_same_v<RETURN_TYPE, void>) {
    // user didn't want us to check whether the return type is something
    // particular
    return true;
  }
  // user specified a RETURN_TYPE -> Check that std::vector<RETURN_TYPE> is
  // returned.
  return std::is_same_v<VectorElement_t<MeshFunctionReturnType_t<T>>,
                        RETURN_TYPE>;
}

template <class T, class RETURN_TYPE>
constexpr bool IsMeshFunctionCallable(long /*unused*/) {
  return false;
}

}  // namespace internal

/**
 * @brief A MeshFunction is a function object that can be evaluated at any point
 * on the mesh.
 * @tparam MF The type which should fulfill the MeshFunction concept.
 * @tparam R If specified, check additionally, that the \ref mesh_function
 * returns objects of type `R`. If `R` is `void`, the MeshFunction can return
 * any type.
 *
 * ## Description
 *
 * Conceptually, a mesh function assigns to every point on the mesh an object of
 * type `R` (e.g. `double` or an `Eigen::Matrix2d`).
 *
 * For efficiency reasons, a mesh function is normally evaluated at a number of
 * points at once. Hence a mesh function must overload the bracket operator as
 * follows:
 * ```
 * std::vector<R> operator()(
 *   const lf::mesh::Entity& e, const Eigen::MatrixXd& local) const
 * ```
 * Here
 * - `e` is a \ref lf::mesh::Entity "mesh entity"
 * - `local` is a matrix of size `(e.RefEl().Dimension()) x NumPoints` and lists
 * the local coordinates of the evaluation points (with respect to the reference
 * element `e.RefEl()`).
 * - `R` is more or less arbitrary type such as a `double` or a
 * `Eigen::Matrix2d`
 *
 * The return type of `operator()` is a `std::vector<R>` with `NumPoints`
 * length.
 *
 * ## Requirements
 *
 * The type `MF` satisfies the concept `MeshFunction` if
 *
 * Given
 * - `R`, the type of objects returned by the mesh function
 * - `e`, a mesh entity of type `lf::mesh::Entity`,
 * - `local`, a set of local coordinates of type `Eigen::MatrixXd`
 * - `mf` object of type `const MF`
 *
 * the following expressions are valid:
 *
 * <table>
 * <tr>
 *   <th>expression
 *   <th>return type
 *   <th>semantics
 * <tr>
 *   <td> `MF(a)`
 *   <td> `MF`
 *   <td> Creates a copy of `a`
 * <tr>
 *   <td> `MF(std::move(a))`
 *   <td> `MF`
 *   <td> "steals" `a` to create a new MeshFunction
 * <tr>
 *   <td> `a(e, local)`
 *   <td> `std::vector<R>`
 *   <td>Evaluates mesh function at points `local` on the entity `e`
 * </table>
 *
 * ## Usage scenarios
 * The concept of a MeshFunction is used widely in the `lf::uscalfe` module:
 * - Assembler classes such as lf::uscalfe::ScalarLoadEdgeVectorProvider or
 * lf::uscalfe::ReactionDiffusionElementMatrixProvider accept MeshFunctions that
 * describe coefficients or source functions.
 *
 * ## Archetype
 * - lf::mesh::utils::MeshFunctionAT
 *
 * ## See also
 * - lf::mesh::utils::MeshFunctionReturnType to retrieve the type `R` of a mesh
 * function.
 *
 * ## Classes modelling the MeshFunction concept
 * - lf::fe::MeshFunctionFE
 * - lf::fe::MeshFunctionGradFE
 * - lf::mesh::utils::MeshFunctionBinary with it's operator overloads
 * - lf::mesh::utils::MeshFunctionConstant
 * - lf::mesh::utils::MeshFunctionGlobal
 * - lf::mesh::utils::MeshFunctionUnary with it's operator overloads
 * - lf::refinement::MeshFunctionTransfer
 *
 *
 */
template <class MF, class R = void>
concept MeshFunction =
    // @note We use IsMeshFunctionCallable here because otherwise e.g.
    // Eigen::Matrix also fulfills the concept of a MeshFunction (even though it
    // leads to compile errors during instantiation)
    std::is_object_v<MF> && std::copy_constructible<MF> &&
    internal::IsMeshFunctionCallable<MF, R>(0);

/**
 * @brief Archetype for the MeshFunction concept
 * @tparam R The return type of the mesh function, i.e. the value it takes at
 * every point on the mesh.
 *
 * @sa \ref archetypes for more information about archetypes.
 */
template <class R>
struct MeshFunctionAT {
  /// Copy constructor
  MeshFunctionAT(const MeshFunctionAT&) noexcept = default;

  /// Move constructor
  MeshFunctionAT(MeshFunctionAT&&) noexcept = default;

  /// Copy assignment operator deleted because not part of the concept
  auto operator=(const MeshFunctionAT&) noexcept -> MeshFunctionAT& = delete;

  /// Move assignment operator deleted because not part of the concept
  auto operator=(MeshFunctionAT&&) noexcept -> MeshFunctionAT& = delete;

  /// Destructor
  ~MeshFunctionAT() noexcept = default;

  /**
   * @brief Evaluates the mesh function at the given points.
   * @param e The entity on which the mesh function should be evaluated.
   * @param local The local coordinates of the points at which the mesh function
   * should be evaluated. Has dimension `e.RefEl().Dimension() x NumPoints`.
   * @return A vector of length `local.cols()` containing the values of the mesh
   * function at the given points.
   */
  std::vector<R> operator()(const Entity& /*e*/,
                            const Eigen::MatrixXd& /*local*/) const {
    LF_VERIFY_MSG(false, "Should never be called, this is an archetype");
    return {};
  }
};

/**
 * @brief Determine the type of objects returned by a MeshFunction
 * @tparam R The type of the mesh function (as defined in \ref mesh_function).
 */
template <class R>
using MeshFunctionReturnType =
    internal::VectorElement_t<internal::MeshFunctionReturnType_t<R>>;

}  // namespace lf::mesh::utils

#endif  // INCG7e4ffaa81e244723acbfbbaea68e03b1
