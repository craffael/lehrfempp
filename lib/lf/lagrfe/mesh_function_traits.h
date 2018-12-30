/**
 * @file
 * @brief Defines a compile time constant that can be used to check if a given
 *        type fulfills the MeshFunction concept + a way to determine the
 *        type of objects returned by a mesh function.
 * @author Raffael Casagrande
 * @date   2018-12-15 03:49:01
 * @copyright MIT License
 */

#ifndef __7e4ffaa81e244723acbfbbaea68e03b1
#define __7e4ffaa81e244723acbfbbaea68e03b1

#include <lf/mesh/mesh.h>
#include <type_traits>

namespace lf::lagrfe {
namespace internal {
template <class T>
using MeshFunctionReturnType_t =
    std::invoke_result_t<T, const lf::mesh::Entity, const Eigen::MatrixXd>;

template <class T>
auto getVectorType(const std::vector<T>& a, int) -> T {
  return a[0];
}

template <class T>
void getVectorType(const T& a, long) {}

template <class T>
using VectorElement_t = decltype(getVectorType(std::declval<T>(), 0));

template <class T, class RETURN_TYPE,
          class = typename std::enable_if<!std::is_same_v<
              VectorElement_t<MeshFunctionReturnType_t<T>>, void>>::type>
constexpr bool IsMeshFunctionCallable(int) {
  if constexpr (std::is_same_v<RETURN_TYPE, void>) {
    // user didn't want us to check whether the return type is something
    // particular
    return true;
  } else {
    // user specified a RETURN_TYPE -> Check that std::vector<RETURN_TYPE> is
    // returned.
    return std::is_same_v<VectorElement_t<MeshFunctionReturnType_t<T>>,
                          RETURN_TYPE>;
  }
}

template <class T, class RETURN_TYPE>
constexpr bool IsMeshFunctionCallable(long) {
  return false;
}

}  // namespace internal

/**
 * @brief Determine whether a given type fulfills the concept \ref
 * mesh_function.
 * @tparam T The type to check
 * @tparam R If specified, check additionally, that the \ref mesh_function
 * returns objects of type `R`
 */
template <class T, class R = void>
constexpr bool isMeshFunction =
    !std::is_reference_v<T> && std::is_copy_constructible_v<T> &&
    std::is_move_constructible_v<T> &&
    internal::IsMeshFunctionCallable<T, R>(0);

/**
 * @brief Determine the type of objects returned by a MeshFunction
 * @tparam T The type of the mesh function.
 */
template <class T>
using MeshFunctionReturnType =
    internal::VectorElement_t<internal::MeshFunctionReturnType_t<T>>;

}  // namespace lf::lagrfe

#endif  // __7e4ffaa81e244723acbfbbaea68e03b1
