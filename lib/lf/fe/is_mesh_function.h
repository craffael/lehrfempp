/**
 * @file
 * @brief Defines a compile time constant that can be used to check if a given
 *        type fulfills the MeshFunction concept.
 * @author Raffael Casagrande
 * @date   2018-12-15 03:49:01
 * @copyright MIT License
 */

#ifndef __7e4ffaa81e244723acbfbbaea68e03b1
#define __7e4ffaa81e244723acbfbbaea68e03b1

#include <lf/mesh/mesh.h>
#include <type_traits>

namespace lf::fe {
namespace internal {
template <class T>
using MeshFunctionReturnType_t = decltype(std::declval<T>()(
    std::declval<const lf::mesh::Entity>(), std::declval<Eigen::MatrixXd>()));

template <
    class T,
    class = typename std::enable_if<
        std::is_convertible_v<MeshFunctionReturnType_t<T>, Eigen::MatrixXd> ||
        std::is_convertible_v<MeshFunctionReturnType_t<T>,
                              Eigen::MatrixXcd>>::type>
constexpr bool IsMeshFunctionCallable(int) {
  return true;
}

template <class T>
constexpr bool IsMeshFunctionCallable(long) {
  return false;
}

}  // namespace internal

template <class T>
constexpr bool isMeshFunction =
    !std::is_reference_v<T> && std::is_copy_constructible_v<T> &&
    std::is_move_constructible_v<T> && internal::IsMeshFunctionCallable<T>(0);

}  // namespace lf::fe

#endif  // __7e4ffaa81e244723acbfbbaea68e03b1
