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
using MeshFunctionReturnType_t =
    decltype(std::declval<T>()(std::declval<const lf::mesh::Entity>(),
                               std::declval<const Eigen::MatrixXd>()));

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

template <class T, class RETURN_TYPE = void>
constexpr bool isMeshFunction =
    !std::is_reference_v<T> && std::is_copy_constructible_v<T> &&
    std::is_move_constructible_v<T> &&
    internal::IsMeshFunctionCallable<T, RETURN_TYPE>(0);

}  // namespace lf::fe

#endif  // __7e4ffaa81e244723acbfbbaea68e03b1
