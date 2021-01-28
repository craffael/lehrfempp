/**
 * @file
 * @brief meta-programming utility to determine if a given type is a "scalar
 * type"
 * @author Raffael Casagrande
 * @date   2021-01-01 03:45:04
 * @copyright MIT License
 */

#ifndef __e384f4f7485d447c8d3f33406a0f5f7e
#define __e384f4f7485d447c8d3f33406a0f5f7e

#include <complex>
#include <type_traits>

namespace lf::base {

/**
 * @brief Base Traits class which can be used to determine if a type `T` is a
 * scalar value.
 *
 * This class template is specialized below for types that are actually scalars.
 * @sa is_scalar<T>
 */
template <class T, typename = void>
struct IsScalar {
  static constexpr bool value = false;
};

/**
 * @brief any type that is a `std::is_arithmetic` is also a scalar type.
 */
template <class T>
struct IsScalar<T, std::enable_if_t<std::is_arithmetic_v<T>>> {
  static constexpr bool value = true;
};

/**
 * @brief Also `std::complex` is a "scalar" type if the underlying type is a
 * `std::is_arithmetic` type.
 */
template <class T>
struct IsScalar<std::complex<T>, std::enable_if_t<std::is_arithmetic_v<T>>> {
  static constexpr bool value = true;
};

/**
 * @brief Variable template that determines if a type `T` is a scalar type, i.e.
 * if it is a "field" in the mathematical sense.
 * @tparam T The type that should be tested.
 *
 * @note This method is e.g. used to determine if a given type can be multiplied
 * to a Eigen::Matrix.
 */
template <class T>
inline constexpr bool is_scalar = IsScalar<T>::value;

}  // namespace lf::base

#endif  // __e384f4f7485d447c8d3f33406a0f5f7e
