/**
 * @file
 * @brief Utilities to deal with Eigen related peculiarities
 * @author Raffael Casagrande
 * @date   2019-01-19 07:21:07
 * @copyright MIT License
 */

#ifndef __1cc1076600024d7ea537871be7fc1fc0
#define __1cc1076600024d7ea537871be7fc1fc0

#include <Eigen/Eigen>
#include <utility>

namespace lf::base {

namespace internal {
template <class T>
constexpr T&& constexpr_declval() noexcept;
struct IsEigenTester {
  template <class SCALAR, int ROWS, int COLS, int OPTIONS, int MAX_ROWS,
            int MAX_COLS>
  static bool IsEigen(
      const Eigen::Matrix<SCALAR, ROWS, COLS, OPTIONS, MAX_ROWS, MAX_COLS>&,
      int);

  template <class T>
  static float IsEigen(const T&, long);
};

}  // namespace internal

/**
 * @brief Check if a given type T is an Eigen::Matrix
 */
template <class T>
constexpr bool is_eigen_matrix = std::is_same_v<
    decltype(internal::IsEigenTester::IsEigen(std::declval<T>(), 0)), bool>;

}  // namespace lf::base

#endif  // __1cc1076600024d7ea537871be7fc1fc0
