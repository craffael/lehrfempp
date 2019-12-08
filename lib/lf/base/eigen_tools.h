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

struct IsEigenMatrixTester {
  template <class SCALAR, int ROWS, int COLS, int OPTIONS, int MAX_ROWS,
            int MAX_COLS>
  static bool Test(
      const Eigen::Matrix<SCALAR, ROWS, COLS, OPTIONS, MAX_ROWS, MAX_COLS>&,
      int);

  template <class T>
  static float Test(const T&, long);
};

struct IsEigenArrayTester {
  template <class SCALAR, int ROWS, int COLS, int OPTIONS, int MAX_ROWS,
            int MAX_COLS>
  static bool Test(
      const Eigen::Array<SCALAR, ROWS, COLS, OPTIONS, MAX_ROWS, MAX_COLS>&,
      int);

  template <class T>
  static float Test(const T&, long);
};

}  // namespace internal

/**
 * @brief Check if a given type T is an Eigen::Matrix
 */
template <class T>
constexpr bool is_eigen_matrix = std::is_same_v<
    decltype(internal::IsEigenMatrixTester::Test(std::declval<T>(), 0)), bool>;

/**
 * @brief Check if a given type T is an Eigen::Array
 */
template <class T>
constexpr bool is_eigen_array = std::is_same_v<
    decltype(internal::IsEigenArrayTester::Test(std::declval<T>(), 0)), bool>;

}  // namespace lf::base

#endif  // __1cc1076600024d7ea537871be7fc1fc0
