/**
 * @file
 * @brief Utilities to deal with Eigen related peculiarities
 * @author Raffael Casagrande
 * @date   2019-01-19 07:21:07
 * @copyright MIT License
 */

#ifndef INCG1cc1076600024d7ea537871be7fc1fc0
#define INCG1cc1076600024d7ea537871be7fc1fc0

// clang-format off
#include "lf_assert.h"  // must be included before eigen!
// clang-format on

#include <Eigen/Core>
#include <utility>
#include <spdlog/fmt/ranges.h>

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
inline constexpr bool is_eigen_matrix = std::is_same_v<
    decltype(internal::IsEigenMatrixTester::Test(std::declval<T>(), 0)), bool>;

/**
 * @brief Check if a given type T is an Eigen::Array
 */
template <class T>
inline constexpr bool is_eigen_array = std::is_same_v<
    decltype(internal::IsEigenArrayTester::Test(std::declval<T>(), 0)), bool>;

}  // namespace lf::base

// The following is needed to prohibit fmt to treat Eigen matrices/arrays as
// ranges.
template <class MATRIX>
  requires(std::is_base_of_v<Eigen::DenseBase<MATRIX>, MATRIX>)
struct fmt::is_range<MATRIX, char> {
  static constexpr const bool value = false;
};

/**
 * @brief this is the fmt::formatter which is used to format eigen
 * matrices/arrays. You can change the precision using a format string such as
 * `{.5}` to get 5 significant digits.
 */
template <class MATRIX>
  requires(std::is_base_of_v<Eigen::DenseBase<MATRIX>, MATRIX>)
struct fmt::formatter<MATRIX> {
  constexpr auto parse(const format_parse_context& ctx) {
    const auto* it = ctx.begin();
    const auto* end = ctx.end();

    if (it != end && *it == '.') {
      // Support precision:
      ++it;
      auto start = it;
      while (is_digit(*it) && it != end) {
        ++it;
      }
      precision = parse_unsigned_int(start, it);
    } else {
      precision = 4;
    }

    if (it != end && *it != '}') {
      throw format_error("invalid format");
    }

    return it;
  }

  template <typename FormatContext>
  auto format(const MATRIX& m, FormatContext& ctx) const {
    std::stringstream ss;
    Eigen::IOFormat format(precision, 0, ", ", "\n", "[", "]");
    ss << m.format(format);

    auto it = ctx.out();
    auto str = ss.str();
    std::copy(str.begin(), str.end(), it);
    return it;
  }

 private:
  int precision = 4;

  constexpr bool is_digit(char c) { return c <= '9' && c >= '0'; }
  constexpr int parse_unsigned_int(const char* begin, const char* end) {
    int result = 0;
    while (begin != end) {
      if (!is_digit(*begin)) {
        throw "compile time error, this is not a digit.";
      }
      result += (*begin - '0');
      result *= 10;
      ++begin;
    }
    result = result / 10;
    return result;
  }
};

#endif  // INCG1cc1076600024d7ea537871be7fc1fc0
