/**
 * @file
 * @brief Implementation of basic validation and diagnostics facilities
 * @author Raffael Casagrande
 * @date   2018-06-29 09:18:48
 * @copyright MIT License
 */

#ifndef c3c605c9e48646758bf03fab65d52836
#define c3c605c9e48646758bf03fab65d52836

#include <boost/assert.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lf::base {

void AssertionFailed(const std::string& expr, const std::string& file,
                     long line, const std::string& msg);

}  // namespace lf::base

/**
 * @brief `LF_UNREACHABLE` is the same as `LF_VERIFY(false)` but its easier
 * for the compiler to analyze and helps to avoid warnings about "control
 * reached end of non-void function"
 */
#define LF_UNREACHABLE                                              \
  {                                                                 \
    ::lf::base::AssertionFailed("false", __FILE__, __LINE__, ""); \
    std::abort();                                                   \
  }

/**
 * @brief LF_VERIFY_MSG(expr, msg) aborts execution of the code if
 * `expr` evaluates to false.
 *
 * @Remark At the moment this macro just forwards everything to
 * BOOST_VERIFY_MSG() and throws an exception if the condition is not satisfied
 * in order to help IDE tools detect that
 * `BOOST_VERIFY_MSG(false, "message")` always aborts execution.
 */
#define LF_VERIFY_MSG(expr, msg)                                        \
  {                                                                     \
    if (!(expr)) {                                                      \
      std::stringstream ss;                                             \
      ss << msg; /* NOLINT */                                           \
      ::lf::base::AssertionFailed(#expr, __FILE__, __LINE__, ss.str()); \
      throw std::runtime_error("this code should not be reached");      \
    }                                                                   \
  }

#ifdef NDEBUG
#define LF_ASSERT_MSG_CONSTEXPR(expr, msg) ((void)0)
#define LF_ASSERT_MSG(expr, msg) ((void)0)
#else
#define LF_ASSERT_MSG_CONSTEXPR(expr, msg)      \
  {                                             \
    if (!(expr)) throw std::runtime_error(msg); \
  }

#define LF_ASSERT_MSG(expr, msg)                                        \
  {                                                                     \
    if (!(expr)) {                                                      \
      std::stringstream ss;                                             \
      ss << msg; /* NOLINT */                                           \
      ::lf::base::AssertionFailed(#expr, __FILE__, __LINE__, ss.str()); \
      LF_UNREACHABLE;                                                 \
    }                                                                   \
  }
#endif

// And now we will redefine eigen_assert if needed:
#ifdef LF_REDIRECT_ASSERTS
#ifdef eigen_assert
#ifndef eigen_assert_redirected
// eigen_assert has already been defined, but not by us!
#ifdef _MSC_VER
#pragma message( \
    "WARNING: Eigen has been included before all LehrFEM++ headers but LF_REDIRECT_ASSERTS=On in cmake. Not all Eigen Asserts may print a stacktrace! https://craffael.github.io/lehrfempp/eigen_stacktrace_warning.html")
#else
#warning \
    "Eigen has been included before all LehrFEM++ headers but LF_REDIRECT_ASSERTS=On in cmake. Not all Eigen Asserts may print a stacktrace! https://craffael.github.io/lehrfempp/eigen_stacktrace_warning.html"
#endif
#undef eigen_assert
#define eigen_assert(x) LF_ASSERT_MSG(x, "")
#define eigen_assert_redirected
#endif
#else
#define eigen_assert(x) LF_ASSERT_MSG(x, "")
#define eigen_assert_redirected
#endif
#endif
#endif  // c3c605c9e48646758bf03fab65d52836
