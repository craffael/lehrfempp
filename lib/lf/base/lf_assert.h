/**
 * @file
 * @brief Implementation of basic validation and diagnostics facilities
 * @author Raffael Casagrande
 * @date   2018-06-29 09:18:48
 * @copyright MIT License
 */

#ifndef __c3c605c9e48646758bf03fab65d52836
#define __c3c605c9e48646758bf03fab65d52836

#include <sstream>
#include <stdexcept>
#include <string>

namespace lf::base {

void AssertionFailed(const std::string& expr, const std::string& file, int line,
                     const std::string& msg);

}  // namespace lf::base

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
      throw std::runtime_error("this code should not be reached");      \
    }                                                                   \
  }
#endif

// the following is needed to redirect BOOST_ASSERT(_MSG)/BOOST_VERIFY (_MSG)
namespace boost {

inline void assertion_failed_msg(char const* expr, char const* msg,
                                 char const* function, char const* file,
                                 long line) {
  lf::base::AssertionFailed(expr, file, line, msg);
}

inline void assertion_failed(char const* expr, char const* function,
                             char const* file, long line) {
  lf::base::AssertionFailed(expr, file, line, "");
}

}  // namespace boost

// And now we will redefine eigen_assert if needed:
#ifdef LF_REDIRECT_ASSERTS
#ifdef eigen_assert
#ifndef eigen_assert_redirected
#warning \
    "Eigen has been included before all LehrFEM++ headers but LF_REDIRECT_ASSERTS=On in cmake. Not all Eigen Asserts may print a stacktrace! https://craffael.github.io/lehrfempp/eigen_stacktrace_warning.html"
#undef eigen_assert
#define eigen_assert(x) LF_ASSERT_MSG(x, "")
#define eigen_assert_redirected
#endif
#else
#define eigen_assert(x) LF_ASSERT_MSG(x, "")
#define eigen_assert_redirected
#endif
#endif
#endif  // __c3c605c9e48646758bf03fab65d52836
