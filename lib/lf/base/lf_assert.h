/** @file lf_assert.h 
  *
  * Implementation of basic validation and diagnostics facilities
  */ 

#ifndef __c3c605c9e48646758bf03fab65d52836
#define __c3c605c9e48646758bf03fab65d52836

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <map>

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
      ss << msg;                                                        \
      ::lf::base::AssertionFailed(#expr, __FILE__, __LINE__, ss.str()); \
      std::abort();                                                     \
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
      ss << msg;                                                        \
      ::lf::base::AssertionFailed(#expr, __FILE__, __LINE__, ss.str()); \
      std::abort();                                                     \
      throw std::runtime_error("this code should not be reached");      \
    }                                                                   \
  }
#endif
  
#endif  // __c3c605c9e48646758bf03fab65d52836
