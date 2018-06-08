

#ifndef __c3c605c9e48646758bf03fab65d52836
#define __c3c605c9e48646758bf03fab65d52836

#include <sstream>
#include <stdexcept>
#include <string>
#include "boost/assert.hpp"

/**
 * @brief HL_VERIFY_MSG(expr, msg) aborts execution of the code if
 * `expr` evaluates to false.
 *
 * @Remark At the moment this macro just forwards everything to
 * BOOST_VERIFY_MSG() and throws an exception if the condition is not satisfied
 * in order to help IDE tools detect that
 * `BOOST_VERIFY_MSG(false, "message")` always aborts execution.
 */
#define LF_VERIFY_MSG(expr, msg)                                      \
  {                                                                   \
    BOOST_VERIFY_MSG(expr, std::string(msg).c_str());                 \
    if (!(expr))                                                      \
      throw std::runtime_error("This exception should never throw."); \
  }

#ifdef NDEBUG
#define LF_ASSERT_MSG_CONSTEXPR(expr, msg) ((void)0)
#else
#define LF_ASSERT_MSG_CONSTEXPR(expr, msg)      \
  {                                             \
    if (!(expr)) throw std::runtime_error(msg); \
  }
#endif

#define LF_ASSERT_MSG(expr, msg)                \
  {                                             \
    if (!(expr)) {                              \
      std::stringstream ss;                     \
      ss << msg;                                \
      BOOST_ASSERT_MSG(expr, ss.str().c_str()); \
    }                                           \
  }

#endif  // __c3c605c9e48646758bf03fab65d52836
