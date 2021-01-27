/** @file lf_assert.cc */

#include "lf_assert.h"

#include <boost/stacktrace.hpp>
#include <iostream>

namespace lf::base {

// Output for assertions
void AssertionFailed(const std::string& expr, const std::string& file,
                     long line, const std::string& msg) {
  std::cerr << "***** Internal Program Error - assertion (" << expr
            << ") failed:\n"
            << file << '(' << line << "): " << msg << std::endl;
  std::cerr << "Backtrace:\n" << boost::stacktrace::stacktrace() << '\n';
  std::abort();
}

}  // end namespace lf::base

// the following is needed to redirect BOOST_ASSERT(_MSG)/BOOST_VERIFY (_MSG)
void boost::assertion_failed_msg(char const* expr, char const* msg,
                                 char const* /*function*/, char const* file,
                                 long line) {
  lf::base::AssertionFailed(expr, file, line, msg);
}

void boost::assertion_failed(char const* expr, char const* /*function*/,
                             char const* file, long line) {
  lf::base::AssertionFailed(expr, file, line, "");
}
