/** @file lf_assert.cc */

#include "lf_assert.h"

#include <boost/stacktrace.hpp>
#include <iostream>

namespace lf::base {

// Output for assertions
void AssertionFailed(const std::string& expr, const std::string& file, int line,
                     const std::string& msg) {
  std::cerr << "***** Internal Program Error - assertion (" << expr
            << ") failed:\n"
            << file << '(' << line << "): " << msg << std::endl;
  std::cerr << "Backtrace:\n" << boost::stacktrace::stacktrace() << '\n';
  std::abort();
}

}  // end namespace lf::base
