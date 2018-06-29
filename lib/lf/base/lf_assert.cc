/** @file lf_assert.cc */

// Preprocessor flag to prevent double definition of global handler object
// for control variables
#define CTRLVARROOT

#include "lf_assert.h"

namespace lf::base {
  
  // Output for assertions
  void AssertionFailed(const std::string& expr, const std::string& file,
		       int line, const std::string& msg) {
    std::cerr << "***** Internal Program Error - assertion (" << expr
	      << ") failed:\n"
	      << file << '(' << line << "): " << msg << std::endl;
  }

} // end namespace lf::base
