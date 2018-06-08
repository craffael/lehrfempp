#include "lf_assert.h"
#include <iostream>

void lf::base::AssertionFailed(const std::string& expr, const std::string& file,
                               int line, const std::string& msg) {
  std::cerr << "***** Internal Program Error - assertion (" << expr
            << ") failed:\n"
            << file << '(' << line << "): " << msg << std::endl;
}
