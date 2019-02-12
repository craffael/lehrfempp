/** @file static_vars.cc */

#include "static_vars.h"
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

namespace lf::base {

// Handling of global control variables
// Global variable !
StaticVar *ctrl_root = nullptr;

}  // end namespace lf::base
