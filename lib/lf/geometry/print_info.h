#ifndef PRINT_INFO_H
#define PRINT_INFO_H

// Include something?
#include "lf/geometry/geometry.h"


namespace lf::geometry {


/**
 * @brief diagnostic output operator
 */

// Add print function here
// Only the declaration
// Need different functions for each type of geometry? Or for single geometry.
// A: Think for single geometry only
void PrintInfo(const Geometry &geom, std::ostream &o);




} // namespace lf::geometry

#endif // PRINT_INFO_H
