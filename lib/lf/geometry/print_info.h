#ifndef PRINT_INFO_H
#define PRINT_INFO_H

#include "lf/geometry/geometry.h"


namespace lf::geometry {


/**
 * @brief diagnostic output operator for Geometry
 */
void PrintInfo(const Geometry &geom, std::ostream &o);




} // namespace lf::geometry

#endif // PRINT_INFO_H
