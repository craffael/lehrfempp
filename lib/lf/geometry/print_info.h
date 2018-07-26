#ifndef PRINT_INFO_H
#define PRINT_INFO_H

#include "lf/geometry/geometry.h"


namespace lf::geometry {


/**
 * @brief Diagnostic output operator for Geometry
 * @param &geom The geometry to print info about
 * @param &o The stream to which this function should output
 *
 */
void PrintInfo(const Geometry &geom, std::ostream &o);


} // namespace lf::geometry

#endif // PRINT_INFO_H
