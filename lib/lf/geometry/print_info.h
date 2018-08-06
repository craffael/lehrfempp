#ifndef PRINT_INFO_H
#define PRINT_INFO_H

#include "lf/geometry/geometry.h"


namespace lf::geometry {


/**
 * @brief Diagnostic output operator for Geometry
 * @param &geom The geometry to print info about
 * @param &o The stream to which this function should output
 *
 * #### Output levels
 * - Geometry::output_ctrl_ == 0: Reference element type of geometry is printed
 * - Geometry::output_ctrl_ > 10: The above and global and local dimension and derived type of reference element is printed.
 * - Geometry::output_ctrl_ > 90: The above and coordinates og the geometry is printed.
 *
 */
void PrintInfo(const Geometry &geom, std::ostream &o);


/**
 * @brief Operator overload to print a `Geometry` to a stream, such as `std::cout`
 * @param stream The stream to which this function should output
 * @param entity The geometry to write to `stream`.
 * @return The stream itself.
 *
 * - If Geometry::output_ctrl_ == 0, type reference element of geometry is sent as output to stream
 * - If Geometry::output_ctrl_ > 0, then lf::geometry::PrintInfo(const Geometry &geom, std::ostream &o) is called.
 *
 */
 std::ostream& operator<<(std::ostream& stream, const Geometry& geom);

} // namespace lf::geometry

#endif // PRINT_INFO_H
