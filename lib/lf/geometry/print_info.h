#ifndef PRINT_INFO_H
#define PRINT_INFO_H

#include "lf/geometry/geometry.h"

namespace lf::geometry {

/**
 * @brief Diagnostic output operator for Geometry
 * @param geom The geometry to print info about
 * @param o The stream to which this function should output
 * @param output_ctrl Controls the level of detail of the written output (see
 * below)
 *
 * #### Output levels
 * - output_ctrl == 0: Reference element type of geometry is printed
 * - output_ctrl > 10: The above and global and local dimension and
 * derived type of reference element is printed.
 * - output_ctrl > 90: The above and coordinates og the geometry is
 * printed.
 *
 * @relates Geometry
 */
void PrintInfo(std::ostream& o, const Geometry& geom, int output_ctrl = 0);

/**
 * @brief Operator overload to print a `Geometry` to a stream, such as
 * `std::cout`
 * @param stream The stream to which this function should output
 * @param geom The geometry to write to `stream`.
 * @return The stream itself.
 *
 * @note At the moment this will just write the type of the reference element to
 * the stream.
 *
 * @relates Geometry
 */
std::ostream& operator<<(std::ostream& stream, const Geometry& geom);

}  // namespace lf::geometry

/**
 * @brief Make lf::geometry::Geometry formattable by fmt
 * (https://fmt.dev/latest/api.html#ostream-api)
 */
template <>
struct fmt::formatter<lf::geometry::Geometry> : ostream_formatter {};

#endif  // PRINT_INFO_H
