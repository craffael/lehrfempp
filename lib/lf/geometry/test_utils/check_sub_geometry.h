/**
 * @file
 * @brief Declaration of SubGeometry() test for geometry objects
 * @author Anian Ruoss
 * @date   2019-02-30 18:06:17
 * @copyright MIT License
 */

#ifndef LEHRFEMPP_CHECK_SUB_GEOMETRY_H
#define LEHRFEMPP_CHECK_SUB_GEOMETRY_H

#include <lf/geometry/geometry.h>
#include <lf/quad/quad.h>

namespace lf::geometry::test_utils {

/**
 * @brief Checks that geometry objects obtained through SubGeometry() map the
 * same points on the reference elements to the same global points
 * @param geom geometry object to be evaluated
 * @param qrProvider function returning a QuadRule object associated with the
 * given RefEl
 */
void checkSubGeometry(
    const lf::geometry::Geometry &geom,
    const std::function<lf::quad::QuadRule(lf::base::RefEl)> &qrProvider);

}  // namespace lf::geometry::test_utils

#endif  // LEHRFEMPP_CHECK_SUB_GEOMETRY_H
