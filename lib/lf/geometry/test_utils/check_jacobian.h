/**
 * @file
 * @brief Declaration of Jacobian() test for geometry objects
 * @author Anian Ruoss
 * @date   2019-02-30 17:30:17
 * @copyright MIT License
 */

#ifndef LEHRFEMPP_CHECK_JACOBIAN_H
#define LEHRFEMPP_CHECK_JACOBIAN_H

#include <lf/geometry/geometry.h>

namespace lf::geometry::test_utils {

/**
 * @brief Checks if Jacobian() is implemented correctly by comparing it to
 * the symmetric difference quotient approximation
 * @param geom geometry object to be evaluated
 * @param eval_points points at which Jacobian should be checked
 * @param tolerance tolerance for floating point equality check
 */
void checkJacobian(const lf::geometry::Geometry &geom,
                   const Eigen::MatrixXd &eval_points, const double &tolerance);

}  // namespace lf::geometry::test_utils

#endif  // LEHRFEMPP_CHECK_JACOBIAN_H
