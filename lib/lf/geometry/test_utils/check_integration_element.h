/**
 * @file
 * @brief Declaration of IntegrationElement() test for geometry objects
 * @author Anian Ruoss
 * @date   2019-02-30 17:59:17
 * @copyright MIT License
 */

#ifndef LEHRFEMPP_CHECK_INTEGRATION_ELEMENT_H
#define LEHRFEMPP_CHECK_INTEGRATION_ELEMENT_H

#include <lf/geometry/geometry.h>

namespace lf::geometry::test_utils {

/**
 * @brief Checks if IntegrationElement() is implemented correctly under the
 * assumption that Jacobian() is correct
 * @param geom geometry object to be evaluated
 * @param eval_points points at which Jacobian should be checked
 */
void checkIntegrationElement(const lf::geometry::Geometry &geom,
                             const Eigen::MatrixXd &eval_points);

}  // namespace lf::geometry::test_utils

#endif  // LEHRFEMPP_CHECK_INTEGRATION_ELEMENT_H
