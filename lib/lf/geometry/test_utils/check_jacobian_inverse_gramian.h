/**
 * @file
 * @brief Declaration of JacobianInverseGramian() test for geometry objects
 * @author Anian Ruoss
 * @date   2019-02-30 17:53:17
 * @copyright MIT License
 */

#ifndef LEHRFEMPP_CHECK_JACOBIAN_INVERSE_GRAMIAN_H
#define LEHRFEMPP_CHECK_JACOBIAN_INVERSE_GRAMIAN_H

#include <lf/geometry/geometry.h>

namespace lf::geometry::test_utils {

/**
 * @brief Checks if JacobianInverseGramian() is implemented correctly
 * assuming that Jacobian() is correct
 * @param geom geometry object to be evaluated
 * @param eval_points points ar which JacobianInverseGramian should be checked
 * @param precision The precision with which the equivalence should be checked.
 */
void checkJacobianInverseGramian(const lf::geometry::Geometry &geom,
                                 const Eigen::MatrixXd &eval_points,
                                 double precision = 1e-12);

}  // namespace lf::geometry::test_utils

#endif  // LEHRFEMPP_CHECK_JACOBIAN_INVERSE_GRAMIAN_H
