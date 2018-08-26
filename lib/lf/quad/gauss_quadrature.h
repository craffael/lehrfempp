/**
 * @file
 * @brief Provides functions for different gauss-type quadrature rules.
 * @author Raffael Casagrande
 * @date   2018-08-11 09:09:08
 * @copyright MIT License
 */

#ifndef __cd0bee55c3bb44e8bda00e45e61f2001
#define __cd0bee55c3bb44e8bda00e45e61f2001
#include <tuple>
#include "quad_rule.h"

namespace lf::quad {

/**
 * @brief Computes the quadrature points and weights for the interval [0,1]
 *  of a Gauss-Legendre quadrature rule using a newton algorithm.
 * @param num_points The number of points/weights that should be used.
 * @return The nodes (first) and weights(second) of the quadrature rule.
 */
std::tuple<Eigen::VectorXd, Eigen::VectorXd> GaussLegendre(
    unsigned int num_points);

/**
 * @brief Computes the quadrature points and weights for the interval [-1,1]
 *  of a Gauss-Jacobi quadrature rule with weight function \f$ (1-x)^\alpha
 * (1+x)^\beta \f$
 * @param num_points The number of points/weights that should be used.
 * @param alpha  The exponent of the weight function
 * @param beta   The exponent of the weight function
 * @return The nodes(first) and weights(second) of the quadrature rule
 */
std::tuple<Eigen::VectorXd, Eigen::VectorXd> GaussJacobi(quadOrder_t num_points,
                                                         double alpha,
                                                         double beta);

}  // namespace lf::quad

#endif  // __cd0bee55c3bb44e8bda00e45e61f2001
