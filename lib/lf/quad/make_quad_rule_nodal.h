/**
 * @file
 * @brief A function to create nodal quadrature rules, i.e. quadrature rules
 * that only evaluate the function on nodes of the reference element
 * @author Raffael Casagrande
 * @date   2018-12-26 10:21:58
 * @copyright MIT License
 */

#ifndef INCG6c48015148ae47f384263808689bdaf3
#define INCG6c48015148ae47f384263808689bdaf3
#include "quad_rule.h"

namespace lf::quad {

/**
 * @brief Create a quadrature rule that evaluates the quadrand only at the nodes
 *        of the reference element
 * @param ref_el The reference element for which the quadrule is.
 * @return A quadrature rule of order 1.
 */
QuadRule make_QuadRuleNodal(base::RefEl ref_el);

}  // namespace lf::quad

#endif  // INCG6c48015148ae47f384263808689bdaf3
