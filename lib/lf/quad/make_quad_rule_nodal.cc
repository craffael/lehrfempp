/**
 * @file
 * @brief Implementation of the make_QuadRuleNodal.
 * @author Raffael Casagrande
 * @date   2018-12-26 10:25:52
 * @copyright MIT License
 */

#include "make_quad_rule_nodal.h"
#include <boost/core/ref.hpp>

namespace lf::quad {
QuadRule make_QuadRuleNodal(base::RefEl ref_el) {
  double ref_el_vol = 1.0;
  if (ref_el == base::RefEl::kTria()) {
    ref_el_vol = 0.5;
  }
  return QuadRule(ref_el, ref_el.NodeCoords(),
                  Eigen::VectorXd::Constant(ref_el.NumNodes(),
                                            ref_el_vol / ref_el.NumNodes()),
                  1);
}
}  // namespace lf::quad
