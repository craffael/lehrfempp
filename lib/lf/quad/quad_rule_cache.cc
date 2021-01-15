/**
 * @file
 * @brief Implementation of QuadRuleCache
 * @author Raffael Casagrande
 * @date   2021-01-13 06:27:46
 * @copyright MIT License
 */

#include "quad_rule_cache.h"
#include "make_quad_rule.h"

namespace lf::quad {
const QuadRule& QuadRuleCache::Get(base::RefEl ref_el, unsigned degree) const {
  LF_ASSERT_MSG(degree >= 0, "degree must be non-negative.");
  auto& vector = cache_[ref_el.Id()];
  if (vector.size() < degree + 1) {
    vector.resize(degree + 1);
  }
  auto& qr = vector[degree];
  if (qr.NumPoints() == 0) {
    qr = make_QuadRule(ref_el, degree);
  }
  LF_ASSERT_MSG(qr.Degree() >= degree, "something is horribly wrong.");
  return qr;
}
}  // namespace lf::quad
