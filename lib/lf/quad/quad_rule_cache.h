/**
 * @file
 * @brief A class which returns a quadrature rule for the given reference
 * element and order that supports caching.
 * @author Raffael Casagrande
 * @date   2021-01-13 06:15:05
 * @copyright MIT License
 */

#ifndef __5b965019a2a74585a86f055191a31def
#define __5b965019a2a74585a86f055191a31def

#include <deque>
#include "quad_rule.h"

namespace lf::quad {

/**
 * @brief A cache for make_QuadRule()
 *
 * This class has one important method: QuadRuleCache::Get()
 * which returns the same result as make_QuadRule() but uses
 * caching to construct the QuadRule only the first time.
 */
class QuadRuleCache {
 public:
  QuadRuleCache() = default;

  /**
   * @brief Copy constructor is deleted to avoid accidental copy.
   */
  QuadRuleCache(const QuadRuleCache&) = delete;

  /**
   * @brief Move construction is allowed.
   */
  QuadRuleCache(QuadRuleCache&& other) noexcept
      : cache_(std::move(other.cache_)) {}

  /**
   * @brief copy assignment is delete to avoid accidental copy.
   */
  QuadRuleCache& operator=(const QuadRuleCache&) = delete;

  /**
   * @brief Move assignment is allowed.
   */
  QuadRuleCache& operator=(QuadRuleCache&&) noexcept = default;

  ~QuadRuleCache() = default;

  /**
   * @brief Retrieve a quadrature rule for reference element `ref_el` with
   * degree `d`. Is identical to calling `make_QuadRule()`.
   * @param ref_el The reference element for which the quadrature rule is.
   * @param degree The degree of the quadrature rule
   * @return A quadrature rule for reference element `ref_el` with degree `d`.
   *
   * @see make_QuadRule()
   *
   * @note It is guaranteed that the returned reference will remain valid for
   * the entire lifetime of the QuadRuleCache.
   */
  [[nodiscard]] const QuadRule& Get(base::RefEl ref_el, unsigned degree) const;

 private:
  // cache_[i][j] contains quadrature rule for reference element with id = i and
  // degree = j
  mutable std::array<std::deque<QuadRule>, 5> cache_;
};

}  // namespace lf::quad

#endif  // __5b965019a2a74585a86f055191a31def
