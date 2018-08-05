/**
 * @file
 * @brief Defines a predicate that always returns true.
 * @author Raffael Casagrande
 * @date   2018-08-04 03:26:32
 * @copyright MIT License
 */

#ifndef __7f5cb415f32d480ba60fcefa5f351637
#define __7f5cb415f32d480ba60fcefa5f351637

namespace lf::base {

/**
 * @brief A <a href="https://en.wikipedia.org/wiki/Function_object">Function
 * Object</a> that can be invoked with any arguments and that always returns the
 *        value true.
 */
class PredicateTrue {
 public:
  template <class... T>
  bool operator()(T&&...) const {
    return true;
  }
};

}  // namespace lf::base

#endif  // __7f5cb415f32d480ba60fcefa5f351637
