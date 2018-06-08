#ifndef __661aa07c70a04907a44d828a353e6537
#define __661aa07c70a04907a44d828a353e6537

#include <iterator>
#include "forward_iterator.h"

namespace lf::base {

/**
 * @brief A pair of ForwardIterator's that make up a set of elements of type T
 * @tparam T The type of elements contained in the range.
 *
 * ### Motivation
 */
template <class T>
class ForwardRange {
 protected:
  ForwardIterator<T> begin_;
  ForwardIterator<T> end_;

 public:
  ForwardRange(const ForwardRange&) = default;
  ForwardRange(ForwardRange&&) = default;

  template <class C, typename = typename std::enable_if<std::is_same<
                         typename std::iterator_traits<decltype(
                             std::declval<C>().begin())>::iterator_category,
                         std::forward_iterator_tag>::value>::type>
  ForwardRange(const C& forward_range)
      : begin_(forward_range.begin()), end_(forward_range.end()) {}

  ForwardRange(std::initializer_list<T> initializer_list)
      : begin_(initializer_list.begin()), end_(initializer_list.end()) {}

  ForwardRange& operator=(const ForwardRange&) = default;
  ForwardRange& operator=(ForwardRange&&) = default;

  ~ForwardRange() = default;

  ForwardRange(ForwardIterator<T> begin, ForwardIterator<T> end)
      : begin_(begin), end_(end) {}

  ForwardIterator<T> begin() const { return begin_; }
  ForwardIterator<T> end() const { return end_; }
};

}  // namespace lf::base

#endif  // __661aa07c70a04907a44d828a353e6537
