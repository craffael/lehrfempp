#ifndef __661aa07c70a04907a44d828a353e6537
#define __661aa07c70a04907a44d828a353e6537

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
  ForwardIterator<T> begin_;
  ForwardIterator<T> end_;

public:

  ForwardRange(const ForwardRange& ) = default;
  ForwardRange(ForwardRange&&) = default;
  ForwardRange& operator=(const ForwardRange&) = default;
  ForwardRange& operator=(ForwardRange&&) = default;

  ForwardRange(ForwardIterator<T> begin, ForwardIterator<T> end)
    : begin_(begin),
      end_(end) {
  }

  ForwardIterator<T> begin() const { return begin_; }
  ForwardIterator<T> end() const {return end_; }
};


}


#endif // __661aa07c70a04907a44d828a353e6537
