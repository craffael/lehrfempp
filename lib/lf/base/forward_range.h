#ifndef __661aa07c70a04907a44d828a353e6537
#define __661aa07c70a04907a44d828a353e6537

#include "forward_iterator.h"

namespace lf::base {

/**
 * @brief A pair of ForwardIterator
 * @tparam T 
 */
template <class T>
class Range {
  ForwardIterator<T> begin_;
  ForwardIterator<T> end_;

public:

  Range(const Range& ) = default;
  Range(Range&&) = default;
  Range& operator=(const Range&) = default;
  Range& operator=(Range&&) = default;

  Range(ForwardIterator<T> begin, ForwardIterator<T> end)
    : begin_(begin),
      end_(end) {
  }

  ForwardIterator<T> begin() const { return begin_; }
  ForwardIterator<T> end() const {return end_; }
};


}


#endif // __661aa07c70a04907a44d828a353e6537
