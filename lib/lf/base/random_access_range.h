#ifndef __c23d69f309b14fbbbf43f5f0ed9a97fa
#define __c23d69f309b14fbbbf43f5f0ed9a97fa
#include "random_access_iterator.h"
#include "forward_range.h"

namespace lf::base {

template <class T>
class RandomAccessRange {
  RandomAccessIterator<T> begin_;
  RandomAccessIterator<T> end_;

public:
  RandomAccessRange(const RandomAccessRange&) = default;
  RandomAccessRange(RandomAccessRange&&) = default;
  RandomAccessRange& operator=(const RandomAccessRange&) = default;
  RandomAccessRange& operator=(RandomAccessRange&&) = default;
  ~RandomAccessRange() = default;

  RandomAccessRange(RandomAccessIterator<T> begin, RandomAccessIterator<T> end)
    : begin_(begin),
      end_(end) {
  }

  RandomAccessIterator<T> begin() const { return begin_; }
  RandomAccessIterator<T> end() const { return end_; }

  operator ForwardRange<T>() const {
    return new ForwardRange<T>(begin_, end_);
  }
};

}


#endif // __c23d69f309b14fbbbf43f5f0ed9a97fa
