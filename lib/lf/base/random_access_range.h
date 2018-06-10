#ifndef __c23d69f309b14fbbbf43f5f0ed9a97fa
#define __c23d69f309b14fbbbf43f5f0ed9a97fa
#include <utility>
#include "forward_range.h"
#include "random_access_iterator.h"

namespace lf::base {

template <class T>
class RandomAccessRange {
  RandomAccessIterator<T> begin_;
  RandomAccessIterator<T> end_;

 public:
  RandomAccessRange(const RandomAccessRange&) = default;
  RandomAccessRange(RandomAccessRange&&) noexcept = default;
  RandomAccessRange& operator=(const RandomAccessRange&) = default;
  RandomAccessRange& operator=(RandomAccessRange&&) noexcept = default;
  ~RandomAccessRange() = default;

  RandomAccessRange(RandomAccessIterator<T> begin, RandomAccessIterator<T> end)
      : begin_(std::move(begin)), end_(std::move(end)) {}

  RandomAccessIterator<T> begin() const { return begin_; }
  RandomAccessIterator<T> end() const { return end_; }

  T& operator[](size_t i) const { return begin_[i]; }

  // NOLINTNEXTLINE(google-explicit-constructor, hicpp-explicit-conversions)
  operator ForwardRange<T>() const { return new ForwardRange<T>(begin_, end_); }
};

}  // namespace lf::base

#endif  // __c23d69f309b14fbbbf43f5f0ed9a97fa
