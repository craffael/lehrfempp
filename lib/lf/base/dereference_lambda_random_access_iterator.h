#ifndef __ae2586d329d6457aa2449e5e097a2ed4
#define __ae2586d329d6457aa2449e5e097a2ed4
#include <functional>
#include <type_traits>
#include "random_access_iterator.h"

namespace lf::base {

template <class Iterator, class T,
          typename = typename std::enable_if<std::is_same<
              typename std::iterator_traits<Iterator>::iterator_category,
              std::random_access_iterator_tag>::value>::type>
class DereferenceLambdaRandomAccessIterator {
  using reference_type = const T&;
  using lambda_t = std::function<reference_type(const Iterator&)>;

 public:
  DereferenceLambdaRandomAccessIterator(Iterator&& iterator, lambda_t lambda)
      : iterator_(iterator), lambda_(std::move(lambda)) {}

  DereferenceLambdaRandomAccessIterator(
      const DereferenceLambdaRandomAccessIterator&) = default;
  DereferenceLambdaRandomAccessIterator(
      DereferenceLambdaRandomAccessIterator&&) noexcept = default;
  DereferenceLambdaRandomAccessIterator& operator=(
      DereferenceLambdaRandomAccessIterator&&) noexcept = default;
  DereferenceLambdaRandomAccessIterator& operator=(
      const DereferenceLambdaRandomAccessIterator&) = default;
  ~DereferenceLambdaRandomAccessIterator() = default;

  auto& operator*() const { return lambda_(iterator_); }

  auto& operator[](std::ptrdiff_t i) const { return lambda_(iterator_ + i); }

  auto operator++() {
    ++iterator_;
    return *this;
  }

  auto operator++(int) {
    return DereferenceLambdaRandomAccessIterator(iterator_++, lambda_);
  }

  auto operator+=(std::ptrdiff_t n) {
    iterator_ += n;
    return *this;
  }

  auto operator--() {
    --iterator_;
    return *this;
  }

  auto operator--(int) {
    return DereferenceLambdaRandomAccessIterator(iterator_--, lambda_);
  }

  auto operator-=(std::ptrdiff_t n) {
    iterator_ -= n;
    return *this;
  }

  auto operator+(std::ptrdiff_t n) const {
    return DereferenceLambdaRandomAccessIterator(iterator_ + n, lambda_);
  }

  auto operator-(std::ptrdiff_t n) const {
    return DereferenceLambdaRandomAccessIterator(iterator_ - n, lambda_);
  }

  auto operator-(const DereferenceLambdaRandomAccessIterator rhs) const {
    return iterator_ - rhs.iterator_;
  }

  auto operator==(const DereferenceLambdaRandomAccessIterator& rhs) const {
    return iterator_ == rhs.iterator_;
  }

  auto operator!=(const DereferenceLambdaRandomAccessIterator& rhs) const {
    return iterator_ != rhs.iterator_;
  }

  auto operator<(const DereferenceLambdaRandomAccessIterator& rhs) const {
    return iterator_ < rhs.iterator_;
  }

  auto operator<=(const DereferenceLambdaRandomAccessIterator& rhs) const {
    return iterator_ <= rhs.iterator_;
  }

  auto operator>(const DereferenceLambdaRandomAccessIterator& rhs) const {
    return iterator_ > rhs.iterator_;
  }

  auto operator>=(const DereferenceLambdaRandomAccessIterator& rhs) const {
    return iterator_ >= rhs.iterator_;
  }

 private:
  Iterator iterator_;
  lambda_t lambda_;
};

template <class Iterator, class Lambda>
DereferenceLambdaRandomAccessIterator<
    Iterator, typename std::remove_reference<decltype(
                  std::declval<Lambda>()(std::declval<Iterator>()))>::type>
make_DereferenceLambdaRandomAccessIterator(Iterator&& i, Lambda l) {
  return {std::forward<Iterator>(i), l};
}  // namespace lf::base

}  // namespace lf::base

namespace std {

template <class Iterator, class T>
struct iterator_traits<
    lf::base::DereferenceLambdaRandomAccessIterator<Iterator, T>> {
  using difference_type = std::ptrdiff_t;
  using value_type = T;

  using pointer = value_type*;
  using reference = value_type&;
  using iterator_category = std::random_access_iterator_tag;
};

}  // namespace std

#endif  // __ae2586d329d6457aa2449e5e097a2ed4
