#ifndef __ae2586d329d6457aa2449e5e097a2ed4
#define __ae2586d329d6457aa2449e5e097a2ed4
#include <type_traits>
#include "random_access_iterator.h"

namespace lf::base {

template <class Iterator, class Lambda, typename = typename std::enable_if<
            std::is_same<
              typename std::iterator_traits<Iterator>::iterator_category,
              std::random_access_iterator_tag>::value>::type>
class DereferenceLambdaRandomAccessIterator {

public:

  DereferenceLambdaRandomAccessIterator(Iterator&& iterator, Lambda lambda)
    : iterator_(iterator),
      lambda_(lambda) {
  }

  DereferenceLambdaRandomAccessIterator() = default;
  DereferenceLambdaRandomAccessIterator(
    const DereferenceLambdaRandomAccessIterator&) = default;
  DereferenceLambdaRandomAccessIterator& operator=(
    const DereferenceLambdaRandomAccessIterator&) = default;
  DereferenceLambdaRandomAccessIterator(DereferenceLambdaRandomAccessIterator&&)
  = default;
  DereferenceLambdaRandomAccessIterator& operator=(
    DereferenceLambdaRandomAccessIterator&&) = default;
  ~DereferenceLambdaRandomAccessIterator() = default;


  auto& operator*() const {
    return lambda_(iterator_);
  }

  auto& operator[](std::ptrdiff_t i) const {
    return lambda_(iterator_ + i);
  }

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
  Lambda lambda_;

};

template<class Iterator, class Lambda>
DereferenceLambdaRandomAccessIterator<Iterator, Lambda> make_DereferenceLambdaRandomAccessIterator(Iterator&& i, Lambda l) {
  return {std::move(i), l};
}

}

namespace std {

template <class Iterator, class Lambda>
struct iterator_traits<lf::base::DereferenceLambdaRandomAccessIterator<
    Iterator, Lambda>> {
  using difference_type = std::ptrdiff_t;
  using value_type = typename std::remove_reference<decltype(std::declval<Lambda>()(std::declval<Iterator>()))>::type;
  ;
  using pointer = value_type*;
  using reference = value_type&;
  using iterator_category = std::random_access_iterator_tag;
};

}


#endif // __ae2586d329d6457aa2449e5e097a2ed4
