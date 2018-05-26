#ifndef __ae2586d329d6457aa2449e5e097a2ed4
#define __ae2586d329d6457aa2449e5e097a2ed4
#include <type_traits>
#include "random_access_iterator.h"

namespace lf::base {

template <class Iterator, typename = typename std::enable_if<std::is_same<
            typename std::iterator_traits<Iterator>::iterator_category, 
            std::random_access_iterator_tag>::value>::type>
class DereferenceRandomAccessIterator {

public:

  DereferenceRandomAccessIterator(Iterator&& iterator)
    : iterator_(iterator) {
  }

  DereferenceRandomAccessIterator() = default;
  DereferenceRandomAccessIterator(const DereferenceRandomAccessIterator&) = default;
  DereferenceRandomAccessIterator& operator=(const DereferenceRandomAccessIterator& ) = default;
  DereferenceRandomAccessIterator(DereferenceRandomAccessIterator&&) = default;
  DereferenceRandomAccessIterator& operator=(DereferenceRandomAccessIterator&&) = default;
  ~DereferenceRandomAccessIterator() = default;


  auto& operator*() const {
    return **iterator_;
  }

  auto& operator[](std::ptrdiff_t i) const {
    return *iterator_[i];
  }

  auto operator++() {
    ++iterator_;
    return *this;
  }

  auto operator++(int) {
    return DereferenceRandomAccessIterator(iterator_++);
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
    return DereferenceRandomAccessIterator(iterator_--);
  }

  auto operator-=(std::ptrdiff_t n) {
    iterator_ -= n;
    return *this;
  }

  auto operator+(std::ptrdiff_t n) const {
    return DereferenceRandomAccessIterator(iterator_ + n);
  }

  auto operator-(std::ptrdiff_t n) const {
    return DereferenceRandomAccessIterator(iterator_ - n);
  }

  auto operator-(const DereferenceRandomAccessIterator rhs) const {
    return iterator_ - rhs.iterator_;
  }

  auto operator==(const DereferenceRandomAccessIterator& rhs) const {
    return iterator_ == rhs.iterator_;
  }

  auto operator!=(const DereferenceRandomAccessIterator& rhs) const {
    return iterator_ != rhs.iterator_;
  }

  auto operator<(const DereferenceRandomAccessIterator& rhs) const {
    return iterator_ < rhs.iterator_;
  }

  auto operator<=(const DereferenceRandomAccessIterator& rhs) const {
    return iterator_ <= rhs.iterator_;
  }

  auto operator>(const DereferenceRandomAccessIterator& rhs) const {
    return iterator_ > rhs.iterator_;
  }

  auto operator>=(const DereferenceRandomAccessIterator& rhs) const {
    return iterator_ >= rhs.iterator_;
  }

private:
  Iterator iterator_;

};

}

namespace std {

template<class Iterator>
struct iterator_traits<lf::base::DereferenceRandomAccessIterator<Iterator>> {
  using difference_type = std::ptrdiff_t;
  using value_type = std::remove_reference<decltype(**std::declval<Iterator>())>;
  using pointer = value_type*;
  using reference = value_type&;
  using iterator_category = std::random_access_iterator_tag;
};

}


#endif // __ae2586d329d6457aa2449e5e097a2ed4
