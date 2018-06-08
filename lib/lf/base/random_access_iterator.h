#ifndef __2ed3a5566eea47a19b103fe7d0f69aca
#define __2ed3a5566eea47a19b103fe7d0f69aca
#include <list>
#include "forward_iterator.h"
#include "invalid_type_exception.h"

namespace lf::base {

template <class T>
class RandomAccessIterator : public ForwardIterator<T> {
  using base_t = ForwardIterator<T>;

 public:
  /**
   * @brief The type that is obtained when two pointers are subtracted from each
   * other.
   */
  using difference_type = std::ptrdiff_t;

 protected:
  using ForwardWrapperInterface = typename ForwardIterator<T>::WrapperInterface;

  class WrapperInterface : public virtual ForwardWrapperInterface {
   public:
    using ForwardWrapperInterface =
        typename ForwardIterator<T>::WrapperInterface;

    virtual T& operator[](difference_type i) const = 0;
    virtual WrapperInterface* operator+=(difference_type n) = 0;
    virtual WrapperInterface* operator-=(difference_type n) = 0;
    virtual std::unique_ptr<WrapperInterface> operator+(
        difference_type n) const = 0;
    virtual std::unique_ptr<WrapperInterface> operator-(
        difference_type n) const = 0;
    virtual difference_type operator-(
        const ForwardWrapperInterface* rhs) const = 0;
    virtual bool operator<(const ForwardWrapperInterface* other) const = 0;
    virtual bool operator<=(const ForwardWrapperInterface* other) const = 0;
    virtual bool operator>(const ForwardWrapperInterface* other) const = 0;
    virtual bool operator>=(const ForwardWrapperInterface* other) const = 0;
  };

  template <class InnerIterator,
            typename = typename std::enable_if<
                !std::is_reference<InnerIterator>::value>::type,
            typename = typename std::enable_if<std::is_same<
                typename std::iterator_traits<InnerIterator>::iterator_category,
                std::random_access_iterator_tag>::value>::type>
  class WrapperImpl
      : public virtual WrapperInterface,
        public ForwardIterator<T>::template WrapperImpl<InnerIterator> {
    using base_t =
        typename ForwardIterator<T>::template WrapperImpl<InnerIterator>;

   public:
    using ForwardWrapperInterface =
        typename ForwardIterator<T>::WrapperInterface;

    WrapperImpl(const InnerIterator& iterator) : base_t(iterator) {}

    WrapperImpl(InnerIterator&& iterator) : base_t(iterator) {}

    T& operator[](difference_type i) const override {
      return base_t::iterator_[i];
    }

    WrapperInterface* operator+=(difference_type n) override {
      base_t::iterator_ += n;
      return this;
    }

    WrapperInterface* operator-=(difference_type n) override {
      base_t::iterator_ -= n;
      return this;
    }

    std::unique_ptr<WrapperInterface> operator+(
        difference_type n) const override {
      return std::make_unique<WrapperImpl>(base_t::iterator_ + n);
    }

    std::unique_ptr<WrapperInterface> operator-(
        difference_type n) const override {
      return std::make_unique<WrapperImpl>(base_t::iterator_ - n);
    }

    difference_type operator-(
        const ForwardWrapperInterface* rhs) const override {
      const base_t* rhs_casted = dynamic_cast<const base_t*>(rhs);
      if (rhs) {
        return base_t::iterator_ - rhs_casted->iterator_;
      } else {
        throw InvalidTypeException(
            std::string("Cannot compare an random access iterator of type ") +
            typeid(*this).name() + " with iterator of type " +
            typeid(*rhs).name());
      }
    }

    bool operator<(const ForwardWrapperInterface* other) const override {
      const base_t* other_casted = dynamic_cast<const base_t*>(other);
      if (other_casted) {
        return base_t::iterator_ < other_casted->iterator_;
      } else {
        throw InvalidTypeException(
            std::string("Cannot compare an random access iterator of type ") +
            typeid(*this).name() + " with iterator of type " +
            typeid(*other).name());
      }
    }

    bool operator<=(const ForwardWrapperInterface* other) const override {
      const base_t* other_casted = dynamic_cast<const base_t*>(other);
      if (other_casted) {
        return base_t::iterator_ <= other_casted->iterator_;
      } else {
        throw InvalidTypeException(
            std::string("Cannot compare an random access iterator of type ") +
            typeid(*this).name() + " with iterator of type " +
            typeid(*other).name());
      }
    }

    bool operator>(const ForwardWrapperInterface* other) const override {
      const base_t* other_casted = dynamic_cast<const base_t*>(other);
      if (other_casted) {
        return base_t::iterator_ > other_casted->iterator_;
      } else {
        throw InvalidTypeException(
            std::string("Cannot compare an random access iterator of type ") +
            typeid(*this).name() + " with iterator of type " +
            typeid(*other).name());
      }
    }

    bool operator>=(const ForwardWrapperInterface* other) const override {
      const base_t* other_casted = dynamic_cast<const base_t*>(other);
      if (other_casted) {
        return base_t::iterator_ >= other_casted->iterator_;
      } else {
        throw InvalidTypeException(
            std::string("Cannot compare an random access iterator of type ") +
            typeid(*this).name() + " with iterator of type " +
            typeid(*other).name());
      }
    }

    virtual std::unique_ptr<typename ForwardIterator<T>::WrapperInterface>
    Clone() const override {
      return std::make_unique<WrapperImpl>(base_t::iterator_);
    }
  };

  class WrapperNull : public base_t::WrapperNull, WrapperInterface {
   public:
    using base_t = typename ForwardIterator<T>::WrapperNull;
    using ForwardWrapperInterface = typename base_t::ForwardWrapperInterface;

    T& operator[](difference_type i) const override {
      LF_VERIFY_MSG(
          false,
          "Cannot dereference an iterator that has been default constructed.");
    }

    WrapperInterface* operator+=(difference_type n) override {
      LF_VERIFY_MSG(
          false,
          "Cannot add to an iterator that has been default constructed.");
    }

    WrapperInterface* operator-=(difference_type n) override {
      LF_VERIFY_MSG(false,
                    "Cannot subtract from an iterator that has been default "
                    "constructed.");
    }

    std::unique_ptr<WrapperInterface> operator+(
        difference_type n) const override {
      LF_VERIFY_MSG(false, "Cannot add to a default-constructed iterator.");
    }

    std::unique_ptr<WrapperInterface> operator-(
        difference_type n) const override {
      LF_VERIFY_MSG(false,
                    "Cannot subtract from a default-constructed iterator.");
    }

    difference_type operator-(
        const ForwardWrapperInterface* rhs) const override {
      const base_t* rhs_casted = dynamic_cast<const base_t*>(rhs);
      if (rhs_casted) {
        return 0;
      } else {
        LF_VERIFY_MSG(false,
                      "cannot subtract an non-default-constructed iterator "
                      "from a default-constructed iterator.");
      }
    }

    bool operator<(const ForwardWrapperInterface* other) const override {
      const base_t* other_casted = dynamic_cast<const base_t*>(other);
      if (other_casted) {
        return false;
      } else {
        LF_VERIFY_MSG(false,
                      "operator< not defined for a default-constructed "
                      "iterator and a a non-default-constructed iterator.");
      }
    }

    bool operator<=(const ForwardWrapperInterface* other) const override {
      const base_t* other_casted = dynamic_cast<const base_t*>(other);
      if (other_casted) {
        return true;
      } else {
        LF_VERIFY_MSG(false,
                      "operator<= not defined for a default-constructed "
                      "iterator and a a non-default-constructed iterator.");
      }
    }

    bool operator>(const ForwardWrapperInterface* other) const override {
      const base_t* other_casted = dynamic_cast<const base_t*>(other);
      if (other_casted) {
        return false;
      } else {
        LF_VERIFY_MSG(false,
                      "operator> not defined for a default-constructed "
                      "iterator and a a non-default-constructed iterator.");
      }
    }

    bool operator>=(const ForwardWrapperInterface* other) const override {
      const base_t* other_casted = dynamic_cast<const base_t*>(other);
      if (other_casted) {
        return true;
      } else {
        LF_VERIFY_MSG(false,
                      "operator>= not defined for a default-constructed "
                      "iterator and a a non-default-constructed iterator.");
      }
    }

    virtual std::unique_ptr<typename ForwardIterator<T>::WrapperInterface>
    Clone() const override {
      return std::make_unique<WrapperNull>(base_t::iterator_);
    }
  };

  RandomAccessIterator(std::unique_ptr<WrapperInterface>&& ptr)
      : base_t(std::move(ptr)) {}

 public:
  template <class IteratorImpl,
            typename = typename std::enable_if<
                !std::is_reference<IteratorImpl>::value>::type,
            typename = typename std::enable_if<std::is_same<
                typename std::iterator_traits<IteratorImpl>::iterator_category,
                std::random_access_iterator_tag>::value>::type>
  RandomAccessIterator(const IteratorImpl& iterator)
      : base_t(std::make_unique<WrapperImpl<IteratorImpl>>(iterator)) {}

  template <
      class IteratorImpl,
      typename = typename std::iterator_traits<IteratorImpl>::difference_type,
      typename = std::enable_if_t<
          !std::is_base_of<RandomAccessIterator, IteratorImpl>::value &&
          !std::is_reference<IteratorImpl>::value>>
  RandomAccessIterator(IteratorImpl&& iterator)
      : base_t(std::unique_ptr<ForwardWrapperInterface>(
            new WrapperImpl<IteratorImpl>(iterator))) {}

  RandomAccessIterator() : base_t(std::make_unique<WrapperNull>()) {}

  RandomAccessIterator(const RandomAccessIterator& other)
      : base_t(other.wrapper_->Clone()) {}

  RandomAccessIterator(RandomAccessIterator&&) = default;

  RandomAccessIterator& operator=(const RandomAccessIterator& rhs) {
    base_t::wrapper_ = rhs.wrapper_->Clone();
    return *this;
  }

  RandomAccessIterator& operator+=(difference_type n) {
    dynamic_cast<WrapperInterface&>(*base_t::wrapper_) += n;
    return *this;
  }

  RandomAccessIterator& operator-=(difference_type n) {
    dynamic_cast<WrapperInterface&>(*base_t::wrapper_) -= n;
    return *this;
  }

  RandomAccessIterator operator+(difference_type n) const {
    return dynamic_cast<const WrapperInterface&>(*base_t::wrapper_) + n;
  }

  RandomAccessIterator operator-(difference_type n) const {
    return dynamic_cast<const WrapperInterface&>(*base_t::wrapper_) - n;
  }

  difference_type operator-(const RandomAccessIterator& rhs) const {
    return dynamic_cast<const WrapperInterface&>(*base_t::wrapper_) -
           &dynamic_cast<const WrapperInterface&>(*rhs.wrapper_);
  }

  T& operator[](difference_type i) const {
    return dynamic_cast<const WrapperInterface&>(*base_t::wrapper_)[i];
  }

  bool operator<(const RandomAccessIterator& rhs) const {
    return dynamic_cast<const WrapperInterface&>(*base_t::wrapper_) <
           &dynamic_cast<const WrapperInterface&>(*rhs.wrapper_);
  }

  bool operator<=(const RandomAccessIterator& rhs) const {
    return dynamic_cast<const WrapperInterface&>(*base_t::wrapper_) <=
           &dynamic_cast<const WrapperInterface&>(*rhs.wrapper_);
  }

  bool operator>(const RandomAccessIterator& rhs) const {
    return dynamic_cast<const WrapperInterface&>(*base_t::wrapper_) >
           &dynamic_cast<const WrapperInterface&>(*rhs.wrapper_);
  }

  bool operator>=(const RandomAccessIterator& rhs) const {
    return dynamic_cast<const WrapperInterface&>(*base_t::wrapper_) >=
           &dynamic_cast<const WrapperInterface&>(*rhs.wrapper_);
  }
};

}  // namespace lf::base

/// \cond
namespace std {
template <class T>
struct iterator_traits<lf::base::RandomAccessIterator<T>> {
  using difference_type = std::ptrdiff_t;
  using value_type = T;
  using pointer = T*;
  using reference = T&;
  using iterator_category = std::random_access_iterator_tag;
};
}  // namespace std

/// \endcond

#endif  // __2ed3a5566eea47a19b103fe7d0f69aca
