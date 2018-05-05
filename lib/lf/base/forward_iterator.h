#ifndef __20db9042f1c04f5c8c173f3e354a915c
#define __20db9042f1c04f5c8c173f3e354a915c
#include <memory>
#include "lf_assert.h"


namespace lf::base {

template <class T>
class ForwardIterator {

  class WrapperInterface {
  public:
    virtual std::unique_ptr<WrapperInterface> Clone() const = 0;
    virtual bool Compare(const WrapperInterface* other) const = 0;
    virtual T& Dereference() const = 0;
    virtual WrapperInterface* operator++() = 0;
    virtual std::unique_ptr<WrapperInterface> operator++(int) = 0;


    virtual ~WrapperInterface() {
    }
  };

  template <class InnerIterator>
  class WrapperImpl : public WrapperInterface {
  private:
    InnerIterator iterator_;

  public:

    WrapperImpl(const InnerIterator& iterator)
      : iterator_(iterator) {
    }

    std::unique_ptr<WrapperInterface> Clone() const override {
      return std::make_unique<WrapperImpl>(iterator_);
    }

    bool Compare(const WrapperInterface* other) const override {
      const WrapperImpl* other_casted = dynamic_cast<const WrapperImpl*>(other);
      if (other_casted) {
        return other_casted->iterator_ == iterator_;
      } else {
        return false;
      }
    }

    T& Dereference() const override {
      return *iterator_;
    }

    WrapperInterface* operator++() override {
      ++iterator_;
      return this;
    }

    std::unique_ptr<WrapperInterface> operator++(int) override {
      return std::make_unique<WrapperImpl>(iterator_++);
    }

  };

  class WrapperNull : public WrapperInterface {

  public:
    std::unique_ptr<WrapperInterface> Clone() const override {
      return std::make_unique<WrapperNull>();
    }

    bool Compare(const WrapperInterface* other) const override {
      return dynamic_cast<const WrapperNull*>(other);
    }

    T& Dereference() const override {
      LF_VERIFY_MSG(false,
        "Cannot dereference a ForwardIterator that has been default constructed."
      );
    }

    WrapperInterface* operator++() override {
      LF_VERIFY_MSG(false,
        "Cannot increment a ForwardIterator that has been default constructed."
      );
    }

    std::unique_ptr<WrapperInterface> operator++(int) override {
      LF_VERIFY_MSG(false,
        "Cannot increment a ForwardIterator that has been default constructed."
      );
    }
  };


  std::unique_ptr<WrapperInterface> wrapper_;
  // needed by operator++()
  ForwardIterator(std::unique_ptr<WrapperInterface>&& ptr)
    : wrapper_(std::move(ptr)) {
  }

public:

  template <class IteratorImpl>
  ForwardIterator(const IteratorImpl& iterator)
    : wrapper_(std::make_unique<WrapperImpl<IteratorImpl>>(iterator)) {
  }

  template <class IteratorImpl, typename = std::enable_if_t<!std::
              is_convertible_v<IteratorImpl, ForwardIterator&>>>
  ForwardIterator(IteratorImpl&& iterator)
    : wrapper_(std::make_unique<WrapperImpl<IteratorImpl>>(iterator)) {
  }

  ForwardIterator()
    : wrapper_(std::make_unique<WrapperNull>()) {
  }

  ForwardIterator(const ForwardIterator& other)
    : wrapper_(other.wrapper_->Clone()) {
  }

  ForwardIterator(ForwardIterator&&) = default;
  ForwardIterator& operator=(ForwardIterator&&) = default;

  ForwardIterator& operator=(const ForwardIterator& rhs) {
    wrapper_ = rhs.wrapper_->Clone();
    return *this;
  }

  bool operator==(const ForwardIterator& rhs) const {
    return wrapper_->Compare(rhs.wrapper_.get());
  }

  bool operator!=(const ForwardIterator& rhs) const {
    return !operator==(rhs);
  }

  T& operator*() const {
    return wrapper_->Dereference();
  }

  T* operator->() const {
    return &wrapper_->Dereference();
  }

  ForwardIterator& operator++() {
    wrapper_->operator++();
    return *this;
  }

  ForwardIterator operator++(int) {
    return ForwardIterator(wrapper_->operator++(0));
  }

};

}

namespace std {
  template<class T>
  struct iterator_traits<lf::base::ForwardIterator<T>> {
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = T*;
    using reference = T&;
    using iterator_category = std::forward_iterator_tag;
  };
}


#endif // __20db9042f1c04f5c8c173f3e354a915c
