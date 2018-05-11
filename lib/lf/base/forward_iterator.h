#ifndef __20db9042f1c04f5c8c173f3e354a915c
#define __20db9042f1c04f5c8c173f3e354a915c
#include <memory>
#include "lf_assert.h"


namespace lf::base {

/**
 * @brief A wrapper around any <a href="http://en.cppreference.com/w/cpp/concept/ForwardIterator">Forward Iterator</a>
 * @tparam T The type to which the Iterator points.
 * 
 * ## Motivation
 * 
 * Simply put, a class satisfies the concept of a forward iterator if it
 * - can be dereferenced: `*fi` and `fi->` are valid expressions
 * - can be incremented: `++fi`, `fi++`
 * - is default constructible, copyable, moveable, (move-)assignable.
 * - is comparable
 * 
 * Many Standard Library containers provide forward iterators:
 * - `std::vector<T>::%begin()` and `std::vector<T>::%end()` return forward iterators.
 * - `std::list<T>::%begin()` and `std::list<T>::%end()` return forward iterators.
 * - and many more...
 * 
 * Sometimes it is necessary to expose a (pair of) forward iterators in
 * a virtual interface, e.g. to give the user of the interface access to a 
 * set of elements. Consider for example the following virtual interface class
 * that represents a recipe:
 * @snippet forward_iterator.cc recipe
 * You can see that this virtual interface class doesn't specify whether the
 * iterators are `std::vector` or `std::list` forward iterators because the 
 * interface shouldn't prescribe the implementation how it should store the
 * ingredients. Instead this implementation uses the ForwardIterator wrapper
 * class which can wrap a `std::vector` forward iterator, a `std::list` forward
 * iterator or any other forward iterator. Thats exactly the reason that this
 * class exists.
 * 
 * ## Usage
 * @snippet forward_iterator.cc usage
 * 
 * ## Thoughts about implementation
 * Maybe you've asked yourself why the ForwardIterator class isn't a virtual 
 * interface itself. The reason for this is mostly that we cannot declare 
 * - default constructibility
 * - copyability
 * - moveability
 * in a virtual interface. However, all iterators of the standard library
 * fulfill this criteria. Therefore the ForwardIterator class is not a virtual
 * interface itself but a wrapper around any forward iterator. Behind the 
 * scenes type erasure is used to implement this behavior.
 * 
 * @see http://en.cppreference.com/w/cpp/concept/ForwardIterator for an exact
 * description of the forward iterator concept.
 * 
 * @see lf::base::ForwardRange A forward iterator is often combined with this class.
 */
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

  /**
   * @brief Construct wrapper from a copy of a forward iterator.
   * @tparam IteratorImpl The type of the forward iterator.
   * @param iterator The iterator that should be wrapped.
   */
  template <class IteratorImpl, typename = std::iterator_traits<IteratorImpl>>
  ForwardIterator(const IteratorImpl& iterator)
    : wrapper_(std::make_unique<WrapperImpl<IteratorImpl>>(iterator)) {
  }


  /**
   * @brief Construct wrapper from a RValue of a forward iterator.
   * @tparam IteratorImpl  The type of the forward iterator.
   * @param iterator The iterator that should be wrapped.
   */
  template <class IteratorImpl, typename = std::iterator_traits<IteratorImpl>,
  typename = std::enable_if_t<!std::is_convertible<IteratorImpl, ForwardIterator&>::value>>
  ForwardIterator(IteratorImpl&& iterator)
    : wrapper_(std::make_unique<WrapperImpl<IteratorImpl>>(iterator)) {
  }

  
  /**
   * @brief Default constructor of the Forward iterator.
   * 
   * The constructed iterator cannot be dereferenced and is not equal to any
   * other iterator (except for another default constructed iterator).
   */
  ForwardIterator()
    : wrapper_(std::make_unique<WrapperNull>()) {
  }


  /**
   * @brief Copy constructor 
   */
  ForwardIterator(const ForwardIterator& other)
    : wrapper_(other.wrapper_->Clone()) {
  }


  /**
   * @brief Move constructor
   */
  ForwardIterator(ForwardIterator&&) = default;


  /**
   * @brief Move-assignment operator
   */
  ForwardIterator& operator=(ForwardIterator&&) = default;


  /**
   * @brief Standard assignment operator
   */
  ForwardIterator& operator=(const ForwardIterator& rhs) {
    wrapper_ = rhs.wrapper_->Clone();
    return *this;
  }


  /**
   * @brief Compare this Forward iterator with another forward iterator.
   * @param rhs The other iterator to which this iterator should be compared against.
   * @return true if the two iterators are equal, i.e. if they point to the same element.
   * 
   * @note Because this class uses type erasure internally, it is unfortunately
   * not possible to compare two iterators of different types, event if
   * they would be equal:
   * @snippet forward_iterator.cc equality
   */
  bool operator==(const ForwardIterator& rhs) const {
    return wrapper_->Compare(rhs.wrapper_.get());
  }


  /**
   * @brief Shortcut for `!operator==(rhs)`
   */
  bool operator!=(const ForwardIterator& rhs) const {
    return !operator==(rhs);
  }


  /**
   * @brief Dereference this iterator
   */
  T& operator*() const {
    return wrapper_->Dereference();
  }


  /**
   * @brief Dereference this iterator
   */
  T* operator->() const {
    return &wrapper_->Dereference();
  }


  /**
   * @brief (Pre-) increment this iterator
   * @return A reference to the incremented iterator
   */
  ForwardIterator& operator++() {
    wrapper_->operator++();
    return *this;
  }

  /**
   * @brief (Post-) increment this iterator
   * @return A reference to the original (not incremented) iterator.
   */
  ForwardIterator operator++(int) {
    return ForwardIterator(wrapper_->operator++(0));
  }

};

}

/// \cond
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
/// \endcond


#endif // __20db9042f1c04f5c8c173f3e354a915c
