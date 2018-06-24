#ifndef __661aa07c70a04907a44d828a353e6537
#define __661aa07c70a04907a44d828a353e6537

#include <iterator>
#include <vector>
#include "forward_iterator.h"

namespace lf::base {

namespace internal {}

/**
 * @brief A pair of ForwardIterator's that make up a set of elements of type T
 * @tparam T The type of elements contained in the range.
 *
 * ### Motivation
 */
template <class T>
class ForwardRange {
 protected:
  class WrapperInterface {
   protected:
    WrapperInterface() = default;

   public:
    WrapperInterface(const WrapperInterface&) = delete;
    WrapperInterface(WrapperInterface&&) = delete;
    WrapperInterface& operator=(const WrapperInterface&) = delete;
    WrapperInterface& operator=(WrapperInterface&&) = delete;

    virtual std::unique_ptr<WrapperInterface> Clone() const = 0;
    virtual ForwardIterator<T> begin() const = 0;
    virtual ForwardIterator<T> end() const = 0;

    virtual ~WrapperInterface() = default;
  };

  template <class Inner,
            class = typename std::enable_if<
                std::is_convertible<decltype(std::declval<Inner>().begin()),
                                    ForwardIterator<T>>::value>::type,
            class = typename std::enable_if<
                std::is_convertible<decltype(std::declval<Inner>().end()),
                                    ForwardIterator<T>>::value>::type,
            class =
                typename std::enable_if<!std::is_reference<Inner>::value>::type>
  class OwningImpl : public virtual WrapperInterface {
   protected:
    Inner inner_;

   public:
    explicit OwningImpl(Inner&& inner) : inner_(inner) {}

    std::unique_ptr<WrapperInterface> Clone() const override {
      Inner copy = inner_;
      return std::unique_ptr<WrapperInterface>(new OwningImpl(std::move(copy)));
    }

    ForwardIterator<T> begin() const override { return inner_.begin(); }

    ForwardIterator<T> end() const override { return inner_.end(); }
  };

  template <class Inner,
            class = typename std::enable_if<!std::is_const<Inner>::value>::type>
  class ConstReferenceImpl : public virtual WrapperInterface {
   protected:
    const Inner& inner_;

   public:
    explicit ConstReferenceImpl(const Inner& inner) : inner_(inner) {}

    std::unique_ptr<WrapperInterface> Clone() const override {
      return std::unique_ptr<WrapperInterface>(new ConstReferenceImpl(inner_));
    }

    ForwardIterator<T> begin() const override { return inner_.begin(); }
    ForwardIterator<T> end() const override { return inner_.end(); }
  };

  class InitializerListImpl : public virtual WrapperInterface {
   protected:
    std::vector<std::remove_const_t<T>> list_{};

    InitializerListImpl(const InitializerListImpl& other)
        : list_{other.list_} {}

   public:
    InitializerListImpl(
        const std::initializer_list<std::remove_const_t<T>>& list)
        : list_(list) {}
    InitializerListImpl(InitializerListImpl&&) = delete;
    InitializerListImpl& operator=(const InitializerListImpl&) = delete;
    InitializerListImpl& operator=(InitializerListImpl&&) = delete;
    ~InitializerListImpl() override = default;

    std::unique_ptr<WrapperInterface> Clone() const override {
      return std::unique_ptr<WrapperInterface>(new InitializerListImpl(*this));
    }

    ForwardIterator<T> begin() const override { return list_.begin(); }
    ForwardIterator<T> end() const override { return list_.end(); }
  };

  class IteratorPairImpl : public virtual WrapperInterface {
   protected:
    ForwardIterator<T> begin_;
    ForwardIterator<T> end_;

   public:
    explicit IteratorPairImpl(ForwardIterator<T> begin, ForwardIterator<T> end)
        : begin_(std::move(begin)), end_(std::move(end)) {}

    std::unique_ptr<WrapperInterface> Clone() const override {
      return std::unique_ptr<WrapperInterface>(
          new IteratorPairImpl(begin_, end_));
    }

    ForwardIterator<T> begin() const override { return begin_; }
    ForwardIterator<T> end() const override { return end_; }
  };

  std::unique_ptr<WrapperInterface> wrapper_;

 public:
  ForwardRange(const ForwardRange& rhs) : wrapper_(rhs->Clone()) {}
  ForwardRange(ForwardRange&&) noexcept = default;

  template <class C, typename = typename std::enable_if<std::is_same<
                         typename std::iterator_traits<decltype(
                             std::declval<C>().begin())>::iterator_category,
                         std::forward_iterator_tag>::value>::type>
  explicit ForwardRange(const C& forward_range)
      : wrapper_(new ConstReferenceImpl(forward_range)) {}

  ForwardRange(std::initializer_list<std::remove_const_t<T>> initializer_list)
      : wrapper_(new InitializerListImpl(std::move(initializer_list))) {}

  ForwardRange& operator=(const ForwardRange& rhs) {
    wrapper_ = rhs.wrapper_->Clone();
  }
  ForwardRange& operator=(ForwardRange&&) noexcept = default;

  ~ForwardRange() = default;

  ForwardRange(ForwardIterator<T> begin, ForwardIterator<T> end)
      : wrapper_(new IteratorPairImpl(std::move(begin), std::move(end))) {}

  ForwardIterator<T> begin() const { return wrapper_->begin(); }
  ForwardIterator<T> end() const { return wrapper_->end(); }
};

}  // namespace lf::base

#endif  // __661aa07c70a04907a44d828a353e6537
