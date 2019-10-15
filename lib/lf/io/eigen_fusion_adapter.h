/**
 * @file
 * @brief Defines boost fusion adapters for eigen. Is used by the boost
 *        spirit and karma libraries from GmshReader and VtkWriter
 * @author Raffael Casagrande
 * @date   2018-07-14 09:45:23
 * @copyright MIT License
 */

#ifndef __a4e6283ee2844d93bc9772c830f33b2d
#define __a4e6283ee2844d93bc9772c830f33b2d

#include <Eigen/Eigen>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/iterator.hpp>
#include <boost/fusion/support/category_of.hpp>
#include <boost/fusion/support/iterator_base.hpp>
#include <boost/fusion/support/tag_of.hpp>
#include <boost/fusion/support/tag_of_fwd.hpp>

// BOOST Fusion Adaption for Eigen fixed size matrices:
//////////////////////////////////////////////////////////////////////////
// The following allows us to use Eigen matrices directly with boost spirit to
// read tuples of double values into fixed size eigen matrices (see
// gmsh_reader.cpp for example usage)

namespace Eigen {
template <class MATRIX>
struct FusionMatrixTag;

template <class MATRIX>
struct FusionIteratorTag;

template <class STRUCT, int POS>
struct FusionIterator
    : boost::fusion::iterator_base<FusionIterator<STRUCT, POS>> {
  static_assert(STRUCT::SizeAtCompileTime != Eigen::Dynamic,
                "Dynamic matrices are not yet supported.");
  static_assert(POS >= 0 & POS <= STRUCT::SizeAtCompileTime,
                "POS has wrong value.");
  using struct_type = STRUCT;
  using index = boost::mpl::int_<POS>;
  using category = boost::fusion::random_access_traversal_tag;

  explicit FusionIterator(STRUCT& str) : struct_(str) {}
  STRUCT& struct_;  // NOLINT
};
}  // namespace Eigen

namespace boost::fusion::traits {
template <class STRUCT, int POS>
struct tag_of<Eigen::FusionIterator<STRUCT, POS>> {
  using type = Eigen::FusionIteratorTag<STRUCT>;
};

template <class T, int ROWS, int COLS, int OPTIONS, int MAX_ROWS, int MAX_COLS>
struct tag_of<Eigen::Matrix<T, ROWS, COLS, OPTIONS, MAX_ROWS, MAX_COLS>> {
  using type = Eigen::FusionMatrixTag<
      Eigen::Matrix<T, ROWS, COLS, OPTIONS, MAX_ROWS, MAX_COLS>>;
};
}  // namespace boost::fusion::traits

namespace boost::fusion::extension {
template <class MATRIX>
struct value_of_impl<Eigen::FusionIteratorTag<MATRIX>> {
  template <class ITERATOR>
  struct apply;
  template <class STRUCT, int N>
  struct apply<Eigen::FusionIterator<STRUCT, N>> {
    using type = typename STRUCT::Scalar;
  };
};

template <class MATRIX>
struct deref_impl<Eigen::FusionIteratorTag<MATRIX>> {
  template <class ITERATOR>
  struct apply;

  template <class STRUCT, int N>
  struct apply<Eigen::FusionIterator<STRUCT, N>> {
    using type =
        typename mpl::if_<is_const<STRUCT>, const typename STRUCT::Scalar&,
                          typename STRUCT::Scalar&>::type;
    static type call(Eigen::FusionIterator<STRUCT, N> const& it) {
      return it.struct_(N);
    }
  };
};

template <class MATRIX>
struct next_impl<Eigen::FusionIteratorTag<MATRIX>> {
  template <class ITERATOR>
  struct apply {
    using struct_type = typename ITERATOR::struct_type;
    using index = typename ITERATOR::index;
    using type = Eigen::FusionIterator<struct_type, index::value + 1>;

    static type call(ITERATOR const& i) { return type(i.struct_); }
  };
};

template <class MATRIX>
struct prior_impl<Eigen::FusionIteratorTag<MATRIX>> {
  template <class ITERATOR>
  struct apply {
    using struct_type = typename ITERATOR::struct_type;
    using index = typename ITERATOR::index;
    using type = Eigen::FusionIterator<struct_type, index::value - 1>;

    static type call(ITERATOR const& i) { return type(i.struct_); }
  };
};

template <class MATRIX>
struct advance_impl<Eigen::FusionIteratorTag<MATRIX>> {
  template <class ITERATOR, class N>
  struct apply {
    using struct_type = typename ITERATOR::struct_type;
    using index = typename ITERATOR::index;
    using type = Eigen::FusionIterator<struct_type, index::value + N::value>;

    static type call(ITERATOR const& i) { return type(i.struct_); }
  };
};

template <class MATRIX>
struct distance_impl<Eigen::FusionIteratorTag<MATRIX>> {
  template <class FIRST, class LAST>
  struct apply : mpl::minus<typename LAST::index, typename FIRST::index> {
    using self = apply<FIRST, LAST>;

    static typename self::type call(FIRST const& /*first*/,
                                    LAST const& /*last*/) {
      return typename self::type();
    }
  };
};

template <class MATRIX>
struct equal_to_impl<Eigen::FusionIteratorTag<MATRIX>> {
  template <class IT1, class IT2>
  struct apply : mpl::equal_to<typename IT1::index, typename IT2::index> {};
};

template <class MATRIX>
struct begin_impl<Eigen::FusionMatrixTag<MATRIX>> {
  template <class SEQUENCE>
  struct apply {
    using type = Eigen::FusionIterator<SEQUENCE, 0>;

    static type call(SEQUENCE& seq) { return type(seq); }
  };
};

template <class MATRIX>
struct end_impl<Eigen::FusionMatrixTag<MATRIX>> {
  template <class SEQUENCE>
  struct apply {
    using type = Eigen::FusionIterator<SEQUENCE, MATRIX::SizeAtCompileTime>;
    static type call(SEQUENCE& seq) { return type(seq); }
  };
};

template <class MATRIX>
struct at_impl<Eigen::FusionMatrixTag<MATRIX>> {
  template <class SEQUENCE, class KEY>
  struct apply;
  template <class SEQUENCE, int N>
  struct apply<SEQUENCE, mpl::int_<N>> {
    using type =
        typename mpl::if_<is_const<SEQUENCE>, const typename MATRIX::Scalar&,
                          typename MATRIX::Scalar&>::type;

    static type call(SEQUENCE& seq) { return seq(N); }
  };
};

template <class MATRIX>
struct value_at_impl<Eigen::FusionMatrixTag<MATRIX>> {
  template <typename Sequence, typename N>
  struct apply;

  template <class SEQUENCE, int N>
  struct apply<SEQUENCE, mpl::int_<N>> {
    using type = typename MATRIX::Scalar;
  };
};

template <class MATRIX>
struct size_impl<Eigen::FusionMatrixTag<MATRIX>> {
  template <class SEQUENCE>
  struct apply : mpl::int_<MATRIX::SizeAtCompileTime> {};
};

template <class MATRIX>
struct category_of_impl<Eigen::FusionMatrixTag<MATRIX>> {
  template <class SEQUENCE>
  struct apply {
    struct type : random_access_traversal_tag {};
  };
};

template <class MATRIX>
struct is_sequence_impl<Eigen::FusionMatrixTag<MATRIX>> {
  template <class T>
  struct apply : mpl::true_ {};
};

template <class MATRIX>
struct is_view_impl<Eigen::FusionMatrixTag<MATRIX>> {
  template <class SEQ>
  struct apply : boost::mpl::false_ {};
};
}  // namespace boost::fusion::extension

#endif  // __a4e6283ee2844d93bc9772c830f33b2d
