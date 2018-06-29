#include "gmsh_reader.h"

#include <fstream>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct_named.hpp>
#include <boost/fusion/include/boost_array.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/fusion/iterator.hpp>
#include <boost/fusion/support/category_of.hpp>
#include <boost/fusion/support/iterator_base.hpp>
#include <boost/fusion/support/tag_of.hpp>
#include <boost/fusion/support/tag_of_fwd.hpp>
#include <boost/mpl/minus.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_binary.hpp>

using size_type = lf::mesh::Mesh::size_type;

// Structures that represent the msh_file:
namespace lf::io {

/// Output the element type onto the console:
std::ostream& operator<<(std::ostream& stream, msh_file::ElementType et) {
  switch (et) {
    default:
      break;
    case msh_file::ElementType::EDGE2:
      stream << "EDGE2";
      break;
    case msh_file::ElementType::TRIA3:
      stream << "TRIA3";
      break;
    case msh_file::ElementType::QUAD4:
      stream << "QUAD4";
      break;
    case msh_file::ElementType::TET4:
      stream << "TET4";
      break;
    case msh_file::ElementType::HEX8:
      stream << "HEX8";
      break;
    case msh_file::ElementType::PRISM6:
      stream << "PRISM6";
      break;
    case msh_file::ElementType::PYRAMID5:
      stream << "PYRAMID5";
      break;
    case msh_file::ElementType::EDGE3:
      stream << "EDGE3";
      break;
    case msh_file::ElementType::TRIA6:
      stream << "TRIA6";
      break;
    case msh_file::ElementType::QUAD9:
      stream << "QUAD9";
      break;
    case msh_file::ElementType::TET10:
      stream << "TET10";
      break;
    case msh_file::ElementType::HEX27:
      stream << "HEX27";
      break;
    case msh_file::ElementType::PRISM18:
      stream << "PRISM18";
      break;
    case msh_file::ElementType::PYRAMID14:
      stream << "PYRAMID14";
      break;
    case msh_file::ElementType::POINT:
      stream << "POINT";
      break;
    case msh_file::ElementType::QUAD8:
      stream << "QUAD8";
      break;
    case msh_file::ElementType::HEX20:
      stream << "HEX20";
      break;
    case msh_file::ElementType::PRISM15:
      stream << "PRISM15";
      break;
    case msh_file::ElementType::PYRAMID13:
      stream << "PYRAMID13";
      break;
    case msh_file::ElementType::TRIA9:
      stream << "TRIA9";
      break;
    case msh_file::ElementType::TRIA10:
      stream << "TRIA10";
      break;
    case msh_file::ElementType::TRIA12:
      stream << "TRIA12";
      break;
    case msh_file::ElementType::TRIA15:
      stream << "TRIA15";
      break;
    case msh_file::ElementType::TRIA15_5:
      stream << "TRIA15_5";
      break;
    case msh_file::ElementType::TRIA21:
      stream << "TRIA21";
      break;
    case msh_file::ElementType::EDGE4:
      stream << "EDGE4";
      break;
    case msh_file::ElementType::EDGE5:
      stream << "EDGE5";
      break;
    case msh_file::ElementType::EDGE6:
      stream << "EDGE6";
      break;
    case msh_file::ElementType::TET20:
      stream << "TET20";
      break;
    case msh_file::ElementType::TET35:
      stream << "TET35";
      break;
    case msh_file::ElementType::TET56:
      stream << "TET56";
      break;
    case msh_file::ElementType::HEX64:
      stream << "HEX64";
      break;
    case msh_file::ElementType::HEX125:
      stream << "HEX125";
      break;
  }
  return stream;
}

/// For debugging purposes: Write the msh_file into a stream
std::ostream& operator<<(std::ostream& stream, const msh_file& mf) {
  stream << "GMSH FILE: Ver. " << mf.VersionNumber
         << (mf.IsBinary ? "(Binary)" : "(Text)")
         << ", size of double = " << mf.DoubleSize << std::endl;
  stream << "======================================================="
         << std::endl;
  stream << "PHYSICAL ENTITIES (Dimension, Number, Name):" << std::endl;
  for (auto pe : mf.PhysicalEntities) {
    stream << "  " << pe.Dimension << "\t , " << pe.Number << "\t , " << pe.Name
           << std::endl;
  }
  stream << "NODES (Number, coords)" << std::endl;
  for (auto n : mf.Nodes) {
    stream << "  " << n.first << "\t , " << n.second.transpose() << std::endl;
  }
  stream << "ELEMENTS (Number, Type, PhysicalEntity Nr, ElementaryEntityNr, "
            "Mesh partitions to which it belongs, Node numbers in it)"
         << std::endl;
  for (auto e : mf.Elements) {
    stream << "  " << e.Number << "\t " << e.Type << '\t' << e.PhysicalEntityNr
           << "\t" << e.ElementaryEntityNr << "\t";
    for (auto p : e.MeshPartitions) stream << p << ", ";
    stream << '\t';
    for (auto n : e.NodeNumbers) stream << n << ", ";
    stream << std::endl;
  }
  stream << std::endl;
  stream << "PERIODIC ENTITIES:" << std::endl;
  for (auto pe : mf.Periodic) {
    std::cout << "  dim=" << pe.Dimension
              << ", slaveNr=" << pe.ElementarySlaveNr
              << ", masterNr=" << pe.ElementaryMasterNr << std::endl;
    for (auto nm : pe.NodeMapping) {
      std::cout << "    " << nm.first << " <-> " << nm.second << std::endl;
    }
  }
  return stream;
}

/// Number of nodes that this element type has
size_type numNodes(msh_file::ElementType et) {
  switch (et) {
    default:
      break;
    case msh_file::ElementType::EDGE2:
      return 2;
    case msh_file::ElementType::TRIA3:
      return 3;
    case msh_file::ElementType::QUAD4:
      return 4;
    case msh_file::ElementType::TET4:
      return 4;
    case msh_file::ElementType::HEX8:
      return 8;
    case msh_file::ElementType::PRISM6:
      return 6;
    case msh_file::ElementType::PYRAMID5:
      return 5;
    case msh_file::ElementType::EDGE3:
      return 3;
    case msh_file::ElementType::TRIA6:
      return 6;
    case msh_file::ElementType::QUAD9:
      return 9;
    case msh_file::ElementType::TET10:
      return 10;
    case msh_file::ElementType::HEX27:
      return 27;
    case msh_file::ElementType::PRISM18:
      return 18;
    case msh_file::ElementType::PYRAMID14:
      return 14;
    case msh_file::ElementType::POINT:
      return 1;
    case msh_file::ElementType::QUAD8:
      return 8;
    case msh_file::ElementType::HEX20:
      return 20;
    case msh_file::ElementType::PRISM15:
      return 15;
    case msh_file::ElementType::PYRAMID13:
      return 13;
    case msh_file::ElementType::TRIA9:
      return 9;
    case msh_file::ElementType::TRIA10:
      return 10;
    case msh_file::ElementType::TRIA12:
      return 12;
    case msh_file::ElementType::TRIA15:
      return 15;
    case msh_file::ElementType::TRIA15_5:
      return 15;
    case msh_file::ElementType::TRIA21:
      return 21;
    case msh_file::ElementType::EDGE4:
      return 4;
    case msh_file::ElementType::EDGE5:
      return 5;
    case msh_file::ElementType::EDGE6:
      return 6;
    case msh_file::ElementType::TET20:
      return 20;
    case msh_file::ElementType::TET35:
      return 35;
    case msh_file::ElementType::TET56:
      return 56;
    case msh_file::ElementType::HEX64:
      return 64;
    case msh_file::ElementType::HEX125:
      return 125;
  }
  LF_VERIFY_MSG(false, "unknown Gmsh element type");
  // Make compiler happy:
  return 0;
}

/// Dimension of the GmshElement type
int dimOf(msh_file::ElementType et) {
  switch (et) {
    case msh_file::ElementType::POINT:
      return 0;
    case msh_file::ElementType::EDGE2:
    case msh_file::ElementType::EDGE3:
    case msh_file::ElementType::EDGE4:
    case msh_file::ElementType::EDGE5:
    case msh_file::ElementType::EDGE6:
      return 1;
    case msh_file::ElementType::TRIA3:
    case msh_file::ElementType::QUAD4:
    case msh_file::ElementType::TRIA6:
    case msh_file::ElementType::QUAD9:
    case msh_file::ElementType::QUAD8:
    case msh_file::ElementType::TRIA9:
    case msh_file::ElementType::TRIA10:
    case msh_file::ElementType::TRIA12:
    case msh_file::ElementType::TRIA15:
    case msh_file::ElementType::TRIA15_5:
    case msh_file::ElementType::TRIA21:
      return 2;
    case msh_file::ElementType::TET4:
    case msh_file::ElementType::HEX8:
    case msh_file::ElementType::PRISM6:
    case msh_file::ElementType::PYRAMID5:
    case msh_file::ElementType::TET10:
    case msh_file::ElementType::HEX27:
    case msh_file::ElementType::PRISM18:
    case msh_file::ElementType::PYRAMID14:
    case msh_file::ElementType::HEX20:
    case msh_file::ElementType::PRISM15:
    case msh_file::ElementType::PYRAMID13:
    case msh_file::ElementType::TET20:
    case msh_file::ElementType::TET35:
    case msh_file::ElementType::TET56:
    case msh_file::ElementType::HEX64:
    case msh_file::ElementType::HEX125:
      return 3;
    default:
      LF_VERIFY_MSG(false, "Unknown GmshElement Type.");
  }
  // Make compiler happy:
  return -1;
}
}  // namespace lf::io

// Boost Fusion Adaptions (needed so boost spirit can parse directly into
// msh_file struct)
//////////////////////////////////////////////////////////////////////////
BOOST_FUSION_ADAPT_STRUCT(lf::io::msh_file::PhysicalEntity,
                          (int, Dimension)(int, Number)(std::string, Name));

BOOST_FUSION_ADAPT_STRUCT(
    lf::io::msh_file::Element,
    (size_type, Number)(lf::io::msh_file::ElementType,
                        Type)(int, PhysicalEntityNr)(int, ElementaryEntityNr)(
        std::vector<int>, MeshPartitions)(std::vector<size_type>, NodeNumbers));

/// To circumvent comma in preprocessor invocation
using nodeMapping_t = std::pair<size_type, size_type>;

BOOST_FUSION_ADAPT_STRUCT(lf::io::msh_file::PeriodicEntity,
                          (int, Dimension)(int, ElementarySlaveNr)(
                              int,
                              ElementaryMasterNr)(std::vector<nodeMapping_t>,
                                                  NodeMapping));

/// To circumvent comma in preprocessor invocation
using nodePair_t = std::pair<size_type, Eigen::Vector3d>;

/// To use msh_file Struct with boost spirit, node that we leave away all
/// header information this is set using attributes.
BOOST_FUSION_ADAPT_STRUCT_NAMED(
    lf::io::msh_file, msh_fileAdapted,
    //(double, VersionNumber)
    //(bool, IsBinary)
    //(int, DoubleSize)
    (std::vector<lf::io::msh_file::PhysicalEntity>,
     PhysicalEntities)(std::vector<nodePair_t>,
                       Nodes)(std::vector<lf::io::msh_file::Element>, Elements)(
        std::vector<lf::io::msh_file::PeriodicEntity>, Periodic));

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

  FusionIterator(STRUCT& str) : struct_(str) {}
  STRUCT& struct_;
};
}  // namespace Eigen

namespace boost::spirit::traits {
/*template<>
struct transform_attribute<hydi::io::msh_file::ElementType, int, qi::domain> {
  using type = int&;
  static int& pre(hydi::io::msh_file::ElementType& d) { return (int&)d; }
  static void post(hydi::io::msh_file::ElementType& dval, const int& attr) {}
  static void fail(hydi::io::msh_file::ElementType&) {}
};*/

template <typename Enum, typename RawValue>
struct assign_to_attribute_from_value<
    Enum, RawValue,
    typename std::enable_if<std::is_enum<Enum>::value &&
                            std::is_same<Enum, RawValue>::value ==
                                false>::type> {
  static void call(RawValue const& raw, Enum& cat) {
    cat = static_cast<Enum>(raw);
  }
};

}  // namespace boost::spirit::traits

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

    static typename self::type call(FIRST const& first, LAST const& last) {
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

namespace lf::io {
namespace /*Anonymous*/ {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

/// A lookup table for boost spirit that can parse an element type
struct gmshElementType : qi::symbols<char, unsigned> {
  gmshElementType() {
    for (auto& et : msh_file::AllElementTypes) {
      add(std::to_string((int)et), (int)et);
    }
  }
};

BOOST_PHOENIX_ADAPT_FUNCTION(int, numNodesAdapted, numNodes, 1);

/// Defines the Grammar of a msh file using boost::spirit
template <class ITERATOR>
struct MshGrammarText
    : qi::grammar<ITERATOR, boost::fusion::adapted::msh_fileAdapted(),
                  ascii::space_type> {
  MshGrammarText(
      qi::rule<ITERATOR, std::pair<size_type, Eigen::Vector3d>()> nodeRule,
      qi::rule<ITERATOR, std::vector<msh_file::Element>(),
               qi::locals<size_type, int, int, int, size_type>>
          elementGroup)
      : MshGrammarText::base_type(start_, "Msh File"),
        node_(nodeRule),
        elementGroup_(elementGroup) {
    using phoenix::push_back;
    using phoenix::reserve;
    using phoenix::val;
    using qi::_val;
    using qi::char_;
    using qi::double_;
    using qi::eps;
    using qi::int_;
    using qi::lexeme;
    using qi::lit;
    using qi::omit;
    using qi::repeat;
    using namespace qi::labels;

    // General Parsers:
    quotedString_ %= lexeme['"' >> +(char_ - '"') >> '"'];
    quotedString_.name("string");
    startComment_ %= !lit("$PhysicalNames") >> !lit("$Nodes") >>
                     !lit("$Elements") >> !lit("$Periodic") >>
                     (lit('$') >> (+(char_ - qi::eol)));
    startComment_.name("Start of Comment");
    comment_ %=
        startComment_[_a = qi::_1] > *(char_ - '$') >> "$End" >> qi::string(_a);
    comment_.name("comment");
    qi::on_error<qi::fail>(comment_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));

    // Physical Entities:
    physicalEntity_ %= int_ > int_ > quotedString_;
    physicalEntity_.name("Physical Entity Entry");
    qi::on_error<qi::fail>(physicalEntity_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));
    physicalEntityGroup_ %= "$PhysicalNames" >
                            omit[int_[reserve(_val, qi::_1), _a = qi::_1]] >
                            repeat(_a)[physicalEntity_] > "$EndPhysicalNames";
    physicalEntityGroup_.name("$Physical Entity Section");
    qi::on_error<qi::fail>(physicalEntityGroup_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));

    // Nodes:
    nodeGroup_ %= "$Nodes" > qi::eol >
                  omit[qi::uint_[reserve(_val, qi::_1), _a = qi::_1]] >
                  qi::eol > repeat(_a)[node_] > -qi::eol > "$EndNodes";
    nodeGroup_.name("$Node Section");
    qi::on_error<qi::fail>(node_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));
    qi::on_error<qi::fail>(nodeGroup_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));

    // Elements:
    qi::on_error<qi::fail>(elementGroup_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));

    // Periodic entities:
    periodicEntityNodeMapping_ =
        omit[qi::uint_[reserve(_val, qi::_1), _a = qi::_1]] >
        repeat(_a)[qi::uint_ > qi::uint_];
    periodicEntityNodeMapping_.name("slave-master node mapping");
    qi::on_error<qi::fail>(periodicEntityNodeMapping_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));
    periodicEntity_ = int_ > int_ > int_ > periodicEntityNodeMapping_;
    periodicEntity_.name("periodic entity");
    qi::on_error<qi::fail>(periodicEntity_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));
    periodicEntityGroup_ = "$Periodic" >
                           omit[qi::uint_[reserve(_val, qi::_1), _a = qi::_1]] >
                           repeat(_a)[periodicEntity_] > "$EndPeriodic";
    periodicEntityGroup_.name("periodic entity section");
    qi::on_error<qi::fail>(periodicEntityGroup_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));

    // The whole file:
    start_ %= *comment_ >> -(physicalEntityGroup_ >> *comment_) >> nodeGroup_ >>
              *comment_ >> elementGroup_ >> *comment_ >>
              -(periodicEntityGroup_ >> *comment_);
    start_.name("beginning of file");
    qi::on_error<qi::fail>(start_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));
  }

  qi::rule<ITERATOR, std::string(), ascii::space_type> quotedString_;
  qi::rule<ITERATOR, std::string()> startComment_;
  qi::rule<ITERATOR, qi::locals<std::string>, ascii::space_type> comment_;

  qi::rule<ITERATOR, msh_file::PhysicalEntity(), ascii::space_type>
      physicalEntity_;
  qi::rule<ITERATOR, std::vector<msh_file::PhysicalEntity>(),
           qi::locals<size_type>, ascii::space_type>
      physicalEntityGroup_;

  qi::rule<ITERATOR, std::pair<size_type, Eigen::Vector3d>()> node_;
  qi::rule<ITERATOR, std::vector<std::pair<size_type, Eigen::Vector3d>>(),
           qi::locals<size_type>>
      nodeGroup_;

  /// locals of elementGroup_ are: (# elements, current element type nr, # tags,
  /// # elements read so far)
  qi::rule<ITERATOR, std::vector<msh_file::Element>(),
           qi::locals<size_type, int, int, int, size_type>>
      elementGroup_;

  qi::rule<ITERATOR, std::vector<std::pair<size_type, size_type>>(),
           qi::locals<size_type>, ascii::space_type>
      periodicEntityNodeMapping_;
  qi::rule<ITERATOR, msh_file::PeriodicEntity(), ascii::space_type>
      periodicEntity_;
  qi::rule<ITERATOR, std::vector<msh_file::PeriodicEntity>(),
           qi::locals<size_type>, ascii::space_type>
      periodicEntityGroup_;

  qi::rule<ITERATOR, boost::fusion::adapted::msh_fileAdapted(),
           ascii::space_type>
      start_;

  struct ErrorHandler {
    template <class, class, class, class>
    struct result {
      using type = void;
    };

    template <class FIRST, class LAST, class ERROR_POS, class WHAT>
    void operator()(FIRST first, LAST last, ERROR_POS errorPos,
                    WHAT what) const {
      std::string input(first, last);
      if (input.length() > 40) input = input.substr(0, 40);
      std::cout << "Error in msh_file! Expecting " << what << " here: \""
                << input << "\"" << std::endl;
    }
  };
  phoenix::function<ErrorHandler> errorHandler_;
};

}  // namespace

const std::vector<msh_file::ElementType> msh_file::AllElementTypes{
    ElementType::EDGE2,     ElementType::TRIA3,     ElementType::QUAD4,
    ElementType::TET4,      ElementType::HEX8,      ElementType::PRISM6,
    ElementType::PYRAMID5,  ElementType::EDGE3,     ElementType::TRIA6,
    ElementType::QUAD9,     ElementType::TET10,     ElementType::HEX27,
    ElementType::PRISM18,   ElementType::PYRAMID14, ElementType::POINT,
    ElementType::QUAD8,     ElementType::HEX20,     ElementType::PRISM15,
    ElementType::PYRAMID13, ElementType::TRIA9,     ElementType::TRIA10,
    ElementType::TRIA12,    ElementType::TRIA15,    ElementType::TRIA15_5,
    ElementType::TRIA21,    ElementType::EDGE4,     ElementType::EDGE5,
    ElementType::EDGE6,     ElementType::TET20,     ElementType::TET35,
    ElementType::TET56,     ElementType::HEX64,     ElementType::HEX125};

msh_file readGmsh_file(std::string filename) {
  // Open file and copy into memory:hydi::io::msh_file
  //////////////////////////////////////////////////////////////////////////
  std::ifstream in(filename, std::ios_base::in);
  if (!in) {
    std::string error("Could not open file ");
    error += filename;
    throw new base::LfException(error);
  }
  std::string storage;
  in.unsetf(std::ios::skipws);  // No white space skipping
  std::copy(std::istream_iterator<char>(in), std::istream_iterator<char>(),
            std::back_inserter(storage));

  // Parse header to determine if we are dealing with ASCII format or binary
  // format + little or big endian:
  //////////////////////////////////////////////////////////////////////////
  msh_file result;
  std::string::const_iterator iter = storage.begin();
  std::string::const_iterator end = storage.end();
  using iterator_t = std::string::const_iterator;

  int one;
  bool successful;
  successful = qi::phrase_parse(
      iter, end,
      qi::lit("$MeshFormat") >>
          qi::double_[phoenix::ref(result.VersionNumber) = qi::_1] >>
          ((qi::lit('0')[phoenix::ref(result.IsBinary) = false] >>
            qi::int_[phoenix::ref(result.DoubleSize) = qi::_1]) |
           (qi::lit('1')[phoenix::ref(result.IsBinary) = true] >>
            qi::int_[phoenix::ref(result.DoubleSize) = qi::_1] >>
            qi::little_dword[phoenix::ref(one) = qi::_1])) >>
          "$EndMeshFormat",
      ascii::space);
  LF_VERIFY_MSG(successful, "Could not read header of file " << filename);
  LF_VERIFY_MSG(result.VersionNumber == 2.2,
                "This GMSH Reader supports only version 2.2 of the mesh file.");
  LF_ASSERT_MSG(result.DoubleSize == 8, "Size of double must be 8.");

  // Parse the rest of the document
  //////////////////////////////////////////////////////////////////////////

  // Setup parsers for node/element sections (which are different depending on
  // binary/non-binary files):
  //
  // Note vec3 has no skipper because it may be used inside lexeme and lexeme
  // can only use parsers without skippers!
  // http://boost-spirit.com/home/2010/02/24/parsing-skippers-and-skipping-parsers/
  // (see comment section)
  qi::rule<iterator_t, Eigen::Vector3d> vec3;
  qi::rule<iterator_t, std::pair<size_type, Eigen::Vector3d>()> node;
  qi::rule<iterator_t, msh_file::Element(), qi::locals<int>> elementText;
  qi::rule<iterator_t, msh_file::Element(msh_file::ElementType, int, int)>
      elementBin;
  qi::rule<iterator_t, std::vector<msh_file::Element>(),
           qi::locals<size_type, int, int, int, size_type>>
      elementGroup;

  using namespace qi::labels;
  using phoenix::reserve;
  using qi::omit;
  using qi::repeat;
  if (result.IsBinary == false) {
    // Text file
    vec3 = qi::double_ >> ' ' >> qi::double_ >> ' ' >> qi::double_;
    node = qi::uint_ >> ' ' >> vec3 >> qi::eol;
    elementText %=
        qi::int_ > ' ' > qi::int_ > ' ' > qi::omit[qi::int_[qi::_a = qi::_1]] >
        ' ' > qi::int_ > ' ' > qi::int_ > ' ' >
        ((qi::eps(_a > 2) >> omit[qi::int_] >> ' ') || qi::eps) >
        qi::repeat(qi::_a - 3)[qi::int_ >> ' '] > (qi::uint_ % ' ') > qi::eol;
    elementGroup %= "$Elements" > qi::eol >
                    qi::omit[qi::uint_[phoenix::reserve(qi::_val, qi::_1),
                                       qi::_a = qi::_1]] > qi::eol >
                    qi::repeat(qi::_a)[elementText] > "$EndElements";
  } else if (result.IsBinary && one == 1) {
    // Binary File Little Endian
    // std::cout << "little endian" << std::endl;
    vec3 %=
        qi::little_bin_double >> qi::little_bin_double >> qi::little_bin_double;
    node %= qi::no_skip[qi::little_dword >> vec3];
    elementBin %= qi::little_dword >> qi::attr(_r1) >> qi::little_dword >>
                  qi::little_dword >>
                  ((qi::eps(_r2 > 2) >> omit[qi::little_dword]) || qi::eps) >>
                  qi::repeat(_r2 - 3)[qi::little_dword] >>
                  qi::repeat(_r3)[qi::little_dword];
    elementGroup %=
        "$Elements" >> qi::eol >> qi::eps[_e = 0] >>
        omit[qi::uint_[reserve(_val, qi::_1), _a = qi::_1]] >>
        qi::eol  // # Elements in total
        >> omit[*((qi::eps(_e < _a) >> qi::little_dword[_b = qi::_1] >>
                   qi::little_dword[_c = qi::_1] >>
                   qi::little_dword[_d = qi::_1]  // elements-header-binary
                   >> repeat(_c)[elementBin(
                          phoenix::static_cast_<msh_file::ElementType>(_b), _d,
                          numNodesAdapted(
                              phoenix::static_cast_<msh_file::ElementType>(
                                  _b)))[phoenix::push_back(_val, qi::_1)]]) >>
                  qi::eps[_e += _c])]  // elements-binary
        >> qi::eol >> "$EndElements";
  } else {
    // std::cout << "big endian" << std::endl;
    // Binary File Big Endian
    vec3 %= qi::big_bin_double >> qi::big_bin_double >> qi::big_bin_double;
    node %= qi::no_skip[qi::big_dword >> vec3];
    elementBin %=
        qi::big_dword >> qi::attr(_r1) >> qi::big_dword >> qi::big_dword >>
        ((qi::eps(_r2 > 2) >> omit[qi::big_dword]) || qi::eps) >>
        qi::repeat(_r2 - 3)[qi::big_dword] >> qi::repeat(_r3)[qi::big_dword];
    elementGroup %=
        "$Elements" >> qi::eol >> qi::eps[_e = 0] >>
        omit[qi::uint_[reserve(_val, qi::_1), _a = qi::_1]] >>
        qi::eol  // # Elements in total
        >> omit[*((qi::eps(_e < _a) >> qi::big_dword[_b = qi::_1] >>
                   qi::big_dword[_c = qi::_1] >>
                   qi::big_dword[_d = qi::_1]  // elements-header-binary
                   >> repeat(_c)[elementBin(
                          phoenix::static_cast_<msh_file::ElementType>(_b), _d,
                          numNodesAdapted(
                              phoenix::static_cast_<msh_file::ElementType>(
                                  _b)))[phoenix::push_back(_val, qi::_1)]]) >>
                  qi::eps[_e += _c])]  // elements-binary
        >> qi::eol >> "$EndElements";
  }

  /// Name the elements for better error parsing:
  vec3.name("vec3");
  node.name("node");
  elementText.name("element");
  elementBin.name("element");
  elementGroup.name("ElementSection");

  // Finally parse everything:
  MshGrammarText<iterator_t> mshGrammar(node, elementGroup);
  bool r = qi::phrase_parse(iter, end, mshGrammar, ascii::space, result);

  // if (r && iter == end) std::cout << "Parsing succeeded" << std::endl;
  // else if (r) std::cout << "Parsing partially succeeded" << std::endl;
  // std::cout << result << std::endl;

  LF_VERIFY_MSG(r, "Could not parse file " << filename);
  LF_VERIFY_MSG(iter == end, "Could not parse all of file " << filename);

  return result;
}

GmshReader::GmshReader(std::unique_ptr<mesh::MeshFactory> factory,
                       const ::lf::io::MshFile& msh_file)
    : mesh_factory_(std::move(factory)) {
  // 1) Check Gmsh_file and initialize
  //////////////////////////////////////////////////////////////////////////

  /*if (msh_file.Nodes.size() == 0)
    LOGGER_ENTRY(logger_, "Warning: The GMSH meshfile " << " does not contain
  any nodes.", 3); if (msh_file.Elements.size() == 0) LOGGER_ENTRY(logger_,
  "Warning: The GMSH meshfile " << " does not contain any elements.", 3);*/

  /// GmshNodeNr2InsertionIndex[i] = j means: msh_file.Nodes[j].first = i.
  std::vector<size_type> GmshNodeNr2InsertionIndex(msh_file.Nodes.size() + 1,
                                                   size_type(-1));

  /// InsertionIndex2GmshElementNr[i] = j means: msh_file.Elements[j] is the
  /// i-th element inserted into the mesh.
  std::vector<size_type> InsertionIndex2GmshElementNr;
  InsertionIndex2GmshElementNr.reserve(
      msh_file.Elements.size());  // educated guess

  /// BoundarySegmentInsertionIndex2GmshElementNr[i] = j means:
  /// msh_file.Elements[j] is the i-th boundary segment inserted by the grid
  /// factory.
  std::vector<size_type> BoundarySegmentInsertionIndex2GmshElementNr;

  /// nodeInsertionIndex2GmshElementNr[i] = j means: Node with insertion index i
  /// is the same as msh_file.Elements[j]
  // if j==size_type(-1) the node with insertion index i has no
  // corresponding element.
  std::vector<size_type> nodeInsertionIndex2GmshElementNr;

  // 2) Add the nodes:
  //////////////////////////////////////////////////////////////////////////
  bool loggedWarning = false;
  for (int i = 0; i < msh_file.Nodes.size(); ++i) {
    if (msh_file.Nodes[i].first >= GmshNodeNr2InsertionIndex.size()) {
      GmshNodeNr2InsertionIndex.resize(msh_file.Nodes[i].first + 1,
                                       size_type(-1));
    }
    size_type insertionIndex;
    insertionIndex =
        gridFactory.insertNode(stripVector<dimWorld>(msh_file.Nodes[i].second));

    GmshNodeNr2InsertionIndex[msh_file.Nodes[i].first] = i;
    HYDI_ASSERT_MSG(insertionIndex == i,
                    "Insertion Index return by grid is not consecutive.");
    if (dimWorld < 3 && loggedWarning == false &&
        msh_file.Nodes[i].second.z() != 0) {
      LOGGER_ENTRY(logger_,
                   "Warning: The GMSH file "
                       << " contains nodes where the z-component is not zero. "
                       << "The z-component will just be ignored because the "
                          "mesh into which we read has dimWorld!= 3",
                   2);
      loggedWarning = true;
    }
  }

  nodeInsertionIndex2GmshElementNr.resize(msh_file.Nodes.size(), size_type(-1));

  // 3) Add the elements, boundary segments and physical entity numbers of
  // nodes:
  //////////////////////////////////////////////////////////////////////////
  std::vector<size_type> nodes;
  nodes.reserve(8);
  loggedWarning = false;
  bool thereExistsNodeWithPhysicalEntityNr = false;
  for (int i = 0; i < msh_file.Elements.size(); ++i) {
    auto& element = msh_file.Elements[i];
    HYDI_ASSERT_MSG(element.NodeNumbers.size() == numNodes(element.Type),
                    "Element " << element.Number << " is of type "
                               << element.Type << " but has "
                               << element.NodeNumbers.size() << " nodes.");

    if (i > 0 && msh_file.Elements[i - 1].NodeNumbers == element.NodeNumbers) {
      // There are two elements with exactly the same node numbers, as far as I
      // understand this only happens when a mesh element belongs to two or more
      // physical entities. So check that all the remaining information is the
      // same:
      HYDI_ASSERT(msh_file.Elements[i - 1].ElementaryEntityNr ==
                  element.ElementaryEntityNr);
      HYDI_ASSERT(msh_file.Elements[i - 1].MeshPartitions ==
                  element.MeshPartitions);
      HYDI_ASSERT(msh_file.Elements[i - 1].Type == element.Type);

      // we ignore all elements that appear for the second/third/... times:
      continue;
    }

    auto dimElement = dimOf(element.Type);
    if (dimElement == 0) {
      nodeInsertionIndex2GmshElementNr
          [GmshNodeNr2InsertionIndex[element.NodeNumbers[0]]] = i;
      thereExistsNodeWithPhysicalEntityNr = true;
      continue;
    }

    if (dimMesh == 3 && dimElement == 1) continue;  // ignore edges in 3D case
    if (dimElement == dimMesh - 1 && element.PhysicalEntityNr == 0)
      continue;  // ignore faces with physicalEntity=0

    HYDI_ASSERT_MSG(dimElement <= dimMesh, "GMSH file contains elements of dim "
                                               << dimElement
                                               << " but dimMesh = " << dimMesh);

    // Get Geometry type:
    base::GeometryType geomType;
    switch (element.Type) {
      case msh_file::ElementType::EDGE2:
        geomType = base::GeometryType::SEGMENT2;
        break;
      case msh_file::ElementType::TRIA3:
        geomType = base::GeometryType::TRIA3;
        break;
      case msh_file::ElementType::QUAD4:
        geomType = base::GeometryType::QUAD4;
        break;
      case msh_file::ElementType::TET4:
        geomType = base::GeometryType::TET4;
        break;
      case msh_file::ElementType::PYRAMID5:
        geomType = base::GeometryType::PYRAMID5;
        break;
      case msh_file::ElementType::PRISM6:
        geomType = base::GeometryType::PRISM6;
        break;
      case msh_file::ElementType::HEX8:
        geomType = base::GeometryType::HEX8;
        break;
      default:
        HYDI_VERIFY_MSG(false,
                        "Gmsh Element type "
                            << element.Type
                            << " is not (yet) supported by Gmsh reader.");
        break;
    }

    // create nodes (using mapping):
    for (int iNode = 0; iNode < element.NodeNumbers.size(); ++iNode) {
      nodes.push_back(GmshNodeNr2InsertionIndex[element.NodeNumbers[iNode]]);
    }
    if (dimElement == dimMesh) {
      // Create Volume element:
      gridFactory.insertElement(geomType, nodes);
      InsertionIndex2GmshElementNr.push_back(i);
    } else if (dimElement == dimMesh - 1) {
      // Create Boundary element:
      gridFactory.insertBoundarySegment(geomType, nodes);
      BoundarySegmentInsertionIndex2GmshElementNr.push_back(i);
    }

    nodes.clear();
  }

  // 4) create grid:
  //////////////////////////////////////////////////////////////////////////
  HYDI_VERIFY_MSG(
      InsertionIndex2GmshElementNr.size() > 0,
      "The GMSH Msh file does not contain any elements of dim=dimMesh="
          << dimMesh);
  grid_ = gridFactory.createGrid();

  // 5) Assign regionID to elements and boundary patch id to boundary
  // intersections.:
  //////////////////////////////////////////////////////////////////////////
  auto gridView = grid_->levelView(0);
  auto gdsRegionID =
      grid::utils::create_SingleCodimGDS<0, size_type>(gridView, 0);
  size_type maxRegionID = 0;
  auto bdsPatchID = grid::utils::create_BoundaryDataSet<size_type>(grid_, 0);
  size_type maxPatchID = 0;
  for (auto& e : gridView.template entities<0>()) {
    size_type regionID =
        msh_file
            .Elements[InsertionIndex2GmshElementNr[gridFactory.insertionIndex(
                e)]]
            .PhysicalEntityNr;
    gdsRegionID->data(e, 0) = regionID;
    if (regionID > maxRegionID) maxRegionID = regionID;
    if (BoundarySegmentInsertionIndex2GmshElementNr.size() > 0) {
      for (auto& i : gridView.intersections(e)) {
        if (i.boundary() == false) continue;
        auto insertionIndex = gridFactory.insertionIndex(i);
        if (insertionIndex >=
            BoundarySegmentInsertionIndex2GmshElementNr.size())
          continue;
        size_type patchID =
            msh_file
                .Elements[BoundarySegmentInsertionIndex2GmshElementNr
                              [insertionIndex]]
                .PhysicalEntityNr;
        bdsPatchID->data(i) = patchID;
        if (patchID > maxPatchID) maxPatchID = patchID;
      }
    }
  }
  egv_ = grid::utils::create_EnhancedGridView(gridView, gdsRegionID, bdsPatchID,
                                              maxRegionID, maxPatchID);

  // 6) Assign the physical entity nr to every node:
  //////////////////////////////////////////////////////////////////////////
  gdsPhysicalNodeNr_ =
      grid::utils::create_SingleCodimGDS<dimMesh, size_type>(gridView, 0);
  size_type maxPhysicalNodeNr = 0;
  if (thereExistsNodeWithPhysicalEntityNr) {
    for (auto& n : gridView.template entities<dimMesh>()) {
      size_type gmshNr =
          nodeInsertionIndex2GmshElementNr[gridFactory.insertionIndex(n)];
      if (gmshNr == size_type(-1)) continue;
      size_type physicalNodeNr = msh_file.Elements[gmshNr].PhysicalEntityNr;
      gdsPhysicalNodeNr_->data(n, 0) = physicalNodeNr;
      if (physicalNodeNr > maxPhysicalNodeNr)
        maxPhysicalNodeNr = physicalNodeNr;
    }
  }

  // 7) Create mapping physicalEntityNr <-> physicalEntityName:
  //////////////////////////////////////////////////////////////////////////

  regionID2Name_.resize(maxRegionID + 1);
  bPatchID2Name_.resize(maxPatchID + 1);
  physicalNodeNr2PhysicalNodeName_.resize(maxPhysicalNodeNr + 1);

  for (auto pe : msh_file.PhysicalEntities) {
    if (pe.Name.empty()) continue;
    if (pe.Dimension == 0) {
      physicalNodeNr2PhysicalNodeName_[pe.Number] = pe.Name;
      physicalNodeName2physicalNodeNr_[pe.Name] = pe.Number;
    } else if (pe.Dimension == dimMesh - 1) {
      bPatchID2Name_[pe.Number] = pe.Name;
      bPatchName2Id_[pe.Name] = pe.Number;
    } else if (pe.Dimension == dimMesh) {
      regionID2Name_[pe.Number] = pe.Name;
      regionName2Id_[pe.Name] = pe.Number;
    }
  }

  if (msh_file.Periodic.size() > 0) {
    LOGGER_ENTRY(logger_,
                 "WARNING: GMSH File  contains periodic boundary relations "
                 "between elements. These are ignored by GmshReader.",
                 3);
  }
}

}  // namespace lf::io
