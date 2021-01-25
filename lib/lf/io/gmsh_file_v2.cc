/**
 * @file
 * @brief Implementation of the logic to parse a GmshFile Version 2.0 into a
 * in-memory data structure
 * @author Raffael Casagrande
 * @date   2019-08-19 03:30:38
 * @copyright MIT License
 */

#include "gmsh_file_v2.h"

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

#include "eigen_fusion_adapter.h"

using size_type = lf::mesh::Mesh::size_type;

namespace lf::io {

/// Output the element type onto the console:
std::ostream& operator<<(std::ostream& stream, GMshFileV2::ElementType et) {
  switch (et) {
    default:
      break;
    case GMshFileV2::ElementType::EDGE2:
      stream << "EDGE2";
      break;
    case GMshFileV2::ElementType::TRIA3:
      stream << "TRIA3";
      break;
    case GMshFileV2::ElementType::QUAD4:
      stream << "QUAD4";
      break;
    case GMshFileV2::ElementType::TET4:
      stream << "TET4";
      break;
    case GMshFileV2::ElementType::HEX8:
      stream << "HEX8";
      break;
    case GMshFileV2::ElementType::PRISM6:
      stream << "PRISM6";
      break;
    case GMshFileV2::ElementType::PYRAMID5:
      stream << "PYRAMID5";
      break;
    case GMshFileV2::ElementType::EDGE3:
      stream << "EDGE3";
      break;
    case GMshFileV2::ElementType::TRIA6:
      stream << "TRIA6";
      break;
    case GMshFileV2::ElementType::QUAD9:
      stream << "QUAD9";
      break;
    case GMshFileV2::ElementType::TET10:
      stream << "TET10";
      break;
    case GMshFileV2::ElementType::HEX27:
      stream << "HEX27";
      break;
    case GMshFileV2::ElementType::PRISM18:
      stream << "PRISM18";
      break;
    case GMshFileV2::ElementType::PYRAMID14:
      stream << "PYRAMID14";
      break;
    case GMshFileV2::ElementType::POINT:
      stream << "POINT";
      break;
    case GMshFileV2::ElementType::QUAD8:
      stream << "QUAD8";
      break;
    case GMshFileV2::ElementType::HEX20:
      stream << "HEX20";
      break;
    case GMshFileV2::ElementType::PRISM15:
      stream << "PRISM15";
      break;
    case GMshFileV2::ElementType::PYRAMID13:
      stream << "PYRAMID13";
      break;
    case GMshFileV2::ElementType::TRIA9:
      stream << "TRIA9";
      break;
    case GMshFileV2::ElementType::TRIA10:
      stream << "TRIA10";
      break;
    case GMshFileV2::ElementType::TRIA12:
      stream << "TRIA12";
      break;
    case GMshFileV2::ElementType::TRIA15:
      stream << "TRIA15";
      break;
    case GMshFileV2::ElementType::TRIA15_5:
      stream << "TRIA15_5";
      break;
    case GMshFileV2::ElementType::TRIA21:
      stream << "TRIA21";
      break;
    case GMshFileV2::ElementType::EDGE4:
      stream << "EDGE4";
      break;
    case GMshFileV2::ElementType::EDGE5:
      stream << "EDGE5";
      break;
    case GMshFileV2::ElementType::EDGE6:
      stream << "EDGE6";
      break;
    case GMshFileV2::ElementType::TET20:
      stream << "TET20";
      break;
    case GMshFileV2::ElementType::TET35:
      stream << "TET35";
      break;
    case GMshFileV2::ElementType::TET56:
      stream << "TET56";
      break;
    case GMshFileV2::ElementType::HEX64:
      stream << "HEX64";
      break;
    case GMshFileV2::ElementType::HEX125:
      stream << "HEX125";
      break;
  }
  return stream;
}

/// For debugging purposes: Write the MshFile into a stream
std::ostream& operator<<(std::ostream& stream, const GMshFileV2& mf) {
  stream << "GMSH FILE: Ver. " << mf.VersionNumber
         << (mf.IsBinary ? "(Binary)" : "(Text)")
         << ", size of double = " << mf.DoubleSize << std::endl;
  stream << "======================================================="
         << std::endl;
  stream << "PHYSICAL ENTITIES (Dimension, Number, Name):" << std::endl;
  for (const auto& pe : mf.PhysicalEntities) {
    stream << "  " << pe.Dimension << "\t , " << pe.Number << "\t , " << pe.Name
           << std::endl;
  }
  stream << "NODES (Number, coords)" << std::endl;
  for (const auto& n : mf.Nodes) {
    stream << "  " << n.first << "\t , " << n.second.transpose() << std::endl;
  }
  stream << "ELEMENTS (Number, Type, PhysicalEntity Nr, ElementaryEntityNr, "
            "Mesh partitions to which it belongs, Node numbers in it)"
         << std::endl;
  for (const auto& e : mf.Elements) {
    stream << "  " << e.Number << "\t " << e.Type << '\t' << e.PhysicalEntityNr
           << "\t" << e.ElementaryEntityNr << "\t";
    for (auto p : e.MeshPartitions) {
      stream << p << ", ";
    }
    stream << '\t';
    for (auto n : e.NodeNumbers) {
      stream << n << ", ";
    }
    stream << std::endl;
  }
  stream << std::endl;
  stream << "PERIODIC ENTITIES:" << std::endl;
  for (const auto& pe : mf.Periodic) {
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
size_type NumNodes(GMshFileV2::ElementType et) {
  switch (et) {
    default:
      break;
    case GMshFileV2::ElementType::EDGE2:
      return 2;
    case GMshFileV2::ElementType::TRIA3:
      return 3;
    case GMshFileV2::ElementType::QUAD4:
    case GMshFileV2::ElementType::TET4:
      return 4;
    case GMshFileV2::ElementType::HEX8:
      return 8;
    case GMshFileV2::ElementType::PRISM6:
      return 6;
    case GMshFileV2::ElementType::PYRAMID5:
      return 5;
    case GMshFileV2::ElementType::EDGE3:
      return 3;
    case GMshFileV2::ElementType::TRIA6:
      return 6;
    case GMshFileV2::ElementType::QUAD9:
      return 9;
    case GMshFileV2::ElementType::TET10:
      return 10;
    case GMshFileV2::ElementType::HEX27:
      return 27;
    case GMshFileV2::ElementType::PRISM18:
      return 18;
    case GMshFileV2::ElementType::PYRAMID14:
      return 14;
    case GMshFileV2::ElementType::POINT:
      return 1;
    case GMshFileV2::ElementType::QUAD8:
      return 8;
    case GMshFileV2::ElementType::HEX20:
      return 20;
    case GMshFileV2::ElementType::PRISM15:
      return 15;
    case GMshFileV2::ElementType::PYRAMID13:
      return 13;
    case GMshFileV2::ElementType::TRIA9:
      return 9;
    case GMshFileV2::ElementType::TRIA10:
      return 10;
    case GMshFileV2::ElementType::TRIA12:
      return 12;
    case GMshFileV2::ElementType::TRIA15:
    case GMshFileV2::ElementType::TRIA15_5:
      return 15;
    case GMshFileV2::ElementType::TRIA21:
      return 21;
    case GMshFileV2::ElementType::EDGE4:
      return 4;
    case GMshFileV2::ElementType::EDGE5:
      return 5;
    case GMshFileV2::ElementType::EDGE6:
      return 6;
    case GMshFileV2::ElementType::TET20:
      return 20;
    case GMshFileV2::ElementType::TET35:
      return 35;
    case GMshFileV2::ElementType::TET56:
      return 56;
    case GMshFileV2::ElementType::HEX64:
      return 64;
    case GMshFileV2::ElementType::HEX125:
      return 125;
  }
  LF_VERIFY_MSG(false, "unknown Gmsh element type");
  // Make compiler happy:
  return 0;
}

base::RefEl RefElOf(GMshFileV2::ElementType et) {
  switch (et) {
    case GMshFileV2::ElementType::POINT:
      return base::RefEl::kPoint();

    case GMshFileV2::ElementType::EDGE2:
    case GMshFileV2::ElementType::EDGE3:
    case GMshFileV2::ElementType::EDGE4:
    case GMshFileV2::ElementType::EDGE5:
    case GMshFileV2::ElementType::EDGE6:
      return base::RefEl::kSegment();

    case GMshFileV2::ElementType::TRIA3:
    case GMshFileV2::ElementType::TRIA6:
    case GMshFileV2::ElementType::TRIA9:
    case GMshFileV2::ElementType::TRIA10:
    case GMshFileV2::ElementType::TRIA12:
    case GMshFileV2::ElementType::TRIA15:
    case GMshFileV2::ElementType::TRIA15_5:
    case GMshFileV2::ElementType::TRIA21:
      return base::RefEl::kTria();

    case GMshFileV2::ElementType::QUAD4:
    case GMshFileV2::ElementType::QUAD8:
    case GMshFileV2::ElementType::QUAD9:
      return base::RefEl::kQuad();

    case GMshFileV2::ElementType::TET4:
    case GMshFileV2::ElementType::HEX8:
    case GMshFileV2::ElementType::PRISM6:
    case GMshFileV2::ElementType::PYRAMID5:
    case GMshFileV2::ElementType::TET10:
    case GMshFileV2::ElementType::HEX27:
    case GMshFileV2::ElementType::PRISM18:
    case GMshFileV2::ElementType::PYRAMID14:
    case GMshFileV2::ElementType::HEX20:
    case GMshFileV2::ElementType::PRISM15:
    case GMshFileV2::ElementType::PYRAMID13:
    case GMshFileV2::ElementType::TET20:
    case GMshFileV2::ElementType::TET35:
    case GMshFileV2::ElementType::TET56:
    case GMshFileV2::ElementType::HEX64:
    case GMshFileV2::ElementType::HEX125:
    default:
      LF_VERIFY_MSG(
          false, "Reference element not supported for GmshElement type " << et);
  }
}

/// Dimension of the GmshElement type
int DimOf(GMshFileV2::ElementType et) {
  switch (et) {
    case GMshFileV2::ElementType::POINT:
      return 0;
    case GMshFileV2::ElementType::EDGE2:
    case GMshFileV2::ElementType::EDGE3:
    case GMshFileV2::ElementType::EDGE4:
    case GMshFileV2::ElementType::EDGE5:
    case GMshFileV2::ElementType::EDGE6:
      return 1;
    case GMshFileV2::ElementType::TRIA3:
    case GMshFileV2::ElementType::QUAD4:
    case GMshFileV2::ElementType::TRIA6:
    case GMshFileV2::ElementType::QUAD9:
    case GMshFileV2::ElementType::QUAD8:
    case GMshFileV2::ElementType::TRIA9:
    case GMshFileV2::ElementType::TRIA10:
    case GMshFileV2::ElementType::TRIA12:
    case GMshFileV2::ElementType::TRIA15:
    case GMshFileV2::ElementType::TRIA15_5:
    case GMshFileV2::ElementType::TRIA21:
      return 2;
    case GMshFileV2::ElementType::TET4:
    case GMshFileV2::ElementType::HEX8:
    case GMshFileV2::ElementType::PRISM6:
    case GMshFileV2::ElementType::PYRAMID5:
    case GMshFileV2::ElementType::TET10:
    case GMshFileV2::ElementType::HEX27:
    case GMshFileV2::ElementType::PRISM18:
    case GMshFileV2::ElementType::PYRAMID14:
    case GMshFileV2::ElementType::HEX20:
    case GMshFileV2::ElementType::PRISM15:
    case GMshFileV2::ElementType::PYRAMID13:
    case GMshFileV2::ElementType::TET20:
    case GMshFileV2::ElementType::TET35:
    case GMshFileV2::ElementType::TET56:
    case GMshFileV2::ElementType::HEX64:
    case GMshFileV2::ElementType::HEX125:
      return 3;
    default:
      LF_VERIFY_MSG(false, "Unknown GmshElement Type.");
  }
  // Make compiler happy:
  return -1;
}
}  // namespace lf::io

// Boost Fusion Adaptions (needed so boost spirit can parse directly into
// GMshFileV2 struct)
//////////////////////////////////////////////////////////////////////////
BOOST_FUSION_ADAPT_STRUCT(lf::io::GMshFileV2::PhysicalEntity,
                          (int, Dimension)(int, Number)(std::string, Name));

BOOST_FUSION_ADAPT_STRUCT(
    lf::io::GMshFileV2::Element,
    (size_type, Number)(lf::io::GMshFileV2::ElementType,
                        Type)(int, PhysicalEntityNr)(int, ElementaryEntityNr)(
        std::vector<int>, MeshPartitions)(std::vector<size_type>, NodeNumbers));

/// To circumvent comma in preprocessor invocation
using nodeMapping_t = std::pair<size_type, size_type>;

BOOST_FUSION_ADAPT_STRUCT(lf::io::GMshFileV2::PeriodicEntity,
                          (int, Dimension)(int, ElementarySlaveNr)(
                              int,
                              ElementaryMasterNr)(std::vector<nodeMapping_t>,
                                                  NodeMapping));

/// To circumvent comma in preprocessor invocation
using nodePair_t = std::pair<size_type, Eigen::Vector3d>;

/// To use MshFile Struct with boost spirit, node that we leave away all
/// header information this is set using attributes.
// NOLINTNEXTLINE
BOOST_FUSION_ADAPT_STRUCT_NAMED(
    lf::io::GMshFileV2, MshFileAdapted,
    //(double, VersionNumber)
    //(bool, IsBinary)
    //(int, DoubleSize)
    (std::vector<lf::io::GMshFileV2::PhysicalEntity>,
     PhysicalEntities)(std::vector<nodePair_t>, Nodes)(
        std::vector<lf::io::GMshFileV2::Element>,
        Elements)(std::vector<lf::io::GMshFileV2::PeriodicEntity>, Periodic));

/// \cond
namespace boost::spirit::traits {
/*template<>
struct transform_attribute<hydi::io::MshFile::ElementType, int, qi::domain> {
  using type = int&;
  static int& pre(hydi::io::MshFile::ElementType& d) { return (int&)d; }
  static void post(hydi::io::MshFile::ElementType& dval, const int& attr) {}
  static void fail(hydi::io::MshFile::ElementType&) {}
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
/// \endcond

namespace lf::io {
namespace /*Anonymous*/ {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

/// A lookup table for boost spirit that can parse an element type
struct gmshElementType : qi::symbols<char, unsigned> {
  gmshElementType() {
    for (const auto& et : GMshFileV2::AllElementTypes) {
      add(std::to_string(static_cast<int>(et)), static_cast<int>(et));
    }
  }
};

// NOLINTNEXTLINE
BOOST_PHOENIX_ADAPT_FUNCTION(int, numNodesAdapted, NumNodes, 1);

/// Defines the Grammar of a msh file using boost::spirit
template <class ITERATOR>
struct MshGrammarText
    : qi::grammar<ITERATOR, boost::fusion::adapted::MshFileAdapted(),
                  ascii::space_type> {
  MshGrammarText(
      qi::rule<ITERATOR, std::pair<size_type, Eigen::Vector3d>()> nodeRule,
      qi::rule<ITERATOR, std::vector<GMshFileV2::Element>(),
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
    using qi::labels::_1;
    using qi::labels::_2;
    using qi::labels::_3;
    using qi::labels::_4;
    using qi::labels::_a;

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
    physicalEntity_ %= int_ > int_ > quotedString_;  // NOLINT
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
        repeat(_a)[qi::uint_ > qi::uint_];  // NOLINT
    periodicEntityNodeMapping_.name("slave-master node mapping");
    qi::on_error<qi::fail>(periodicEntityNodeMapping_,
                           errorHandler_(qi::_1, qi::_2, qi::_3, qi::_4));
    periodicEntity_ =
        int_ > int_ > int_ > periodicEntityNodeMapping_;  // NOLINT
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

 private:
  qi::rule<ITERATOR, std::string(), ascii::space_type> quotedString_;
  qi::rule<ITERATOR, std::string()> startComment_;
  qi::rule<ITERATOR, qi::locals<std::string>, ascii::space_type> comment_;

  qi::rule<ITERATOR, GMshFileV2::PhysicalEntity(), ascii::space_type>
      physicalEntity_;
  qi::rule<ITERATOR, std::vector<GMshFileV2::PhysicalEntity>(),
           qi::locals<size_type>, ascii::space_type>
      physicalEntityGroup_;

  qi::rule<ITERATOR, std::pair<size_type, Eigen::Vector3d>()> node_;
  qi::rule<ITERATOR, std::vector<std::pair<size_type, Eigen::Vector3d>>(),
           qi::locals<size_type>>
      nodeGroup_;

  /// locals of elementGroup_ are: (# elements, current element type nr, # tags,
  /// # elements read so far)
  qi::rule<ITERATOR, std::vector<GMshFileV2::Element>(),
           qi::locals<size_type, int, int, int, size_type>>
      elementGroup_;

  qi::rule<ITERATOR, std::vector<std::pair<size_type, size_type>>(),
           qi::locals<size_type>, ascii::space_type>
      periodicEntityNodeMapping_;
  qi::rule<ITERATOR, GMshFileV2::PeriodicEntity(), ascii::space_type>
      periodicEntity_;
  qi::rule<ITERATOR, std::vector<GMshFileV2::PeriodicEntity>(),
           qi::locals<size_type>, ascii::space_type>
      periodicEntityGroup_;

  qi::rule<ITERATOR, boost::fusion::adapted::MshFileAdapted(),
           ascii::space_type>
      start_;

  struct ErrorHandler {
    template <class, class, class, class>
    struct result {
      using type = void;
    };

    template <class FIRST, class LAST, class ERROR_POS, class WHAT>
    void operator()(FIRST first, LAST last, ERROR_POS /*errorPos*/,
                    WHAT what) const {
      std::string input(first, last);
      if (input.length() > 40) {
        input = input.substr(0, 40);
      }
      std::cout << "Error in MshFile! Expecting " << what << " here: \""
                << input << "\"" << std::endl;
    }
  };
  phoenix::function<ErrorHandler> errorHandler_;
};

}  // namespace

const std::vector<GMshFileV2::ElementType> GMshFileV2::AllElementTypes{
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

GMshFileV2 readGmshFileV2(std::string::const_iterator begin,
                          std::string::const_iterator end,
                          const std::string& version, bool is_binary,
                          int size_t_size, int one,
                          const std::string& filename) {
  LF_VERIFY_MSG(version == "2.2",
                "Version " << version << " not supported by readGmshFileV2");
  LF_ASSERT_MSG(size_t_size == 8, "Size of std::size_t must be 8.");

  GMshFileV2 result;
  result.IsBinary = is_binary;
  result.VersionNumber = version;
  result.DoubleSize = size_t_size;

  // Parse the rest of the document
  //////////////////////////////////////////////////////////////////////////

  // Setup parsers for node/element sections (which are different depending on
  // binary/non-binary files):
  //
  // Note vec3 has no skipper because it may be used inside lexeme and lexeme
  // can only use parsers without skippers!
  // http://boost-spirit.com/home/2010/02/24/parsing-skippers-and-skipping-parsers/
  // (see comment section)
  using iterator_t = std::string::const_iterator;
  qi::rule<iterator_t, Eigen::Vector3d> vec3;
  qi::rule<iterator_t, std::pair<size_type, Eigen::Vector3d>()> node;
  qi::rule<iterator_t, GMshFileV2::Element(), qi::locals<int>> elementText;
  qi::rule<iterator_t, GMshFileV2::Element(GMshFileV2::ElementType, int, int)>
      elementBin;
  qi::rule<iterator_t, std::vector<GMshFileV2::Element>(),
           qi::locals<size_type, int, int, int, size_type>>
      elementGroup;

  using phoenix::reserve;
  using qi::omit;
  using qi::repeat;
  using qi::labels::_a;
  using qi::labels::_b;
  using qi::labels::_c;
  using qi::labels::_d;
  using qi::labels::_e;
  using qi::labels::_r1;
  using qi::labels::_r2;
  using qi::labels::_r3;
  using qi::labels::_val;

  if (!is_binary) {
    // Text file
    vec3 = qi::double_ >> ' ' >> qi::double_ >> ' ' >> qi::double_;
    node = qi::uint_ >> ' ' >> vec3 >> qi::eol;
    elementText %= qi::int_ > ' ' > qi::int_ > ' ' >
                   qi::omit[qi::int_[qi::_a = qi::_1]] > ' ' > qi::int_ > ' ' >
                   qi::int_ > ' ' >
                   ((qi::eps(_a > 2) >> omit[qi::int_] >> ' ') || qi::eps) >
                   qi::repeat(qi::_a - 3)[qi::int_ >> ' '] > (qi::uint_ % ' ') >
                   *(qi::blank) > qi::eol;
    elementGroup %= "$Elements" > qi::eol >
                    qi::omit[qi::uint_[phoenix::reserve(qi::_val, qi::_1),
                                       qi::_a = qi::_1]] > qi::eol >
                    qi::repeat(qi::_a)[elementText] > "$EndElements";
  } else if (is_binary && one == 1) {
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
        >>
        omit[*((qi::eps(_e < _a) >> qi::little_dword[_b = qi::_1] >>
                qi::little_dword[_c = qi::_1] >>
                qi::little_dword[_d = qi::_1]  // elements-header-binary
                >> repeat(_c)[elementBin(
                       phoenix::static_cast_<GMshFileV2::ElementType>(_b), _d,
                       numNodesAdapted(
                           phoenix::static_cast_<GMshFileV2::ElementType>(
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
        >>
        omit[*((qi::eps(_e < _a) >> qi::big_dword[_b = qi::_1] >>
                qi::big_dword[_c = qi::_1] >>
                qi::big_dword[_d = qi::_1]  // elements-header-binary
                >> repeat(_c)[elementBin(
                       phoenix::static_cast_<GMshFileV2::ElementType>(_b), _d,
                       numNodesAdapted(
                           phoenix::static_cast_<GMshFileV2::ElementType>(
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
  bool r = qi::phrase_parse(begin, end, mshGrammar, ascii::space, result);

  // if (r && iter == end) std::cout << "Parsing succeeded" << std::endl;
  // else if (r) std::cout << "Parsing partially succeeded" << std::endl;
  // std::cout << result << std::endl;

  LF_VERIFY_MSG(r, "Could not parse file " << filename);
  LF_VERIFY_MSG(begin == end, "Could not parse all of file " << filename);

  return result;
}

}  // namespace lf::io
