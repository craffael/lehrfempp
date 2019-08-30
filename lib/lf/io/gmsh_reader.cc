#include "gmsh_reader.h"
#include "eigen_fusion_adapter.h"

#include <lf/geometry/geometry.h>
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

// Structures that represent the MshFile:
namespace lf::io {

/// Output the element type onto the console:
std::ostream& operator<<(std::ostream& stream, MshFile::ElementType et) {
  switch (et) {
    default:
      break;
    case MshFile::ElementType::EDGE2:
      stream << "EDGE2";
      break;
    case MshFile::ElementType::TRIA3:
      stream << "TRIA3";
      break;
    case MshFile::ElementType::QUAD4:
      stream << "QUAD4";
      break;
    case MshFile::ElementType::TET4:
      stream << "TET4";
      break;
    case MshFile::ElementType::HEX8:
      stream << "HEX8";
      break;
    case MshFile::ElementType::PRISM6:
      stream << "PRISM6";
      break;
    case MshFile::ElementType::PYRAMID5:
      stream << "PYRAMID5";
      break;
    case MshFile::ElementType::EDGE3:
      stream << "EDGE3";
      break;
    case MshFile::ElementType::TRIA6:
      stream << "TRIA6";
      break;
    case MshFile::ElementType::QUAD9:
      stream << "QUAD9";
      break;
    case MshFile::ElementType::TET10:
      stream << "TET10";
      break;
    case MshFile::ElementType::HEX27:
      stream << "HEX27";
      break;
    case MshFile::ElementType::PRISM18:
      stream << "PRISM18";
      break;
    case MshFile::ElementType::PYRAMID14:
      stream << "PYRAMID14";
      break;
    case MshFile::ElementType::POINT:
      stream << "POINT";
      break;
    case MshFile::ElementType::QUAD8:
      stream << "QUAD8";
      break;
    case MshFile::ElementType::HEX20:
      stream << "HEX20";
      break;
    case MshFile::ElementType::PRISM15:
      stream << "PRISM15";
      break;
    case MshFile::ElementType::PYRAMID13:
      stream << "PYRAMID13";
      break;
    case MshFile::ElementType::TRIA9:
      stream << "TRIA9";
      break;
    case MshFile::ElementType::TRIA10:
      stream << "TRIA10";
      break;
    case MshFile::ElementType::TRIA12:
      stream << "TRIA12";
      break;
    case MshFile::ElementType::TRIA15:
      stream << "TRIA15";
      break;
    case MshFile::ElementType::TRIA15_5:
      stream << "TRIA15_5";
      break;
    case MshFile::ElementType::TRIA21:
      stream << "TRIA21";
      break;
    case MshFile::ElementType::EDGE4:
      stream << "EDGE4";
      break;
    case MshFile::ElementType::EDGE5:
      stream << "EDGE5";
      break;
    case MshFile::ElementType::EDGE6:
      stream << "EDGE6";
      break;
    case MshFile::ElementType::TET20:
      stream << "TET20";
      break;
    case MshFile::ElementType::TET35:
      stream << "TET35";
      break;
    case MshFile::ElementType::TET56:
      stream << "TET56";
      break;
    case MshFile::ElementType::HEX64:
      stream << "HEX64";
      break;
    case MshFile::ElementType::HEX125:
      stream << "HEX125";
      break;
  }
  return stream;
}

/// For debugging purposes: Write the MshFile into a stream
std::ostream& operator<<(std::ostream& stream, const MshFile& mf) {
  stream << "GMSH FILE: Ver. " << mf.VersionNumber
         << (mf.IsBinary ? "(Binary)" : "(Text)")
         << ", size of double = " << mf.DoubleSize << std::endl;
  stream << "======================================================="
         << std::endl;
  stream << "PHYSICAL ENTITIES (Dimension, Number, Name):" << std::endl;
  for (auto& pe : mf.PhysicalEntities) {
    stream << "  " << pe.Dimension << "\t , " << pe.Number << "\t , " << pe.Name
           << std::endl;
  }
  stream << "NODES (Number, coords)" << std::endl;
  for (auto& n : mf.Nodes) {
    stream << "  " << n.first << "\t , " << n.second.transpose() << std::endl;
  }
  stream << "ELEMENTS (Number, Type, PhysicalEntity Nr, ElementaryEntityNr, "
            "Mesh partitions to which it belongs, Node numbers in it)"
         << std::endl;
  for (auto& e : mf.Elements) {
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
  for (auto& pe : mf.Periodic) {
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
size_type NumNodes(MshFile::ElementType et) {
  switch (et) {
    default:
      break;
    case MshFile::ElementType::EDGE2:
      return 2;
    case MshFile::ElementType::TRIA3:
      return 3;
    case MshFile::ElementType::QUAD4:
      return 4;
    case MshFile::ElementType::TET4:
      return 4;
    case MshFile::ElementType::HEX8:
      return 8;
    case MshFile::ElementType::PRISM6:
      return 6;
    case MshFile::ElementType::PYRAMID5:
      return 5;
    case MshFile::ElementType::EDGE3:
      return 3;
    case MshFile::ElementType::TRIA6:
      return 6;
    case MshFile::ElementType::QUAD9:
      return 9;
    case MshFile::ElementType::TET10:
      return 10;
    case MshFile::ElementType::HEX27:
      return 27;
    case MshFile::ElementType::PRISM18:
      return 18;
    case MshFile::ElementType::PYRAMID14:
      return 14;
    case MshFile::ElementType::POINT:
      return 1;
    case MshFile::ElementType::QUAD8:
      return 8;
    case MshFile::ElementType::HEX20:
      return 20;
    case MshFile::ElementType::PRISM15:
      return 15;
    case MshFile::ElementType::PYRAMID13:
      return 13;
    case MshFile::ElementType::TRIA9:
      return 9;
    case MshFile::ElementType::TRIA10:
      return 10;
    case MshFile::ElementType::TRIA12:
      return 12;
    case MshFile::ElementType::TRIA15:
      return 15;
    case MshFile::ElementType::TRIA15_5:
      return 15;
    case MshFile::ElementType::TRIA21:
      return 21;
    case MshFile::ElementType::EDGE4:
      return 4;
    case MshFile::ElementType::EDGE5:
      return 5;
    case MshFile::ElementType::EDGE6:
      return 6;
    case MshFile::ElementType::TET20:
      return 20;
    case MshFile::ElementType::TET35:
      return 35;
    case MshFile::ElementType::TET56:
      return 56;
    case MshFile::ElementType::HEX64:
      return 64;
    case MshFile::ElementType::HEX125:
      return 125;
  }
  LF_VERIFY_MSG(false, "unknown Gmsh element type");
  // Make compiler happy:
  return 0;
}

base::RefEl RefElOf(MshFile::ElementType et) {
  switch (et) {
    case MshFile::ElementType::POINT:
      return base::RefEl::kPoint();

    case MshFile::ElementType::EDGE2:
    case MshFile::ElementType::EDGE3:
    case MshFile::ElementType::EDGE4:
    case MshFile::ElementType::EDGE5:
    case MshFile::ElementType::EDGE6:
      return base::RefEl::kSegment();

    case MshFile::ElementType::TRIA3:
    case MshFile::ElementType::TRIA6:
    case MshFile::ElementType::TRIA9:
    case MshFile::ElementType::TRIA10:
    case MshFile::ElementType::TRIA12:
    case MshFile::ElementType::TRIA15:
    case MshFile::ElementType::TRIA15_5:
    case MshFile::ElementType::TRIA21:
      return base::RefEl::kTria();

    case MshFile::ElementType::QUAD4:
    case MshFile::ElementType::QUAD8:
    case MshFile::ElementType::QUAD9:
      return base::RefEl::kQuad();

    case MshFile::ElementType::TET4:
    case MshFile::ElementType::HEX8:
    case MshFile::ElementType::PRISM6:
    case MshFile::ElementType::PYRAMID5:
    case MshFile::ElementType::TET10:
    case MshFile::ElementType::HEX27:
    case MshFile::ElementType::PRISM18:
    case MshFile::ElementType::PYRAMID14:
    case MshFile::ElementType::HEX20:
    case MshFile::ElementType::PRISM15:
    case MshFile::ElementType::PYRAMID13:
    case MshFile::ElementType::TET20:
    case MshFile::ElementType::TET35:
    case MshFile::ElementType::TET56:
    case MshFile::ElementType::HEX64:
    case MshFile::ElementType::HEX125:
    default:
      LF_VERIFY_MSG(
          false, "Reference element not supported for GmshElement type " << et);
  }
}

/// Dimension of the GmshElement type
int DimOf(MshFile::ElementType et) {
  switch (et) {
    case MshFile::ElementType::POINT:
      return 0;
    case MshFile::ElementType::EDGE2:
    case MshFile::ElementType::EDGE3:
    case MshFile::ElementType::EDGE4:
    case MshFile::ElementType::EDGE5:
    case MshFile::ElementType::EDGE6:
      return 1;
    case MshFile::ElementType::TRIA3:
    case MshFile::ElementType::QUAD4:
    case MshFile::ElementType::TRIA6:
    case MshFile::ElementType::QUAD9:
    case MshFile::ElementType::QUAD8:
    case MshFile::ElementType::TRIA9:
    case MshFile::ElementType::TRIA10:
    case MshFile::ElementType::TRIA12:
    case MshFile::ElementType::TRIA15:
    case MshFile::ElementType::TRIA15_5:
    case MshFile::ElementType::TRIA21:
      return 2;
    case MshFile::ElementType::TET4:
    case MshFile::ElementType::HEX8:
    case MshFile::ElementType::PRISM6:
    case MshFile::ElementType::PYRAMID5:
    case MshFile::ElementType::TET10:
    case MshFile::ElementType::HEX27:
    case MshFile::ElementType::PRISM18:
    case MshFile::ElementType::PYRAMID14:
    case MshFile::ElementType::HEX20:
    case MshFile::ElementType::PRISM15:
    case MshFile::ElementType::PYRAMID13:
    case MshFile::ElementType::TET20:
    case MshFile::ElementType::TET35:
    case MshFile::ElementType::TET56:
    case MshFile::ElementType::HEX64:
    case MshFile::ElementType::HEX125:
      return 3;
    default:
      LF_VERIFY_MSG(false, "Unknown GmshElement Type.");
  }
  // Make compiler happy:
  return -1;
}
}  // namespace lf::io

// Boost Fusion Adaptions (needed so boost spirit can parse directly into
// MshFile struct)
//////////////////////////////////////////////////////////////////////////
BOOST_FUSION_ADAPT_STRUCT(lf::io::MshFile::PhysicalEntity,
                          (int, Dimension)(int, Number)(std::string, Name));

BOOST_FUSION_ADAPT_STRUCT(
    lf::io::MshFile::Element,
    (size_type, Number)(lf::io::MshFile::ElementType,
                        Type)(int, PhysicalEntityNr)(int, ElementaryEntityNr)(
        std::vector<int>, MeshPartitions)(std::vector<size_type>, NodeNumbers));

/// To circumvent comma in preprocessor invocation
using nodeMapping_t = std::pair<size_type, size_type>;

BOOST_FUSION_ADAPT_STRUCT(lf::io::MshFile::PeriodicEntity,
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
    lf::io::MshFile, MshFileAdapted,
    //(double, VersionNumber)
    //(bool, IsBinary)
    //(int, DoubleSize)
    (std::vector<lf::io::MshFile::PhysicalEntity>,
     PhysicalEntities)(std::vector<nodePair_t>,
                       Nodes)(std::vector<lf::io::MshFile::Element>, Elements)(
        std::vector<lf::io::MshFile::PeriodicEntity>, Periodic));

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

namespace lf::io {
namespace /*Anonymous*/ {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

/// A lookup table for boost spirit that can parse an element type
struct gmshElementType : qi::symbols<char, unsigned> {
  gmshElementType() {
    for (auto& et : MshFile::AllElementTypes) {
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
      qi::rule<ITERATOR, std::vector<MshFile::Element>(),
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

  qi::rule<ITERATOR, MshFile::PhysicalEntity(), ascii::space_type>
      physicalEntity_;
  qi::rule<ITERATOR, std::vector<MshFile::PhysicalEntity>(),
           qi::locals<size_type>, ascii::space_type>
      physicalEntityGroup_;

  qi::rule<ITERATOR, std::pair<size_type, Eigen::Vector3d>()> node_;
  qi::rule<ITERATOR, std::vector<std::pair<size_type, Eigen::Vector3d>>(),
           qi::locals<size_type>>
      nodeGroup_;

  /// locals of elementGroup_ are: (# elements, current element type nr, # tags,
  /// # elements read so far)
  qi::rule<ITERATOR, std::vector<MshFile::Element>(),
           qi::locals<size_type, int, int, int, size_type>>
      elementGroup_;

  qi::rule<ITERATOR, std::vector<std::pair<size_type, size_type>>(),
           qi::locals<size_type>, ascii::space_type>
      periodicEntityNodeMapping_;
  qi::rule<ITERATOR, MshFile::PeriodicEntity(), ascii::space_type>
      periodicEntity_;
  qi::rule<ITERATOR, std::vector<MshFile::PeriodicEntity>(),
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

const std::vector<MshFile::ElementType> MshFile::AllElementTypes{
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

MshFile readGMshFile(const std::string& filename) {
  // Open file and copy into memory:hydi::io::MshFile
  //////////////////////////////////////////////////////////////////////////
  std::ifstream in(filename, std::ios_base::in);
  if (!in) {
    std::string error("Could not open file ");
    error += filename;
    throw base::LfException(error);
  }
  std::string storage;
  in.unsetf(std::ios::skipws);  // No white space skipping
  std::copy(std::istream_iterator<char>(in), std::istream_iterator<char>(),
            std::back_inserter(storage));

  // Parse header to determine if we are dealing with ASCII format or binary
  // format + little or big endian:
  //////////////////////////////////////////////////////////////////////////
  MshFile result;
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
  qi::rule<iterator_t, MshFile::Element(), qi::locals<int>> elementText;
  qi::rule<iterator_t, MshFile::Element(MshFile::ElementType, int, int)>
      elementBin;
  qi::rule<iterator_t, std::vector<MshFile::Element>(),
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

  if (!result.IsBinary) {
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
        >>
        omit[*((qi::eps(_e < _a) >> qi::little_dword[_b = qi::_1] >>
                qi::little_dword[_c = qi::_1] >>
                qi::little_dword[_d = qi::_1]  // elements-header-binary
                >>
                repeat(_c)[elementBin(
                    phoenix::static_cast_<MshFile::ElementType>(_b), _d,
                    numNodesAdapted(phoenix::static_cast_<MshFile::ElementType>(
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
                >>
                repeat(_c)[elementBin(
                    phoenix::static_cast_<MshFile::ElementType>(_b), _d,
                    numNodesAdapted(phoenix::static_cast_<MshFile::ElementType>(
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

bool GmshReader::IsPhysicalEntity(const mesh::Entity& e,
                                  size_type physical_entity_nr) const {
  auto physical_entities = PhysicalEntityNr(e);
  return std::find(physical_entities.begin(), physical_entities.end(),
                   physical_entity_nr) != physical_entities.end();
}

GmshReader::GmshReader(std::unique_ptr<mesh::MeshFactory> factory,
                       const ::lf::io::MshFile& msh_file)
    : mesh_factory_(std::move(factory)) {
  // 1) Check Gmsh_file and initialize
  //////////////////////////////////////////////////////////////////////////

  dim_t dim_mesh = mesh_factory_->DimMesh();
  dim_t dim_world = mesh_factory_->DimWorld();
  LF_VERIFY_MSG(
      dim_mesh >= 2 && dim_mesh <= 3 && dim_world >= 2 && dim_world <= 3,
      "GmshReader supports only 2D and 3D meshes.");

  // 1) Determine which nodes of gmsh are also nodes of the LehrFEM++ mesh +
  // count number of entities of each codimension (exclude auxilliary nodes).
  // This is necessary because e.g. second order meshes in gmsh need a lot of
  // auxilliary nodes to define their geometry...
  /////////////////////////////////////////////////////////////////////////////
  // is_main_node[i] = true means that gmsh-node i is one of the main nodes of
  // an element
  std::vector<bool> is_main_node(msh_file.Nodes.size(), false);

  // mi2gi = mesh_index_2_gmsh_index
  // mi2gi[c][i] contains the gmsh entities that belong to the mesh entity with
  //             codim = c and mesh index i.
  std::vector<std::vector<std::vector<size_type>>> mi2gi(dim_mesh + 1);

  {
    // count the number of entities for each codimension and reserve space:
    std::vector<size_type> num_entities(mesh_factory_->DimMesh() + 1, 0);

    for (const auto& e : msh_file.Elements) {
      LF_ASSERT_MSG(DimOf(e.Type) <= dim_mesh,
                    "mesh_factory->DimMesh() = "
                        << dim_mesh
                        << ", but msh-file contains entities with dimension "
                        << DimOf(e.Type));

      ++num_entities[DimOf(e.Type)];

      if (DimOf(e.Type) == dim_mesh) {
        // mark main nodes
        auto ref_el = RefElOf(e.Type);
        for (unsigned int i = 0; i < ref_el.NumNodes(); ++i) {
          auto node_number = e.NodeNumbers[i];
          if (is_main_node.size() <= node_number) {
            is_main_node.resize(node_number + 1);
          }
          is_main_node[e.NodeNumbers[i]] = true;
        }
      }
    }

    for (dim_t c = 0; c <= dim_mesh; ++c) {
      mi2gi[c].reserve(num_entities[dim_mesh - c]);
    }

    LF_ASSERT_MSG(num_entities[dim_mesh] > 0,
                  "MshFile contains no elements with dimension " << dim_mesh);
  }

  // 2) Insert main nodes into MeshFactory
  /////////////////////////////////////////////////////////////////////////////

  // gmsh_index_2_mesh_index for nodes:
  // gi2mi[i] = j means that gmsh node with gmsh index i has mesh index j
  std::vector<size_type> gi2mi;
  gi2mi.resize(is_main_node.size(), -1);

  // gi2i[i] = j implies msh_file.Nodes[j].first = i
  std::vector<size_type> gi2i;
  gi2i.resize(is_main_node.size(), -1);

  for (std::size_t i = 0; i < msh_file.Nodes.size(); ++i) {
    auto& n = msh_file.Nodes[i];
    if (gi2i.size() <= n.first) {
      gi2i.resize(n.first + 1, -1);
    }
    gi2i[n.first] = i;

    if (is_main_node.size() <= n.first || !is_main_node[n.first]) {
      continue;
    }

    size_type mi;
    if (dim_world == 2) {
      LF_ASSERT_MSG(
          n.second(2) == 0,
          "In a 2D GmshMesh, the z-coordinate of every node must be zero");
      mi = mesh_factory_->AddPoint(n.second.topRows(2));
    } else {
      mi = mesh_factory_->AddPoint(n.second);
    }
    gi2mi[n.first] = mi;
  }

  // 3) Insert entities (except nodes) into MeshFactory:
  //////////////////////////////////////////////////////////////////////////////

  std::size_t begin = 0;
  for (std::size_t end = 0; end < msh_file.Elements.size(); ++end) {
    auto& begin_element = msh_file.Elements[begin];
    auto& end_element = msh_file.Elements[end];
    auto ref_el = RefElOf(end_element.Type);
    auto codim = dim_mesh - ref_el.Dimension();
    if (begin_element.NodeNumbers == end_element.NodeNumbers && begin != end &&
        begin_element.Type == end_element.Type) {
      // This entity appears more than once
      mi2gi[codim].back().push_back(end);
      continue;
    }

    begin = end;
    if (ref_el == base::RefEl::kPoint()) {
      // special case, this entity is a point (which has already been inserted)
      auto mesh_index = gi2mi[end_element.NodeNumbers[0]];
      if (mi2gi[dim_mesh].size() <= mesh_index) {
        mi2gi[dim_mesh].resize(mesh_index + 1);
      }
      mi2gi[dim_mesh][mesh_index].push_back(end);
    } else {
      // gmsh element is not a point -> insert entity:
      auto num_nodes = end_element.NodeNumbers.size();
      Eigen::MatrixXd node_coords(dim_world, num_nodes);
      for (std::size_t i = 0; i < num_nodes; ++i) {
        auto node_coord =
            msh_file.Nodes[gi2i[end_element.NodeNumbers[i]]].second;
        if (dim_world == 2) {
          node_coords.col(i) = node_coord.topRows(2);
        } else {
          node_coords.col(i) = node_coord;
        }
      }
      std::unique_ptr<geometry::Geometry> geom;

      switch (end_element.Type) {
        case MshFile::ElementType::EDGE2:
          ref_el = base::RefEl::kSegment();
          geom = std::make_unique<geometry::SegmentO1>(node_coords);
          break;
        case MshFile::ElementType::EDGE3:
          ref_el = base::RefEl::kSegment();
          geom = std::make_unique<geometry::SegmentO2>(node_coords);
          break;
        case MshFile::ElementType::TRIA3:
          ref_el = base::RefEl::kTria();
          geom = std::make_unique<geometry::TriaO1>(node_coords);
          break;
        case MshFile::ElementType::TRIA6:
          ref_el = base::RefEl::kTria();
          geom = std::make_unique<geometry::TriaO2>(node_coords);
          break;
        case MshFile::ElementType::QUAD4:
          ref_el = base::RefEl::kQuad();
          geom = std::make_unique<geometry::QuadO1>(node_coords);
          break;
        case MshFile::ElementType::QUAD8:
          ref_el = base::RefEl::kQuad();
          geom = std::make_unique<geometry::QuadO2>(node_coords);
          break;
        case MshFile::ElementType::QUAD9:
          ref_el = base::RefEl::kQuad();
          geom = std::make_unique<geometry::QuadO2>(node_coords.leftCols(8));
          break;
        default:
          LF_VERIFY_MSG(false, "Gmsh element type "
                                   << end_element.Type
                                   << " not (yet) supported by GmshReader.");
      }
      std::vector<size_type> main_nodes(ref_el.NumNodes());
      for (dim_t i = 0; i < ref_el.NumNodes(); ++i) {
        main_nodes[i] = gi2mi[end_element.NodeNumbers[i]];
      }

      mesh_factory_->AddEntity(ref_el, main_nodes, std::move(geom));
      mi2gi[codim].emplace_back(std::vector{static_cast<unsigned int>(end)});
    }
  }

  // 4) Construct mesh
  //////////////////////////////////////////////////////////////////////////////
  mesh_ = mesh_factory_->Build();

  // 5) Build MeshDataSet that assigns the physical entitiies:
  //////////////////////////////////////////////////////////////////////////////
  physical_nrs_ =
      mesh::utils::make_AllCodimMeshDataSet<std::vector<size_type>>(mesh_);

  for (dim_t c = 0; c <= dim_mesh; ++c) {
    for (auto& e : mesh_->Entities(c)) {
      auto mi = mesh_->Index(e);
      if (c == dim_mesh && mi >= mi2gi[dim_mesh].size()) {
        // this point did not appear as a gmsh element in the file -> don't
        // assign any physical entity nr.
        continue;
      }
      if (mi2gi[c].size() > mi) {
        std::vector<size_type> temp;
        for (auto& gmsh_index : mi2gi[c][mi]) {
          temp.push_back(msh_file.Elements[gmsh_index].PhysicalEntityNr);
        }

        physical_nrs_->operator()(e) = std::move(temp);
      }
    }
  }

  // 6) Create mapping physicalEntityNr <-> physicalEntityName:
  //////////////////////////////////////////////////////////////////////////

  for (auto& pe : msh_file.PhysicalEntities) {
    name_2_nr_.insert(
        std::pair{pe.Name, std::pair{pe.Number, dim_mesh - pe.Dimension}});
    nr_2_name_.insert(
        std::pair{pe.Number, std::pair{pe.Name, dim_mesh - pe.Dimension}});
  }

  if (!msh_file.Periodic.empty()) {
    /*LOGGER_ENTRY(logger_,
                 "WARNING: GMSH File  contains periodic boundary relations "
                 "between elements. These are ignored by GmshReader.",
                 3);*/
  }
}

GmshReader::GmshReader(std::unique_ptr<mesh::MeshFactory> factory,
                       const std::string& filename)
    : GmshReader(std::move(factory), readGMshFile(filename)) {}

size_type GmshReader::PhysicalEntityName2Nr(const std::string& name,
                                            dim_t codim) const {
  LF_ASSERT_MSG(!name.empty(), "name is empty");
  auto [begin, end] = name_2_nr_.equal_range(name);  // NOLINT
  if (begin == end) {
    throw base::LfException("No Physical Entity with this name found.");
  }
  auto result = *begin;
  ++begin;
  if (begin == end) {
    if (codim == dim_t(-1) || codim == result.second.second) {
      return result.second.first;
    }
  } else {
    if (codim == dim_t(-1)) {
      throw base::LfException(
          "There are multiple physical entities with the name " + name +
          ", please specify also the codimension.");
    }
    if (result.second.second == codim) {
      return result.second.first;
    }
    while (begin->second.second != codim && begin != end) {
      ++begin;
    }
    if (begin->second.second == codim) {
      return begin->second.first;
    }
  }
  throw base::LfException("Physical Entity with name='" + name +
                          "' and codimension=" + std::to_string(codim) +
                          "' not found.");
}

std::string GmshReader::PhysicalEntityNr2Name(size_type number,
                                              dim_t codim) const {
  auto [begin, end] = nr_2_name_.equal_range(number);  // NOLINT
  if (begin == end) {
    throw base::LfException("Physical entity with number " +
                            std::to_string(number) + " not found.");
  }
  auto result = *begin;
  ++begin;
  if (begin == end) {
    if (codim == dim_t(-1) || result.second.second == codim) {
      return result.second.first;
    }
  } else {
    if (codim == dim_t(-1)) {
      throw base::LfException(
          "There are multiple physical entities with the Number " +
          std::to_string(number) + ", please specify also the codimension");
    }
    if (result.second.second == codim) {
      return result.second.first;
    }
    while (begin->second.second != codim && begin != end) {
      ++begin;
    }
    if (begin->second.second == codim) {
      return begin->second.first;
    }
  }
  throw base::LfException(
      "Physical entity with number=" + std::to_string(number) +
      ", codim=" + std::to_string(codim) + " not found.");
}

std::vector<std::pair<size_type, std::string>> GmshReader::PhysicalEntities(
    dim_t codim) const {
  std::vector<std::pair<size_type, std::string>> result;
  for (auto& p : nr_2_name_) {
    if (p.second.second != codim) {
      continue;
    }
    result.emplace_back(p.first, p.second.first);
  }
  return result;
}

std::vector<size_type> GmshReader::PhysicalEntityNr(
    const mesh::Entity& e) const {
  return physical_nrs_->operator()(e);
}

}  // namespace lf::io
