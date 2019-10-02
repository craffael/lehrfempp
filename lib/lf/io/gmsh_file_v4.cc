/**
 * @file
 * @brief Implementation of all logic related to the class GmshFileV4
 * @author Raffael Casagrande
 * @date   2019-09-02 03:02:19
 * @copyright MIT License
 */

#include "gmsh_file_v4.h"

#include <boost/fusion/adapted.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/adapt_struct_named.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/include/boost_array.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/fusion/iterator.hpp>
#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/support/category_of.hpp>
#include <boost/fusion/support/iterator_base.hpp>
#include <boost/fusion/support/tag_of.hpp>
#include <boost/fusion/support/tag_of_fwd.hpp>
#include <boost/mpl/minus.hpp>
#include <boost/phoenix/function/adapt_function.hpp>
#include <boost/phoenix/fusion.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_binary.hpp>
#include <boost/spirit/include/qi_uint.hpp>
#include "eigen_fusion_adapter.h"

namespace lf::io {

/// Output the element type onto a stream:
std::ostream& operator<<(std::ostream& stream, GMshFileV4::ElementType et) {
  switch (et) {
    default:
      break;
    case GMshFileV4::ElementType::EDGE2:
      stream << "EDGE2";
      break;
    case GMshFileV4::ElementType::TRIA3:
      stream << "TRIA3";
      break;
    case GMshFileV4::ElementType::QUAD4:
      stream << "QUAD4";
      break;
    case GMshFileV4::ElementType::TET4:
      stream << "TET4";
      break;
    case GMshFileV4::ElementType::HEX8:
      stream << "HEX8";
      break;
    case GMshFileV4::ElementType::PRISM6:
      stream << "PRISM6";
      break;
    case GMshFileV4::ElementType::PYRAMID5:
      stream << "PYRAMID5";
      break;
    case GMshFileV4::ElementType::EDGE3:
      stream << "EDGE3";
      break;
    case GMshFileV4::ElementType::TRIA6:
      stream << "TRIA6";
      break;
    case GMshFileV4::ElementType::QUAD9:
      stream << "QUAD9";
      break;
    case GMshFileV4::ElementType::TET10:
      stream << "TET10";
      break;
    case GMshFileV4::ElementType::HEX27:
      stream << "HEX27";
      break;
    case GMshFileV4::ElementType::PRISM18:
      stream << "PRISM18";
      break;
    case GMshFileV4::ElementType::PYRAMID14:
      stream << "PYRAMID14";
      break;
    case GMshFileV4::ElementType::POINT:
      stream << "POINT";
      break;
    case GMshFileV4::ElementType::QUAD8:
      stream << "QUAD8";
      break;
    case GMshFileV4::ElementType::HEX20:
      stream << "HEX20";
      break;
    case GMshFileV4::ElementType::PRISM15:
      stream << "PRISM15";
      break;
    case GMshFileV4::ElementType::PYRAMID13:
      stream << "PYRAMID13";
      break;
    case GMshFileV4::ElementType::TRIA9:
      stream << "TRIA9";
      break;
    case GMshFileV4::ElementType::TRIA10:
      stream << "TRIA10";
      break;
    case GMshFileV4::ElementType::TRIA12:
      stream << "TRIA12";
      break;
    case GMshFileV4::ElementType::TRIA15:
      stream << "TRIA15";
      break;
    case GMshFileV4::ElementType::TRIA15_5:
      stream << "TRIA15_5";
      break;
    case GMshFileV4::ElementType::TRIA21:
      stream << "TRIA21";
      break;
    case GMshFileV4::ElementType::EDGE4:
      stream << "EDGE4";
      break;
    case GMshFileV4::ElementType::EDGE5:
      stream << "EDGE5";
      break;
    case GMshFileV4::ElementType::EDGE6:
      stream << "EDGE6";
      break;
    case GMshFileV4::ElementType::TET20:
      stream << "TET20";
      break;
    case GMshFileV4::ElementType::TET35:
      stream << "TET35";
      break;
    case GMshFileV4::ElementType::TET56:
      stream << "TET56";
      break;
    case GMshFileV4::ElementType::HEX64:
      stream << "HEX64";
      break;
    case GMshFileV4::ElementType::HEX125:
      stream << "HEX125";
      break;
  }
  return stream;
}

/// Number of nodes that this element type has
int NumNodes(GMshFileV4::ElementType et) {
  switch (et) {
    default:
      break;
    case GMshFileV4::ElementType::EDGE2:
      return 2;
    case GMshFileV4::ElementType::TRIA3:
      return 3;
    case GMshFileV4::ElementType::QUAD4:
      return 4;
    case GMshFileV4::ElementType::TET4:
      return 4;
    case GMshFileV4::ElementType::HEX8:
      return 8;
    case GMshFileV4::ElementType::PRISM6:
      return 6;
    case GMshFileV4::ElementType::PYRAMID5:
      return 5;
    case GMshFileV4::ElementType::EDGE3:
      return 3;
    case GMshFileV4::ElementType::TRIA6:
      return 6;
    case GMshFileV4::ElementType::QUAD9:
      return 9;
    case GMshFileV4::ElementType::TET10:
      return 10;
    case GMshFileV4::ElementType::HEX27:
      return 27;
    case GMshFileV4::ElementType::PRISM18:
      return 18;
    case GMshFileV4::ElementType::PYRAMID14:
      return 14;
    case GMshFileV4::ElementType::POINT:
      return 1;
    case GMshFileV4::ElementType::QUAD8:
      return 8;
    case GMshFileV4::ElementType::HEX20:
      return 20;
    case GMshFileV4::ElementType::PRISM15:
      return 15;
    case GMshFileV4::ElementType::PYRAMID13:
      return 13;
    case GMshFileV4::ElementType::TRIA9:
      return 9;
    case GMshFileV4::ElementType::TRIA10:
      return 10;
    case GMshFileV4::ElementType::TRIA12:
      return 12;
    case GMshFileV4::ElementType::TRIA15:
      return 15;
    case GMshFileV4::ElementType::TRIA15_5:
      return 15;
    case GMshFileV4::ElementType::TRIA21:
      return 21;
    case GMshFileV4::ElementType::EDGE4:
      return 4;
    case GMshFileV4::ElementType::EDGE5:
      return 5;
    case GMshFileV4::ElementType::EDGE6:
      return 6;
    case GMshFileV4::ElementType::TET20:
      return 20;
    case GMshFileV4::ElementType::TET35:
      return 35;
    case GMshFileV4::ElementType::TET56:
      return 56;
    case GMshFileV4::ElementType::HEX64:
      return 64;
    case GMshFileV4::ElementType::HEX125:
      return 125;
  }
  LF_VERIFY_MSG(false, "unknown Gmsh element type");
  // Make compiler happy:
  return 0;
}

base::RefEl RefElOf(GMshFileV4::ElementType et) {
  switch (et) {
    case GMshFileV4::ElementType::POINT:
      return base::RefEl::kPoint();

    case GMshFileV4::ElementType::EDGE2:
    case GMshFileV4::ElementType::EDGE3:
    case GMshFileV4::ElementType::EDGE4:
    case GMshFileV4::ElementType::EDGE5:
    case GMshFileV4::ElementType::EDGE6:
      return base::RefEl::kSegment();

    case GMshFileV4::ElementType::TRIA3:
    case GMshFileV4::ElementType::TRIA6:
    case GMshFileV4::ElementType::TRIA9:
    case GMshFileV4::ElementType::TRIA10:
    case GMshFileV4::ElementType::TRIA12:
    case GMshFileV4::ElementType::TRIA15:
    case GMshFileV4::ElementType::TRIA15_5:
    case GMshFileV4::ElementType::TRIA21:
      return base::RefEl::kTria();

    case GMshFileV4::ElementType::QUAD4:
    case GMshFileV4::ElementType::QUAD8:
    case GMshFileV4::ElementType::QUAD9:
      return base::RefEl::kQuad();

    case GMshFileV4::ElementType::TET4:
    case GMshFileV4::ElementType::HEX8:
    case GMshFileV4::ElementType::PRISM6:
    case GMshFileV4::ElementType::PYRAMID5:
    case GMshFileV4::ElementType::TET10:
    case GMshFileV4::ElementType::HEX27:
    case GMshFileV4::ElementType::PRISM18:
    case GMshFileV4::ElementType::PYRAMID14:
    case GMshFileV4::ElementType::HEX20:
    case GMshFileV4::ElementType::PRISM15:
    case GMshFileV4::ElementType::PYRAMID13:
    case GMshFileV4::ElementType::TET20:
    case GMshFileV4::ElementType::TET35:
    case GMshFileV4::ElementType::TET56:
    case GMshFileV4::ElementType::HEX64:
    case GMshFileV4::ElementType::HEX125:
    default:
      LF_VERIFY_MSG(
          false, "Reference element not supported for GmshElement type " << et);
  }
}

/// Dimension of the GmshElement type
int DimOf(GMshFileV4::ElementType et) {
  switch (et) {
    case GMshFileV4::ElementType::POINT:
      return 0;
    case GMshFileV4::ElementType::EDGE2:
    case GMshFileV4::ElementType::EDGE3:
    case GMshFileV4::ElementType::EDGE4:
    case GMshFileV4::ElementType::EDGE5:
    case GMshFileV4::ElementType::EDGE6:
      return 1;
    case GMshFileV4::ElementType::TRIA3:
    case GMshFileV4::ElementType::QUAD4:
    case GMshFileV4::ElementType::TRIA6:
    case GMshFileV4::ElementType::QUAD9:
    case GMshFileV4::ElementType::QUAD8:
    case GMshFileV4::ElementType::TRIA9:
    case GMshFileV4::ElementType::TRIA10:
    case GMshFileV4::ElementType::TRIA12:
    case GMshFileV4::ElementType::TRIA15:
    case GMshFileV4::ElementType::TRIA15_5:
    case GMshFileV4::ElementType::TRIA21:
      return 2;
    case GMshFileV4::ElementType::TET4:
    case GMshFileV4::ElementType::HEX8:
    case GMshFileV4::ElementType::PRISM6:
    case GMshFileV4::ElementType::PYRAMID5:
    case GMshFileV4::ElementType::TET10:
    case GMshFileV4::ElementType::HEX27:
    case GMshFileV4::ElementType::PRISM18:
    case GMshFileV4::ElementType::PYRAMID14:
    case GMshFileV4::ElementType::HEX20:
    case GMshFileV4::ElementType::PRISM15:
    case GMshFileV4::ElementType::PYRAMID13:
    case GMshFileV4::ElementType::TET20:
    case GMshFileV4::ElementType::TET35:
    case GMshFileV4::ElementType::TET56:
    case GMshFileV4::ElementType::HEX64:
    case GMshFileV4::ElementType::HEX125:
      return 3;
    default:
      LF_VERIFY_MSG(false, "Unknown GmshElement Type.");
  }
  // Make compiler happy:
  return -1;
}
}  // namespace lf::io

// Boost Fusion Adaptions (needed so boost spirit can parse directly into
// GmshFileV4 struct)
///////////////////////////////////////////////////////////////////////////////
///

BOOST_FUSION_ADAPT_STRUCT(lf::io::GMshFileV4::PhysicalName,
                          (int, dimension)(int, physical_tag)(std::string,
                                                              name));

BOOST_FUSION_ADAPT_STRUCT(lf::io::GMshFileV4::PointEntity,
                          (int, tag)(Eigen::Vector3d, coord)(std::vector<int>,
                                                             physical_tags));

BOOST_FUSION_ADAPT_STRUCT(
    lf::io::GMshFileV4::Entity,
    (int, tag)(Eigen::Vector3d, min_coord)(Eigen::Vector3d, max_coord)(
        std::vector<int>, physical_tags)(std::vector<int>, bounding_entities));

BOOST_FUSION_ADAPT_STRUCT(lf::io::GMshFileV4::GhostEntity,
                          (int, tag)(int, partition));

BOOST_FUSION_ADAPT_STRUCT(
    lf::io::GMshFileV4::PartitionedPointEntity,
    (int, tag)(int, parent_dim)(int, parent_tag)(std::vector<int>, partitions)(
        Eigen::Vector3d, coord)(std::vector<int>, physical_tags));

BOOST_FUSION_ADAPT_STRUCT(
    lf::io::GMshFileV4::PartitionedEntity,
    (int, tag)(int, parent_dim)(int, parent_tag)(std::vector<int>, partitions)(
        Eigen::Vector3d, min_coord)(Eigen::Vector3d, max_coord)(
        std::vector<int>, physical_tags)(std::vector<int>, bounding_entities));

using partitionedEntities_t =
    std::tuple<std::vector<lf::io::GMshFileV4::PartitionedPointEntity>,
               std::vector<lf::io::GMshFileV4::PartitionedEntity>,
               std::vector<lf::io::GMshFileV4::PartitionedEntity>,
               std::vector<lf::io::GMshFileV4::PartitionedEntity>>;

BOOST_FUSION_ADAPT_STRUCT(
    lf::io::GMshFileV4::PartitionedEntities,
    (std::size_t, num_partitions)(std::vector<lf::io::GMshFileV4::GhostEntity>,
                                  ghost_entities)(partitionedEntities_t,
                                                  partitioned_entities));

// To prevent comma in preprocessor invocation
using nodeMapping_t = std::pair<std::size_t, Eigen::Vector3d>;

BOOST_FUSION_ADAPT_STRUCT(lf::io::GMshFileV4::NodeBlock,
                          (int, entity_dim)(int, entity_tag)(bool, parametric)(
                              std::vector<nodeMapping_t>, nodes));

BOOST_FUSION_ADAPT_STRUCT(
    lf::io::GMshFileV4::Nodes,
    (int, num_nodes)(std::size_t, min_node_tag)(std::size_t, max_node_tag)(
        std::vector<lf::io::GMshFileV4::NodeBlock>, node_blocks));

// To prevent comma in preprocessor invocation
using elementMapping_t = std::tuple<std::size_t, std::vector<std::size_t>>;

BOOST_FUSION_ADAPT_STRUCT(lf::io::GMshFileV4::ElementBlock,
                          (int, dimension)(int, entity_tag)(
                              lf::io::GMshFileV4::ElementType,
                              element_type)(std::vector<elementMapping_t>,
                                            elements));

BOOST_FUSION_ADAPT_STRUCT(
    lf::io::GMshFileV4::Elements,
    (std::size_t, num_elements)(std::size_t, min_element_tag)(std::size_t,
                                                              max_element_tag)(
        std::vector<lf::io::GMshFileV4::ElementBlock>, element_blocks));

// to prevent comma in preprocessor invocation
using periodicLink_t = std::pair<std::size_t, std::size_t>;
BOOST_FUSION_ADAPT_STRUCT(
    lf::io::GMshFileV4::PeriodicLink,
    (int, dimension)(int, entity_tag_slave)(int, entity_tag_master)(
        std::optional<Eigen::Matrix4d>,
        affine_transform)(std::vector<periodicLink_t>, node_mapping));

BOOST_FUSION_ADAPT_STRUCT(lf::io::GMshFileV4::GhostElement,
                          (std::size_t, element_tag)(int, partition_tag)(
                              std::vector<int>, ghost_partition_tags));

// to prevent comma in preprocessor invocation:
using entities_t = std::tuple<std::vector<lf::io::GMshFileV4::PointEntity>,
                              std::vector<lf::io::GMshFileV4::Entity>,
                              std::vector<lf::io::GMshFileV4::Entity>,
                              std::vector<lf::io::GMshFileV4::Entity>>;

// clang-format off
BOOST_FUSION_ADAPT_STRUCT_NAMED(
    lf::io::GMshFileV4, MshFileV4Adapted,
    (std::vector<lf::io::GMshFileV4::PhysicalName>, physical_names)
    (entities_t, entities)
    (lf::io::GMshFileV4::PartitionedEntities, partitioned_entities)
    (lf::io::GMshFileV4::Nodes, nodes)
    (lf::io::GMshFileV4::Elements, elements)
    (std::vector<lf::io::GMshFileV4::PeriodicLink>, periodic_links)
    (std::vector<lf::io::GMshFileV4::GhostElement>, ghost_elements)
)
// clang-format on

namespace boost::spirit::traits {

namespace /*anonymous*/ {
template <class T>
constexpr bool has_value_type(long) {
  return false;
}

template <class T, typename = typename T::value_type>
constexpr bool has_value_type(int) {
  return true;
}

}  // namespace

template <class Enum, class RawValue>
struct assign_to_attribute_from_value<
    Enum, RawValue,
    typename std::enable_if_t<std::is_enum_v<Enum> &&
                              !std::is_same_v<Enum, RawValue>>> {
  static void call(RawValue const& raw, Enum& cat) {
    if constexpr (has_value_type<RawValue>(0)) {
      // specialization for endian::endian
      typename RawValue::value_type value = raw;
      cat = static_cast<Enum>(value);
    } else {
      cat = static_cast<Enum>(raw);
    }
  }
};

}  // namespace boost::spirit::traits

namespace lf::io {
namespace /* anonymous */ {
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

// A lookup table for boost spirit that can parse an element type
struct GmshElementType : qi::symbols<char, unsigned> {
  GmshElementType() {
    for (auto& et : GMshFileV4::AllElementTypes) {
      add(std::to_string(static_cast<int>(et)), static_cast<int>(et));
    }
  }
};

// NOLINTNEXTLINE
BOOST_PHOENIX_ADAPT_FUNCTION(int, numNodesAdapted, NumNodes, 1);

template <class ITERATOR>
struct MshV4GrammarText
    : qi::grammar<ITERATOR, boost::fusion::adapted::MshFileV4Adapted(),
                  ascii::space_type> {
  MshV4GrammarText() : MshV4GrammarText::base_type(entry_, "Msh File") {
    using phoenix::at_c;
    using phoenix::push_back;
    using phoenix::reserve;
    using phoenix::val;
    using qi::_val;
    using qi::attr;
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
    using qi::labels::_b;
    using qi::labels::_c;
    using qi::labels::_d;

    auto size_t_ = qi::ulong_long;

    // General Parsers:
    quoted_string_ %= lexeme['"' >> +(char_ - '"') >> '"'];
    quoted_string_.name("string");
    vec3_ %= double_ > double_ > double_;
    int_vec_ %= omit[size_t_[reserve(_val, _1), _a = _1]] > repeat(_a)[int_];
    start_comment_ %= !lit("$PhysicalNames") >> !lit("$Entities") >>
                      !lit("$PartitionedEntities") >> !lit("$Nodes") >>
                      !lit("$Elements") >> !lit("$Periodic") >>
                      !lit("$GhostElements") >>
                      (lit('$') >> (+(char_ - qi::eol)));
    start_comment_.name("Start of Comment");
    comment_ %= start_comment_[_a = qi::_1] > *(char_ - '$') >> "$End" >>
                qi::string(_a);
    comment_.name("comment");
    qi::on_error<qi::fail>(comment_, error_handler_(_1, _2, _3, _4));

    // Physical names
    physical_name_ %= int_ > int_ > quoted_string_;
    physical_name_.name("Physical Name");
    qi::on_error<qi::fail>(physical_name_, error_handler_(_1, _2, _3, _4));
    physical_name_vector_ %= "$PhysicalNames" >
                             omit[size_t_[reserve(_val, _1), _a = _1]] >
                             repeat(_a)[physical_name_] > "$EndPhysicalNames";
    physical_name_vector_.name("$PhyiscalNames");
    /*qi::on_error<qi::fail>(physical_name_vector_,
                           error_handler_(_1, _2, _3, _4));*/

    // entities:
    point_entity_ %= int_ > vec3_ > int_vec_;
    point_entity_.name("Point entity");

    entity_ %= int_ > vec3_ > vec3_ > int_vec_ > int_vec_;
    entity_.name("entity");

    entities_ %= "$Entities" >
                 omit[size_t_[reserve(phoenix::at_c<0>(_val), _1), _a = _1]] >
                 omit[size_t_[reserve(phoenix::at_c<1>(_val), _1), _b = _1]] >
                 omit[size_t_[reserve(phoenix::at_c<2>(_val), _1), _c = _1]] >
                 omit[size_t_[reserve(phoenix::at_c<3>(_val), _1), _d = _1]] >
                 repeat(_a)[point_entity_] > repeat(_b)[entity_] >
                 repeat(_c)[entity_] > repeat(_d)[entity_] > "$EndEntities";
    entities_.name("$Entities");
    qi::on_error<qi::fail>(entities_, error_handler_(_1, _2, _3, _4));

    // PartitionedEntities
    ghost_entities_ %=
        omit[size_t_[reserve(_val, _1), _a = _1]] > repeat(_a)[int_ > int_];
    ghost_entities_.name("ghost_entities");
    partitioned_point_entity_ %=
        int_ > int_ > int_ > int_vec_ > vec3_ > int_vec_;
    partitioned_point_entity_.name("partitioned_point_entity");

    partitioned_entity_ %=
        int_ > int_ > int_ > int_vec_ > vec3_ > vec3_ > int_vec_ > int_vec_;
    partitioned_entity_.name("partitioned_entity");

    partitioned_entities2_ %=
        omit[size_t_[reserve(at_c<0>(_val), _1), _a = _1]] >
        omit[size_t_[reserve(at_c<1>(_val), _1), _b = _1]] >
        omit[size_t_[reserve(at_c<2>(_val), _1), _c = _1]] >
        omit[size_t_[reserve(at_c<3>(_val), _1), _d = _1]] >
        repeat(_a)[partitioned_point_entity_] >
        repeat(_b)[partitioned_entity_] > repeat(_c)[partitioned_entity_] >
        repeat(_d)[partitioned_entity_];
    partitioned_entities2_.name("partitioned_entities2");

    partitioned_entities_ %= "$PartitionedEntities" > size_t_ >
                             ghost_entities_ > partitioned_entities2_ >
                             "$EndPartitionedEntities";
    partitioned_entities_.name("partitioned_entities");

    // nodes:
    node_block_ %=
        int_ > int_ > int_ >
        omit[size_t_[phoenix::resize(at_c<3>(_val), _1), _a = _1, _b = 0]] >
        omit[repeat(_a)[size_t_[at_c<0>(at_c<3>(_val)[_b++]) = _1]]] >
        eps[_b = 0] >
        omit[repeat(_a)[vec3_[at_c<1>(at_c<3>(_val)[_b++]) = _1]]];
    node_block_.name("node_block");

    nodes_ %= "$Nodes" > omit[size_t_[reserve(at_c<3>(_val), _1), _a = _1]] >
              size_t_ > size_t_ > size_t_ > repeat(_a)[node_block_] >
              "$EndNodes";
    nodes_.name("nodes");

    // elements:
    element_block_ %=
        int_ > int_ > int_ >
        omit[size_t_[reserve(at_c<3>(_val), _1), _a = _1]] >
        repeat(_a)[size_t_ > repeat(numNodesAdapted(at_c<2>(_val)))[size_t_]];
    element_block_.name("element_block");

    elements_ %= "$Elements" >
                 omit[size_t_[reserve(at_c<3>(_val), _1), _a = _1]] > size_t_ >
                 size_t_ > size_t_ > repeat(_a)[element_block_] >
                 "$EndElements";
    elements_.name("elements");

    // periodic link
    matrix4d_ %= double_ > double_ > double_ > double_ > double_ > double_ >
                 double_ > double_ > double_ > double_ > double_ > double_ >
                 double_ > double_ > double_ > double_;
    matrix4d_.name("matrix4d");

    periodic_link_ %= int_ > int_ > int_ > ('0' | ("16" > matrix4d_)) >
                      omit[size_t_[reserve(at_c<4>(_val), _1), _a = _1]] >
                      repeat(_a)[size_t_ > size_t_];
    periodic_link_.name("periodic_link");

    periodic_links_ %= "$Periodic" > omit[size_t_[reserve(_val, _1), _a = _1]] >
                       repeat(_a)[periodic_link_] > "$EndPeriodic";
    periodic_links_.name("periodic_links");

    // ghost elements:
    ghost_element_ %= size_t_ > int_ >
                      omit[size_t_[reserve(at_c<2>(_val), _1), _a = _1]] >
                      repeat(_a)[int_];
    ghost_element_.name("ghost_element");

    ghost_elements_ %= "$GhostElements" >
                       omit[size_t_[reserve(_val, _1), _a = _1]] >
                       repeat(_a)[ghost_element_] > "$EndGhostElements";
    ghost_elements_.name("ghost_elements");

    // The whole file
    entry_ %= *comment_ >> -(physical_name_vector_ >> *comment_) >> entities_ >>
              *comment_ >> -(partitioned_entities_ >> *comment_) >> nodes_ >>
              *comment_ >> elements_ >> *comment_ >>
              -(periodic_links_ >> *comment_) >>
              -(ghost_elements_ >> *comment_);
    entry_.name("entry");

    qi::on_error<qi::fail>(entry_, error_handler_(_1, _2, _3, _4));
  }

  qi::rule<ITERATOR, std::string(), ascii::space_type> quoted_string_;
  qi::rule<ITERATOR, std::string()> start_comment_;
  qi::rule<ITERATOR, qi::locals<std::string>, ascii::space_type> comment_;
  qi::rule<ITERATOR, Eigen::Vector3d(), ascii::space_type> vec3_;
  qi::rule<ITERATOR, std::vector<int>, qi::locals<std::size_t>,
           ascii::space_type>
      int_vec_;

  qi::rule<ITERATOR, GMshFileV4::PhysicalName(), ascii::space_type>
      physical_name_;
  qi::rule<ITERATOR, std::vector<GMshFileV4::PhysicalName>(),
           qi::locals<std::size_t>, ascii::space_type>
      physical_name_vector_;

  qi::rule<ITERATOR, boost::fusion::adapted::MshFileV4Adapted(),
           ascii::space_type>
      entry_;

  qi::rule<ITERATOR, GMshFileV4::PointEntity(), ascii::space_type>
      point_entity_;

  qi::rule<ITERATOR, GMshFileV4::Entity(), ascii::space_type> entity_;

  qi::rule<
      ITERATOR,
      std::tuple<
          std::vector<GMshFileV4::PointEntity>, std::vector<GMshFileV4::Entity>,
          std::vector<GMshFileV4::Entity>, std::vector<GMshFileV4::Entity>>(),
      qi::locals<std::size_t, std::size_t, std::size_t, std::size_t>,
      ascii::space_type>
      entities_;

  qi::rule<ITERATOR, std::vector<GMshFileV4::GhostEntity>(),
           qi::locals<std::size_t>, ascii::space_type>
      ghost_entities_;

  qi::rule<ITERATOR, GMshFileV4::PartitionedPointEntity(), ascii::space_type>
      partitioned_point_entity_;

  qi::rule<ITERATOR, GMshFileV4::PartitionedEntity(), ascii::space_type>
      partitioned_entity_;

  qi::rule<ITERATOR,
           std::tuple<std::vector<GMshFileV4::PartitionedPointEntity>,
                      std::vector<GMshFileV4::PartitionedEntity>,
                      std::vector<GMshFileV4::PartitionedEntity>,
                      std::vector<GMshFileV4::PartitionedEntity>>(),
           qi::locals<std::size_t, std::size_t, std::size_t, std::size_t>,
           ascii::space_type>
      partitioned_entities2_;

  qi::rule<ITERATOR, GMshFileV4::PartitionedEntities(), ascii::space_type>
      partitioned_entities_;

  qi::rule<ITERATOR, GMshFileV4::NodeBlock(),
           qi::locals<std::size_t, std::size_t>, ascii::space_type>
      node_block_;

  qi::rule<ITERATOR, GMshFileV4::Nodes(), qi::locals<std::size_t>,
           ascii::space_type>
      nodes_;

  qi::rule<ITERATOR, GMshFileV4::ElementBlock(), qi::locals<std::size_t>,
           ascii::space_type>
      element_block_;

  qi::rule<ITERATOR, GMshFileV4::Elements(), qi::locals<std::size_t>,
           ascii::space_type>
      elements_;

  qi::rule<ITERATOR, Eigen::Matrix4d, ascii::space_type> matrix4d_;

  qi::rule<ITERATOR, GMshFileV4::PeriodicLink(), qi::locals<std::size_t>,
           ascii::space_type>
      periodic_link_;

  qi::rule<ITERATOR, std::vector<GMshFileV4::PeriodicLink>(),
           qi::locals<std::size_t>, ascii::space_type>
      periodic_links_;

  qi::rule<ITERATOR, GMshFileV4::GhostElement(), qi::locals<std::size_t>,
           ascii::space_type>
      ghost_element_;

  qi::rule<ITERATOR, std::vector<GMshFileV4::GhostElement>(),
           qi::locals<std::size_t>, ascii::space_type>
      ghost_elements_;

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
      std::cout << "Error in MshFileV4! Expecting " << what << " here: \""
                << input << "\"" << std::endl;
    }
  };
  phoenix::function<ErrorHandler> error_handler_;
};

template <class ITERATOR>
struct MshV4GrammarBinary
    : qi::grammar<ITERATOR, boost::fusion::adapted::MshFileV4Adapted(),
                  ascii::space_type> {
  /**
   * @brief Construct a new binary grammar
   * @tparam DWORD type of dword
   * @tparam QWORD type of qword
   * @tparam DOBULE type of double_
   * @param dword The spirit terminal that represents a double word (32bit),
   * either little_dword or big_dword
   * @param qword the binary spirit terminal that represents a qword (64),
   * either little_qword or big_qword
   * @param bin_double The binary spirit terminal that represents a double (64
   * bits), is either a little_bin_double or big_bin_double
   */
  template <class DWORD, class QWORD, class DOUBLE>
  explicit MshV4GrammarBinary(const DWORD& dword, const QWORD& qword,
                              const DOUBLE& bin_double)
      : MshV4GrammarBinary::base_type(entry_, "Msh File") {
    using phoenix::at_c;
    using phoenix::push_back;
    using phoenix::reserve;
    using phoenix::val;
    using qi::_val;
    using qi::attr;
    using qi::char_;
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
    using qi::labels::_b;
    using qi::labels::_c;
    using qi::labels::_d;

    auto size_t_ = qi::ulong_long;

    // General Parsers:
    quoted_string_ %= lexeme['"' >> +(char_ - '"') >> '"'];
    quoted_string_.name("string");
    vec3_ %= bin_double > bin_double > bin_double;
    vec3_.name("vec3");
    int_vec_ %= omit[qword[reserve(_val, _1), _a = _1]] > repeat(_a)[dword];
    int_vec_.name("int_vec");
    start_comment_ %= !lit("$PhysicalNames") >> !lit("$Entities") >>
                      !lit("$PartitionedEntities") >> !lit("$Nodes") >>
                      !lit("$Elements") >> !lit("$Periodic") >>
                      !lit("$GhostElements") >>
                      (lit('$') >> (+(char_ - qi::eol)));
    start_comment_.name("Start of Comment");
    comment_ %= start_comment_[_a = qi::_1] > *(char_ - '$') >> "$End" >>
                qi::string(_a);
    comment_.name("comment");
    qi::on_error<qi::fail>(comment_, error_handler_(_1, _2, _3, _4));

    // Physical names
    physical_name_ %= int_ > int_ > quoted_string_;
    physical_name_.name("Physical Name");
    qi::on_error<qi::fail>(physical_name_, error_handler_(_1, _2, _3, _4));
    physical_name_vector_ %= "$PhysicalNames" >
                             omit[size_t_[reserve(_val, _1), _a = _1]] >
                             repeat(_a)[physical_name_] > "$EndPhysicalNames";
    physical_name_vector_.name("$PhysicalNames");

    // entities:
    point_entity_ %= dword > vec3_ > int_vec_;
    point_entity_.name("Point entity");

    entity_ %= dword > vec3_ > vec3_ > int_vec_ > int_vec_;
    entity_.name("entity");

    entities_ %= "$Entities\n" >
                 omit[qword[reserve(phoenix::at_c<0>(_val), _1), _a = _1]] >
                 omit[qword[reserve(phoenix::at_c<1>(_val), _1), _b = _1]] >
                 omit[qword[reserve(phoenix::at_c<2>(_val), _1), _c = _1]] >
                 omit[qword[reserve(phoenix::at_c<3>(_val), _1), _d = _1]] >
                 repeat(_a)[point_entity_] > repeat(_b)[entity_] >
                 repeat(_c)[entity_] > repeat(_d)[entity_] > "\n$EndEntities";
    entities_.name("$Entities");

    // PartitionedEntities
    ghost_entities_ %=
        omit[qword[reserve(_val, _1), _a = _1]] > repeat(_a)[dword > dword];
    ghost_entities_.name("ghost_entities");
    partitioned_point_entity_ %=
        dword > dword > dword > int_vec_ > vec3_ > int_vec_;
    partitioned_point_entity_.name("partitioned_point_entity");

    partitioned_entity_ %=
        dword > dword > dword > int_vec_ > vec3_ > vec3_ > int_vec_ > int_vec_;
    partitioned_entity_.name("partitioned_entity");

    partitioned_entities2_ %= omit[qword[reserve(at_c<0>(_val), _1), _a = _1]] >
                              omit[qword[reserve(at_c<1>(_val), _1), _b = _1]] >
                              omit[qword[reserve(at_c<2>(_val), _1), _c = _1]] >
                              omit[qword[reserve(at_c<3>(_val), _1), _d = _1]] >
                              repeat(_a)[partitioned_point_entity_] >
                              repeat(_b)[partitioned_entity_] >
                              repeat(_c)[partitioned_entity_] >
                              repeat(_d)[partitioned_entity_];
    partitioned_entities2_.name("partitioned_entities2");

    partitioned_entities_ %= "$PartitionedEntities\n" > qword >
                             ghost_entities_ > partitioned_entities2_ >
                             "\n$EndPartitionedEntities";
    partitioned_entities_.name("partitioned_entities");

    // nodes:
    node_block_ %=
        dword > dword > dword >
        omit[qword[phoenix::resize(at_c<3>(_val), _1), _a = _1, _b = 0]] >
        omit[repeat(_a)[qword[at_c<0>(at_c<3>(_val)[_b++]) = _1]]] >
        eps[_b = 0] >
        omit[repeat(_a)[vec3_[at_c<1>(at_c<3>(_val)[_b++]) = _1]]];
    node_block_.name("node_block");

    nodes_ %= "$Nodes\n" > omit[qword[reserve(at_c<3>(_val), _1), _a = _1]] >
              qword > qword > qword > repeat(_a)[node_block_] > "\n$EndNodes";
    nodes_.name("nodes");

    // elements:
    element_block_ %=
        dword > dword > dword >
        omit[qword[reserve(at_c<3>(_val), _1), _a = _1]] >
        repeat(_a)[qword > repeat(numNodesAdapted(at_c<2>(_val)))[qword]];
    element_block_.name("element_block");

    elements_ %= "$Elements\n" >
                 omit[qword[reserve(at_c<3>(_val), _1), _a = _1]] > qword >
                 qword > qword > repeat(_a)[element_block_] > "\n$EndElements";
    elements_.name("elements");

    // periodic link
    matrix4d_ %= bin_double > bin_double > bin_double > bin_double >
                 bin_double > bin_double > bin_double > bin_double >
                 bin_double > bin_double > bin_double > bin_double >
                 bin_double > bin_double > bin_double > bin_double;
    matrix4d_.name("matrix4d");

    periodic_link_ %= dword > dword > dword >
                      (qword(0) | (qword(16) > matrix4d_)) >
                      omit[qword[reserve(at_c<4>(_val), _1), _a = _1]] >
                      repeat(_a)[qword > qword];
    periodic_link_.name("periodic_link");

    periodic_links_ %= "$Periodic\n" > omit[qword[reserve(_val, _1), _a = _1]] >
                       repeat(_a)[periodic_link_] > "\n$EndPeriodic";
    periodic_links_.name("periodic_links");

    // ghost elements:
    ghost_element_ %= qword > dword >
                      omit[qword[reserve(at_c<2>(_val), _1), _a = _1]] >
                      repeat(_a)[dword];
    ghost_element_.name("ghost_element");

    ghost_elements_ %= "$GhostElements\n" >
                       omit[qword[reserve(_val, _1), _a = _1]] >
                       repeat(_a)[ghost_element_] > "\n$EndGhostElements";
    ghost_elements_.name("ghost_elements");

    // The whole file
    entry_ %= *comment_ >> -(physical_name_vector_ >> *comment_) >> entities_ >>
              *comment_ >> -(partitioned_entities_ >> *comment_) >> nodes_ >>
              *comment_ >> elements_ >> *comment_ >>
              -(periodic_links_ >> *comment_) >>
              -(ghost_elements_ >> *comment_);
    entry_.name("entry");

    qi::on_error<qi::fail>(entry_, error_handler_(_1, _2, _3, _4));
  }

  qi::rule<ITERATOR, std::string(), ascii::space_type> quoted_string_;
  qi::rule<ITERATOR, std::string(), ascii::space_type> start_comment_;
  qi::rule<ITERATOR, qi::locals<std::string>, ascii::space_type> comment_;
  qi::rule<ITERATOR, Eigen::Vector3d()> vec3_;
  qi::rule<ITERATOR, std::vector<int>, qi::locals<std::size_t>> int_vec_;

  qi::rule<ITERATOR, GMshFileV4::PhysicalName(), ascii::space_type>
      physical_name_;
  qi::rule<ITERATOR, std::vector<GMshFileV4::PhysicalName>(),
           qi::locals<std::size_t>, ascii::space_type>
      physical_name_vector_;

  qi::rule<ITERATOR, boost::fusion::adapted::MshFileV4Adapted(),
           ascii::space_type>
      entry_;

  qi::rule<ITERATOR, GMshFileV4::PointEntity()> point_entity_;

  qi::rule<ITERATOR, GMshFileV4::Entity()> entity_;

  qi::rule<
      ITERATOR,
      std::tuple<
          std::vector<GMshFileV4::PointEntity>, std::vector<GMshFileV4::Entity>,
          std::vector<GMshFileV4::Entity>, std::vector<GMshFileV4::Entity>>(),
      qi::locals<std::size_t, std::size_t, std::size_t, std::size_t>>
      entities_;

  qi::rule<ITERATOR, std::vector<GMshFileV4::GhostEntity>(),
           qi::locals<std::size_t>>
      ghost_entities_;

  qi::rule<ITERATOR, GMshFileV4::PartitionedPointEntity()>
      partitioned_point_entity_;

  qi::rule<ITERATOR, GMshFileV4::PartitionedEntity()> partitioned_entity_;

  qi::rule<ITERATOR,
           std::tuple<std::vector<GMshFileV4::PartitionedPointEntity>,
                      std::vector<GMshFileV4::PartitionedEntity>,
                      std::vector<GMshFileV4::PartitionedEntity>,
                      std::vector<GMshFileV4::PartitionedEntity>>(),
           qi::locals<std::size_t, std::size_t, std::size_t, std::size_t>>
      partitioned_entities2_;

  qi::rule<ITERATOR, GMshFileV4::PartitionedEntities()> partitioned_entities_;

  qi::rule<ITERATOR, GMshFileV4::NodeBlock(),
           qi::locals<std::size_t, std::size_t>>
      node_block_;

  qi::rule<ITERATOR, GMshFileV4::Nodes(), qi::locals<std::size_t>> nodes_;

  qi::rule<ITERATOR, GMshFileV4::ElementBlock(), qi::locals<std::size_t>>
      element_block_;

  qi::rule<ITERATOR, GMshFileV4::Elements(), qi::locals<std::size_t>> elements_;

  qi::rule<ITERATOR, Eigen::Matrix4d> matrix4d_;

  qi::rule<ITERATOR, GMshFileV4::PeriodicLink(), qi::locals<std::size_t>>
      periodic_link_;

  qi::rule<ITERATOR, std::vector<GMshFileV4::PeriodicLink>(),
           qi::locals<std::size_t>>
      periodic_links_;

  qi::rule<ITERATOR, GMshFileV4::GhostElement(), qi::locals<std::size_t>>
      ghost_element_;

  qi::rule<ITERATOR, std::vector<GMshFileV4::GhostElement>(),
           qi::locals<std::size_t>>
      ghost_elements_;

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
      std::cout << "Error in MshFileV4! Expecting " << what << " here: \""
                << input << "\"" << std::endl;
    }
  };
  phoenix::function<ErrorHandler> error_handler_;
};

}  // namespace

GMshFileV4 ReadGmshFileV4(std::string::const_iterator begin,
                          std::string::const_iterator end, std::string version,
                          bool is_binary, int size_t_size, int one,
                          std::string filename) {
  LF_VERIFY_MSG(version == "4.1", "Only version 4.1 is supported so far");
  LF_VERIFY_MSG(size_t_size == sizeof(std::size_t),
                "size of size_t must be " << sizeof(std::size_t));

  GMshFileV4 result;
  result.version_number = version;
  result.is_binary = is_binary;
  result.size_t_size = size_t_size;

  bool succesful = false;
  if (!is_binary) {
    // Text file
    MshV4GrammarText<std::string::const_iterator> grammar;
    succesful = qi::phrase_parse(begin, end, grammar, ascii::space, result);
  } else if (one == 1) {
    MshV4GrammarBinary<std::string::const_iterator> grammar(
        qi::little_dword, qi::little_qword, qi::little_bin_double);
    succesful = qi::phrase_parse(begin, end, grammar, ascii::space, result);
  } else {
    MshV4GrammarBinary<std::string::const_iterator> grammar(
        qi::big_dword, qi::big_qword, qi::big_bin_double);
    succesful = qi::phrase_parse(begin, end, grammar, ascii::space, result);
  }

  LF_VERIFY_MSG(succesful, "Could not parse file " << filename);
  // LF_VERIFY_MSG(iter == end, "Could not parse all of file " << filename);

  // transpose all periodic matrices because they are read in column-first mode
  // but gmsh writes them in row-first mode
  for (auto& p : result.periodic_links) {
    if (p.affine_transform) {
      p.affine_transform->transposeInPlace();
    }
  }

  return result;
}

}  // namespace lf::io
