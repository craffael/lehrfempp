/**
 * @file
 * @brief Contains the binary grammar for the GmshFileV4. This is separated from
 * the grammar that is used to read the textual file because g++ would allocate
 * too much memory otherwise.
 * @author Raffael Casagrande
 * @date   2019-10-13 09:42:29
 * @copyright MIT License
 */

#include "gmsh_file_v4_detail.h"

namespace lf::io {

namespace /* anonymous */ {
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

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
    vec3_ %= bin_double > bin_double > bin_double;  // NOLINT
    vec3_.name("vec3");
    int_vec_ %= omit[qword[(reserve(_val, _1), _a = _1)]] > repeat(_a)[dword];
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
    physical_name_ %= int_ > int_ > quoted_string_;  // NOLINT
    physical_name_.name("Physical Name");
    qi::on_error<qi::fail>(physical_name_, error_handler_(_1, _2, _3, _4));
    physical_name_vector_ %= "$PhysicalNames" >
                             omit[size_t_[(reserve(_val, _1), _a = _1)]] >
                             repeat(_a)[physical_name_] > "$EndPhysicalNames";
    physical_name_vector_.name("$PhysicalNames");

    // entities:
    point_entity_ %= dword > vec3_ > int_vec_;
    point_entity_.name("Point entity");

    entity_ %= dword > vec3_ > vec3_ > int_vec_ > int_vec_;
    entity_.name("entity");

    entities_ %= "$Entities\n" >
                 omit[qword[(reserve(phoenix::at_c<0>(_val), _1), _a = _1)]] >
                 omit[qword[(reserve(phoenix::at_c<1>(_val), _1), _b = _1)]] >
                 omit[qword[(reserve(phoenix::at_c<2>(_val), _1), _c = _1)]] >
                 omit[qword[(reserve(phoenix::at_c<3>(_val), _1), _d = _1)]] >
                 repeat(_a)[point_entity_] > repeat(_b)[entity_] >
                 repeat(_c)[entity_] > repeat(_d)[entity_] > "\n$EndEntities";
    entities_.name("$Entities");

    // PartitionedEntities
    ghost_entities_ %= omit[qword[(reserve(_val, _1), _a = _1)]] >
                       repeat(_a)[dword > dword];  // NOLINT
    ghost_entities_.name("ghost_entities");
    partitioned_point_entity_ %=
        dword > dword > dword > int_vec_ > vec3_ > int_vec_;  // NOLINT
    partitioned_point_entity_.name("partitioned_point_entity");

    // NOLINTNEXTLINE
    partitioned_entity_ %= dword > dword > dword > int_vec_ > vec3_ > vec3_ >
                           int_vec_ > int_vec_;  // NOLINT
    partitioned_entity_.name("partitioned_entity");

    partitioned_entities2_ %= omit[qword[(reserve(at_c<0>(_val), _1), _a = _1)]] >
                              omit[qword[(reserve(at_c<1>(_val), _1), _b = _1)]] >
                              omit[qword[(reserve(at_c<2>(_val), _1), _c = _1)]] >
                              omit[qword[(reserve(at_c<3>(_val), _1), _d = _1)]] >
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
        dword > dword > dword >  // NOLINT
        omit[qword[(phoenix::resize(at_c<3>(_val), _1), _a = _1, _b = 0)]] >
        omit[repeat(_a)[qword[at_c<0>(at_c<3>(_val)[_b++]) = _1]]] >
        eps[_b = 0] >
        omit[repeat(_a)[vec3_[at_c<1>(at_c<3>(_val)[_b++]) = _1]]];
    node_block_.name("node_block");

    nodes_ %= "$Nodes\n" > omit[qword[(reserve(at_c<3>(_val), _1), _a = _1)]] >
              qword > qword > qword > repeat(_a)[node_block_] > "\n$EndNodes";
    nodes_.name("nodes");

    // elements:
    element_block_ %=
        dword > dword > dword >  // NOLINT
        omit[qword[(reserve(at_c<3>(_val), _1), _a = _1)]] >
        repeat(_a)[qword > repeat(numNodesAdapted(at_c<2>(_val)))[qword]];
    element_block_.name("element_block");

    elements_ %= "$Elements\n" >
                 omit[qword[(reserve(at_c<3>(_val), _1), _a = _1)]] > qword >
                 qword > qword > repeat(_a)[element_block_] > "\n$EndElements";
    elements_.name("elements");

    // periodic link
    // NOLINTNEXTLINE
    matrix4d_ %= bin_double > bin_double > bin_double > bin_double >
                 bin_double > bin_double > bin_double > bin_double >
                 bin_double > bin_double > bin_double > bin_double >
                 bin_double > bin_double > bin_double > bin_double;
    matrix4d_.name("matrix4d");

    // NOLINTNEXTLINE
    periodic_link_ %= dword > dword > dword >
                      (qword(0) | (qword(16) > matrix4d_)) >
                      omit[qword[(reserve(at_c<4>(_val), _1), _a = _1)]] >
                      repeat(_a)[qword > qword];  // NOLINT
    periodic_link_.name("periodic_link");

    periodic_links_ %= "$Periodic\n" > omit[qword[(reserve(_val, _1), _a = _1)]] >
                       repeat(_a)[periodic_link_] > "\n$EndPeriodic";
    periodic_links_.name("periodic_links");

    // ghost elements:
    ghost_element_ %= qword > dword >
                      omit[qword[(reserve(at_c<2>(_val), _1), _a = _1)]] >
                      repeat(_a)[dword];
    ghost_element_.name("ghost_element");

    ghost_elements_ %= "$GhostElements\n" >
                       omit[qword[(reserve(_val, _1), _a = _1)]] >
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

namespace detail {
bool ParseGmshFileV4Binary(std::string::const_iterator begin,
                           std::string::const_iterator end, int one,
                           GMshFileV4* result) {
  if (one == 1) {
    const MshV4GrammarBinary<std::string::const_iterator> grammar(
        qi::little_dword, qi::little_qword, qi::little_bin_double);
    return qi::phrase_parse(begin, end, grammar, ascii::space, *result);
  }
  const MshV4GrammarBinary<std::string::const_iterator> grammar(
      qi::big_dword, qi::big_qword, qi::big_bin_double);
  return qi::phrase_parse(begin, end, grammar, ascii::space, *result);
}
}  // namespace detail

}  // namespace lf::io
