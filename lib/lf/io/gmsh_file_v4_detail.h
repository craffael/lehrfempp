/**
 * @file
 * @brief Contains some common definition which is included by
 * gmsh_file_v4_text.cc and gmsh_file_v4_binary.cc This file should not be
 * included by code that uses LehrFEM++, it is only needed to compile
 * gmsh_file_v4_binary and gmsh_file_v4_text
 * @author Raffael Casagrande
 * @date   2019-10-13 09:44:03
 * @copyright MIT License
 */

#ifndef __c3071f8127a44f7e8cb57f0b1dd3335a
#define __c3071f8127a44f7e8cb57f0b1dd3335a

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

BOOST_FUSION_ADAPT_STRUCT(lf::io::GMshFileV4::Nodes,
                          (std::size_t, num_nodes)(std::size_t,
                                                   min_node_tag)(std::size_t,
                                                                 max_node_tag)(
                              std::vector<lf::io::GMshFileV4::NodeBlock>,
                              node_blocks));

// To prevent comma in preprocessor invocation
using elementMapping_t = std::pair<std::size_t, std::vector<std::size_t>>;

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
// NOLINTNEXTLINE
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

namespace lf::io {
// NOLINTNEXTLINE
BOOST_PHOENIX_ADAPT_FUNCTION(int, numNodesAdapted, NumNodes, 1);
}  // namespace lf::io

/// \cond
namespace boost::spirit::traits {

template <class Enum, class RawValue>
struct assign_to_attribute_from_value<
    Enum, RawValue,
    typename std::enable_if_t<std::is_enum_v<Enum> &&
                              !std::is_same_v<Enum, RawValue>>> {
  static void call(RawValue const& raw, Enum& cat) {
    if constexpr (detail::has_value_type<RawValue>::value) {  // NOLINT
      // specialization for endian::endian
      typename RawValue::value_type value = raw;
      cat = static_cast<Enum>(value);
    } else {  // NOLINT
      cat = static_cast<Enum>(raw);
    }
  }
};

}  // namespace boost::spirit::traits
/// \endcond
#endif  // __c3071f8127a44f7e8cb57f0b1dd3335a
