#include "mesh.h"
#include "entity.h"

namespace lf::mesh::hybrid2d {

Mesh::Mesh(char dim_world, std::vector<Eigen::VectorXd> nodes,
           std::vector<std::tuple<
             std::vector<size_type>, std::unique_ptr<Geometry>>>
           elements)
  : dim_world_(dim_world),
    entities0_(),
    entities1_(),
    entities2_() {

  // 1) Add all nodes:
  entities2_.reserve(nodes.size());
  for (auto& nodeCoord : nodes) {
    std::array<std::vector<size_type>, 0> sub_entities;
    entities2_.emplace_back(this, entities2_.size(), nullptr, sub_entities);
  }

  // 2) put all edges of all elements into a vector (so there are duplicates)
  // + initialize entities of codim0 with nodes (but don't set edges yet)
  entities0_.reserve(elements.size());
  std::vector<std::tuple<std::array<size_type, 2>, size_type, char>> element_edges{};
  element_edges.reserve(4 * elements.size());

  for (auto& tuple : elements) {
    auto& element_nodes = std::get<0>(tuple);

    base::RefEl ref_el = element_nodes.size() == 3
                           ? base::RefEl::kTria()
                           : base::RefEl::kQuad();
    for (char i = 0; i < ref_el.NumSubEntities(1); ++i) {
      size_type start = element_nodes[ref_el.SubSubEntity2SubEntity(1, i, 1, 0)];
      size_type end = element_nodes[ref_el.SubSubEntity2SubEntity(1, i, 1, 1)];
      std::array<size_type, 2> nodes_ordered = {start, end};
      if (start > end) std::swap(nodes_ordered[0], nodes_ordered[1]);
      element_edges.emplace_back(nodes_ordered, entities0_.size(), i);
    }

    // initialize sub_entities with nodes of element:
    std::array<std::vector<size_type>, 2> sub_entities = {
        {std::vector<size_type>(ref_el.NumSubEntities(1), size_type(-1)), {}}};
    sub_entities[1].reserve(element_nodes.size());
    for (auto node_nr : element_nodes) {
      sub_entities[1].push_back(node_nr);
    }
    entities0_.emplace_back(this, entities0_.size(), nullptr, sub_entities);
  }

  // 3) sort vector -> edges with same start/end node are adjacent
  std::sort(element_edges.begin(), element_edges.end());

  // 4) add edges and register them with elements
  size_type begin = 0;
  entities1_.reserve(element_edges.size()/3*2);
  for (size_type end = 0; end < element_edges.size(); ++end) {
    auto& end_edge = element_edges[end];
    if(begin == end) {
      // add edge
      std::array<std::vector<size_type>, 1> sub_entities;
      // orient the edge always so it agrees with the first (lower) element's sub-entity
      auto& element = entities0_[std::get<1>(element_edges[begin])];
     
      sub_entities[0].push_back(element.sub_entities_[1][element.RefEl().SubSubEntity2SubEntity(1, std::get<2>(element_edges[begin]), 1, 0)]);
      sub_entities[0].push_back(element.sub_entities_[1][element.RefEl().SubSubEntity2SubEntity(1, std::get<2>(element_edges[begin]), 1, 1)]);
      entities1_.emplace_back(this, entities1_.size(), nullptr, sub_entities);
    }

    // register edge at element
    entities0_[std::get<1>(end_edge)].sub_entities_[0][std::get<2>(end_edge)] = entities1_.size() - 1;

    if (end + 1 < element_edges.size() &&
        std::get<0>(element_edges[begin]) != std::get<0>(element_edges[end + 1])
    ) {
      begin = end+1;
    }
  }

}

base::ForwardRange<const mesh::Entity> Mesh::Entities(char codim) const {
  LF_ASSERT_MSG(codim >= 0, "codim negative.");
  LF_ASSERT_MSG(codim <= dim_world_, "codim > dimWorld.");

  switch (codim) {
    case 0:
      return {entities0_.begin(), entities0_.end()};
    case 1:
      return {entities1_.begin(), entities1_.end()};
    case 2:
      return {entities2_.begin(), entities2_.end()};
    default:
    LF_VERIFY_MSG(false, "Something is horribyl wrong, codim = " +
        std::to_string(codim) + " is out of bounds.");
  }
}

Mesh::size_type Mesh::Index(const mesh::Entity& e) const {
  switch (e.Codim()) {
    case 0:
      return dynamic_cast<const Entity<0>&>(e).index_;
    case 1:
      return dynamic_cast<const Entity<1>&>(e).index_;
    case 2:
      return dynamic_cast<const Entity<2>&>(e).index_;
    default:
    LF_VERIFY_MSG(false,
        "Something is horribyl wrong, this entity has codim = " +
        std::to_string(e.Codim()));
  }
}
}
