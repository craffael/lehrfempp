#include "mesh.h"

namespace lf::mesh::hybrid2d {

Mesh::Mesh(char dim_world, std::vector<Eigen::VectorXd> nodes,
           std::vector<std::tuple<
             std::vector<size_t>, std::unique_ptr<Geometry>>>
           elements)
  : dim_world_(dim_world),
    entities0_(elements.size()),
    entities1_(),
    entities2_(nodes.size()) {

  // 1) Add all nodes:
  for (auto& nodeCoord : nodes) {
    std::array<std::vector<mesh::Entity*>, 0> sub_entities;
    entities2_.emplace_back(entities2_.size(), nullptr, sub_entities);
  }

  // 2) put all edges of all elements into a vector (so there are duplicates)
  // + initialize entities of codim0 with nodes (but don't set edges yet)
  std::vector<std::tuple<std::array<size_t, 2>, size_t, char>>
    element_edges(4 * elements.size());

  for (auto& tuple : elements) {
    auto& element_nodes = std::get<0>(tuple);

    base::RefEl ref_el = element_nodes.size() == 3
                           ? base::RefEl::kTria()
                           : base::RefEl::kQuad();
    for (char i = 0; i < ref_el.NumSubEntities(1); ++i) {
      size_t start = element_nodes[ref_el.SubSubEntity2SubEntity(1, i, 1, 0)];
      size_t end = element_nodes[ref_el.SubSubEntity2SubEntity(1, i, 1, 1)];
      std::array<size_t, 2> nodes_ordered = {start, end};
      if (start > end) std::swap(nodes_ordered[0], nodes_ordered[1]);
      element_edges.emplace_back(nodes_ordered, entities0_.size(), i);
    }

    // initialize sub_entities of element:
    std::array<std::vector<mesh::Entity*>, 2> sub_entities = {
        {{ref_el.NumSubEntities(1), nullptr}, {element_nodes.size(), nullptr}}};
    for (auto node_nr : element_nodes) {
      sub_entities[1].push_back(&entities2_[node_nr]);
    }
    entities0_.emplace_back(entities0_.size(), nullptr, sub_entities);
  }

  // 3) sort vector -> edges with same start/end node are adjacent
  std::sort(element_edges.begin(), element_edges.end());

  // 4) add edges and register them with elements
  size_t begin = 0;
  for (size_t end = 0; end < element_edges.size(); ++end) {
    if(begin == end) {
      // add edge
      std::array<std::vector<mesh::Entity*>, 1> sub_entities;
      sub_entities[0].push_back(&entities2_[std::get<0>(element_edges[begin])[0]]);
      sub_entities[0].push_back(&entities2_[std::get<0>(element_edges[begin])[1]]);
      entities1_.emplace_back(entities1_.size(), nullptr, sub_entities);
    }

    // register edge at element
    entities0_[std::get<1>(element_edges.back())].sub_entities_[0][std::get<2>(element_edges.back())] = &entities1_.back();

    if (end + 1 < element_edges.size() &&
        std::get<0>(element_edges[begin]) != std::get<0>(element_edges[end + 1])
    ) {
      begin = end+1;
    }
  }

}
}
