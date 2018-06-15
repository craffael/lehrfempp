#include "mesh.h"
#include <lf/geometry/geometry.h>
#include "entity.h"

namespace lf::mesh::hybrid2d {

Mesh::Mesh(char dim_world, std::vector<Eigen::VectorXd> nodes,
           std::vector<std::tuple<std::vector<size_type>,
                                  std::unique_ptr<geometry::Geometry>>>
               elements)
    : dim_world_(dim_world) {
  // 1) Add all nodes:
  entities2_.reserve(nodes.size());  // allocated memory for all vertex objects
  for (auto& nodeCoord : nodes) {
    // A vertex does not have any sub-entities; dummy argument
    std::array<std::vector<size_type>, 0> sub_entities{};
    // The geometry information for a vertex consists of its coordinates only.
    // Add another vertex to the array of vertices setting its index to the
    // position in the array.
    entities2_.emplace_back(this, entities2_.size(),
                            std::make_unique<geometry::Point>(nodeCoord),
                            sub_entities);
  }

  // Initialize vector of cells
  entities0_.reserve(elements.size());
  // 2) put all edges of all elements into a vector (so there are duplicates)
  // + initialize entities of codim0 with nodes (but don't set edges yet)
  std::vector<std::tuple<std::array<size_type, 2>, size_type, char>>
      element_edges{};
  // Reserve auxiliary vector for edges
  element_edges.reserve(4 * elements.size());
  // Run through all cells, given index arrays of vertices and geometry objects
  for (auto& tuple : elements) {
    // Index array of vertices
    auto& element_nodes = std::get<0>(tuple);
    // The number of vertices determines the cell type: triangle of
    // quadrilateral
    base::RefEl ref_el =
        element_nodes.size() == 3 ? base::RefEl::kTria() : base::RefEl::kQuad();
    // run through the edges of the current cell
    for (char i = 0; i < ref_el.NumSubEntities(1); ++i) {
      // Store the global indices of the endpoints of the current edge
      // in the variables `start` and and `end`
      size_type start =
          element_nodes[ref_el.SubSubEntity2SubEntity(1, i, 1, 0)];
      size_type end = element_nodes[ref_el.SubSubEntity2SubEntity(1, i, 1, 1)];
      // Canonical orientation of an edge: from endpoint with lower
      // index to endpoint with higher index
      std::array<size_type, 2> nodes_ordered = {start, end};
      if (start > end) {
        std::swap(nodes_ordered[0], nodes_ordered[1]);
      }
      // register edge in auxiliary array
      element_edges.emplace_back(nodes_ordered, entities0_.size(), i);
    }

    // Index numbers of sub-entities of current cell
    // Since edges have not been constructed yet, leave their indices
    // unspecified.
    std::array<std::vector<size_type>, 2> sub_entities = {
        {std::vector<size_type>(ref_el.NumSubEntities(1), size_type(-1)), {}}};
    // Register index numbers of vertices with the current cell
    sub_entities[1].reserve(element_nodes.size());
    for (auto node_nr : element_nodes) {
      sub_entities[1].push_back(node_nr);
    }
    // Create new Entity object for current cell
    entities0_.emplace_back(this, entities0_.size(),
                            std::move(std::get<1>(tuple)), sub_entities);
  }

  // 3) sort vector -> edges with same start/end node are adjacent in the
  // array after sorting
  std::sort(element_edges.begin(), element_edges.end());

  // 4) add edges and register them with elements
  size_type begin = 0;
  entities1_.reserve(element_edges.size() / 3 * 2);
  // Traverse auxiliary array containing duplicates of edges
  for (size_type end = 0; end < element_edges.size(); ++end) {
    auto& end_edge = element_edges[end];
    if (begin == end) {
      // We have found a new edge.
      std::array<std::vector<size_type>, 1> sub_entities;
      // Adjacent element with the lowest index number
      // sub-entity
      auto& element = entities0_[std::get<1>(element_edges[begin])];
      // index of edge in element
      auto iie = std::get<2>(element_edges[begin]);
      // Global index of the endpoint 0 of the edge in the element
      auto nodeIndex = element.RefEl().SubSubEntity2SubEntity(1, iie, 1, 0);
      sub_entities[0].push_back(element.sub_entities_[1][nodeIndex]);
      // Global index of the endpoint 1 of the edge in the element
      nodeIndex = element.RefEl().SubSubEntity2SubEntity(1, iie, 1, 1);
      sub_entities[0].push_back(element.sub_entities_[1][nodeIndex]);
      // Note that we have oriented the edge so that it agrees with the first
      // (lower) element's sub-entity

      // The shape of the edge is just inherited from the cell.
      auto geom = element.Geometry()->subGeometry(1, iie);
      //  Create Entity object for the new edge
      entities1_.emplace_back(this, entities1_.size(), std::move(geom),
                              sub_entities);
    }

    // register edge at element. Hmm, fairly intrusive !
    entities0_[std::get<1>(end_edge)].sub_entities_[0][std::get<2>(end_edge)] =
        entities1_.size() - 1;
    // Another new edge is detected
    if (end + 1 < element_edges.size() &&
        std::get<0>(element_edges[begin]) !=
            std::get<0>(element_edges[end + 1])) {
      begin = end + 1;
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

Mesh::size_type Mesh::Size(char codim) const {
  LF_ASSERT_MSG(codim >= 0, "codim negative.");
  LF_ASSERT_MSG(codim <= dim_world_, "codim > dimWorld.");

  switch (codim) {
    case 0:
      return entities0_.size();
    case 1:
      return entities1_.size();
    case 2:
      return entities2_.size();
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
}  // namespace lf::mesh::hybrid2d
