#include "mesh.h"
#include <lf/geometry/geometry.h>
#include "entity.h"

#include <map>
#include <memory>
#include <utility>

namespace lf::mesh::hybrid2d {

Mesh::Mesh(char dim_world, std::vector<Eigen::VectorXd> nodes,
           std::vector<std::tuple<std::array<size_type, 2>,
                                  std::unique_ptr<geometry::Geometry>>>
               edges,
           std::vector<std::tuple<std::array<size_type, 4>,
                                  std::unique_ptr<geometry::Geometry>>>
               elements)
    : dim_world_(dim_world), entities1_(edges.size()) {
  // 1) Add all nodes:
  entities2_.reserve(nodes.size());  // allocated memory for all vertex objects
  for (auto &nodeCoord : nodes) {
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

  // element_edges[i][0] : node indices of the edge, sorted.
  // element_edges[i][1] = false => edge was manually speficied by user
  //                     = true  => edge is sub-entity of an element
  // element_edges[i][2] : if element_edges[i][1]=false, this contains the
  //                       zero-based index of the edge.
  //                     : if element_edges[i][1]=true, this contains the index
  //                       of the element to which the edge belongs
  // element_edges[i][3] : if element_edges[i][1]=true, this contains the index
  //                       of the edge w.r.t. to the element.
  //                       if element_edges[i][1]=false, it is undefined.
  std::vector<std::tuple<std::array<size_type, 2>, bool, size_type, char>>
      element_edges{};
  // Allocate enough space in the vector so it doesn't resize while inserting.
  element_edges.reserve(4 * elements.size() + edges.size());

  // put all edges that have been added explicitly by the user to element_edges:
  for (size_type i = 0; i < edges.size(); ++i) {
    auto &edge_nodes = std::get<0>(edges[i]);
    std::array<size_type, 2> nodes_ordered{{edge_nodes[0], edge_nodes[1]}};
    if (edge_nodes[0] > edge_nodes[1]) {
      std::swap(nodes_ordered[0], nodes_ordered[1]);
    }
    element_edges.emplace_back(nodes_ordered, false, i, -1);
  }

  // Run through all cells, given index arrays of vertices and geometry objects
  for (auto &tuple : elements) {
    // Index array of vertices
    auto &element_nodes = std::get<0>(tuple);
    // The number of vertices determines the cell type: triangle of
    // quadrilateral
    base::RefEl ref_el = element_nodes[3] == size_type(-1)
                             ? base::RefEl::kTria()
                             : base::RefEl::kQuad();
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
      element_edges.emplace_back(nodes_ordered, true, entities0_.size(), i);
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
  entities1_.reserve(std::max(edges.size(), element_edges.size() / 3 * 2));
  // Traverse auxiliary array containing duplicates of edges
  for (size_type end = 0; end < element_edges.size(); ++end) {
    auto &end_edge = element_edges[end];
    std::unique_ptr<geometry::Geometry> geom;
    if (begin == end) {
      // We have found a new edge.
      std::array<std::vector<size_type>, 1> sub_entities;
      if (std::get<1>(element_edges[begin])) {
        // This edge has not been added explictly by the user -> we deduce
        // its properties from the super element with the lower index

        // Adjacent element with the lowest index number
        // sub-entity
        auto &element = entities0_[std::get<2>(element_edges[begin])];
        // index of edge in element
        auto iie = std::get<3>(element_edges[begin]);
        // Global index of the endpoint 0 of the edge in the element
        auto nodeIndex = element.RefEl().SubSubEntity2SubEntity(1, iie, 1, 0);
        sub_entities[0].push_back(element.sub_entities_[1][nodeIndex]);
        // Global index of the endpoint 1 of the edge in the element
        nodeIndex = element.RefEl().SubSubEntity2SubEntity(1, iie, 1, 1);
        sub_entities[0].push_back(element.sub_entities_[1][nodeIndex]);
        // Note that we have oriented the edge so that it agrees with the first
        // (lower) element's sub-entity

        // The shape of the edge is just inherited from the cell.
        geom = element.Geometry()->SubGeometry(1, iie);

        //  Create Entity object for the new edge
        entities1_.emplace_back(this, entities1_.size(), std::move(geom),
                                sub_entities);
      } else {
        // This edge has been added explictly by the user
        auto &original_edge = edges[std::get<2>(element_edges[begin])];
        LF_ASSERT_MSG(begin + 1 < element_edges.size() &&
                          std::get<0>(element_edges[begin]) ==
                              std::get<0>(element_edges[begin + 1]),
                      "The explictly added entity with codim=1 and nodes "
                          << std::get<0>(original_edge)[0] << ", "
                          << std::get<0>(original_edge)[1]
                          << " doesn't belong to an element.");
        LF_ASSERT_MSG(std::get<1>(element_edges[begin + 1]),
                      "The explictly added entity with codim=1 and nodes "
                          << std::get<0>(original_edge)[0] << ", "
                          << std::get<0>(original_edge)[1]
                          << "has been added twice.");
        auto node_index = std::get<0>(original_edge)[0];
        sub_entities[0].push_back(node_index);
        node_index = std::get<0>(original_edge)[1];
        sub_entities[0].push_back(node_index);

        if (std::get<1>(original_edge)) {
          // The user also added a geometry object:
          geom = std::move(std::get<1>(original_edge));
        } else {
          // The user didn't specify a geometry object
          // -> We find the super entity that contains this entity and
          // construct the geometry object from it.

          auto &other_edge = element_edges[begin + 1];
          auto &super_element = entities0_[std::get<2>(other_edge)];
          auto iise = std::get<3>(other_edge);
          auto sub_index =
              super_element.RefEl().SubSubEntity2SubEntity(1, iise, 1, 0);
          auto first_node_index = super_element.sub_entities_[1][sub_index];
          geom = super_element.Geometry()->SubGeometry(1, iise);
          if (first_node_index != std::get<0>(original_edge)[0]) {
            // the edge in the super_element has the opposite orientation as the
            // edge that was specified by the user -> switch node order:
            std::swap(sub_entities[0][0], sub_entities[0][1]);
          }
        }
        // put entity element at the right location:
        entities1_[std::get<2>(end_edge)] = Entity<1>(
            this, std::get<2>(end_edge), std::move(geom), sub_entities);
      }
    }

    if (std::get<1>(element_edges[end])) {
      // register edge at element. Hmm, fairly intrusive !
      size_type edge_index;
      if (!std::get<1>(element_edges[begin])) {
        // the edge was specified manually by the user:
        edge_index = std::get<2>(element_edges[begin]);
      } else {
        edge_index = entities1_.size() - 1;
      }
      entities0_[std::get<2>(end_edge)]
          .sub_entities_[0][std::get<3>(end_edge)] = edge_index;
    }
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

Mesh::size_type Mesh::Index(const mesh::Entity &e) const {
  switch (e.Codim()) {
    case 0:
      return dynamic_cast<const Entity<0> &>(e).index_;
    case 1:
      return dynamic_cast<const Entity<1> &>(e).index_;
    case 2:
      return dynamic_cast<const Entity<2> &>(e).index_;
    default:
      LF_VERIFY_MSG(false,
                    "Something is horribyl wrong, this entity has codim = " +
                        std::to_string(e.Codim()));
  }
}

bool Mesh::Contains(const mesh::Entity &e) const {
  switch (e.Codim()) {
    case 0:
      return &e >= &entities0_.front() && &e <= &entities0_.back();
    case 1:
      return &e >= &entities1_.front() && &e <= &entities1_.back();
    case 2:
      return &e >= &entities2_.front() && &e <= &entities2_.back();
    default:
      return false;
  }
}
}  // namespace lf::mesh::hybrid2d
