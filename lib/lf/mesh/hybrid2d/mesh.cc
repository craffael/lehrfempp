#include "mesh.h"
#include <lf/geometry/geometry.h>
#include "entity.h"

#include <iostream>
#include <map>
#include <memory>
#include <utility>

namespace lf::mesh::hybrid2d {

using GeometryPtr = std::unique_ptr<geometry::Geometry>;

using EdgeList =
    std::vector<std::pair<std::array<Mesh::size_type, 2>, GeometryPtr>>;
using CellList =
    std::vector<std::pair<std::vector<Mesh::size_type>, GeometryPtr>>;

class EndpointIndexPair {
 public:
  using size_type = Mesh::size_type;
  // Constructor ensures ordering of indices
  EndpointIndexPair(size_type p0, size_type p1) {
    LF_ASSERT_MSG(p0 != p1, "No loops allowed");
    if (p1 > p0) {
      p0_ = p0;
      p1_ = p1;
    } else {
      p0_ = p1;
      p1_ = p0;
    }
  }
  size_type first(void) const { return p0_; }
  size_type second(void) const { return p1_; }
  friend bool operator<(const EndpointIndexPair &e1,
                        const EndpointIndexPair &e2) {
    return ((e1.p0_ == e2.p0_) ? (e1.p1_ < e2.p1_) : (e1.p0_ < e2.p0_));
  }

 private:
  size_type p0_, p1_;  // indices of endpoints
};

/**
 * @brief Construction of mesh from information gathered in a MeshFactory
 * @param nodes sequential container of node coordinates
 * @param edges sequential container of pairs of
                (i) vectors of indices of the nodes of an edge
                (ii) pointers to the geometry object describing an edge
 * @param cells sequential container of pairs of
                (i) vectors of indices of the nodes of a cell
                (ii) pointers to the geometry object for the cell
 * @return TBD
 * @note the position of node information the `nodes` array and of cell
 *        information in the `cells` array, respectively,
 *        determines the interpretation of the index numbers,
 *        that is the n-th node in the container has index n-1.
 *
 */
void mesh_from_node_incidence(std::vector<Eigen::VectorXd> nodes,
                              EdgeList &&edges, CellList cells) {
  using size_type = Mesh::size_type;
  // Index of adjacent cell and local index of edge w.r.t. cell
  using AdjCellInfo = std::pair<size_type, size_type>;
  // Information about cells adjacent to an edge
  using AdjCellsList = std::vector<AdjCellInfo>;
  // Information about edge in auxiliary array
  using EdgeData = std::pair<GeometryPtr, AdjCellsList>;
  // Type of associative auxiliary array for edge information
  using EdgeMap = std::map<EndpointIndexPair, EdgeData>;

  // STEP I: Set up and fill array of nodes
  // In the beginning initialize vector of vertices and do not touch it anymore
  const size_type no_of_nodes(nodes.size());
  /*
     Initialize vector for Node entities of size `no_of_nodes`
  */
  size_type node_index = 0;
  for (const auto &v : nodes) {
    const Eigen::VectorXd &node_coordinates(v);
    /*
       Add Node entity object with index `node_index` at the end
       of vector of Node entities.
    */
    node_index++;
  }

  // STEP II: Initialize array of edges using pointers to
  //          entries of the array of nodes

  // Register supplied edges in auxiliary map data structure
  EdgeMap edge_map;
  for (auto &e : edges) {
    // Node indices of endpoints: the KEY
    std::array<size_type, 2> end_nodes(e.first);
    EndpointIndexPair e_endpoint_idx(end_nodes[0], end_nodes[1]);
    // Store provided geometry information; information on adjacent cells
    // not yet available.
    AdjCellsList empty_cells_list{};
    EdgeData edge_data{std::move(e.second), empty_cells_list};
    EdgeMap::value_type edge_info =
        std::make_pair(e_endpoint_idx, std::move(edge_data));
    std::pair<EdgeMap::iterator, bool> insert_status =
        edge_map.insert(std::move(edge_info));
    if (!insert_status.second) {
      std::cerr << "Duplicate edge " << end_nodes[0] << " <-> " << end_nodes[1]
                << std::endl;
    }
  }  // end loop over predefined edges
  // At this point all predefined edges have been stored in the auxiliary
  // associative array, though without information about adjacent cells

  // Run through cells in order to
  // (i) build edges missing in the list of predefined edges
  // (ii) determine cells adjacent to edges
  size_type cell_index = 0;
  size_type no_of_trilaterals = 0;
  size_type no_of_quadrilaterals = 0;
  for (const auto &c : cells) {
    // node indices of corners of cell c
    const std::vector<size_type> &cell_node_list(c.first);
    // Can be either a trilateral or a quadrilateral
    const size_type no_of_vertices = cell_node_list.size();
    LF_ASSERT_MSG((no_of_vertices == 3) || (no_of_vertices == 4),
                  "Cell with invalid number of edges");
    // Count the different cell types
    if (no_of_vertices == 3)
      no_of_trilaterals++;
    else
      no_of_quadrilaterals++;

    // Visit all edges of the current cell
    for (int j = 0; j < no_of_vertices; j++) {
      EndpointIndexPair c_edge_vertex_indices(
          cell_node_list[j], cell_node_list[(j + 1) % no_of_vertices]);
      // Store number of cell and the local index j of the edge
      AdjCellInfo edge_cell_info{cell_index, j};
      // Check whether edge exists already
      EdgeMap::iterator edge_ptr = edge_map.find(c_edge_vertex_indices);
      if (edge_ptr == edge_map.end()) {
        // Edge does not exist: first create its geometry from the cell
        // by invoking subGeometry()
        // Does not work at this point because to type mismatch
        GeometryPtr edge_geo_ptr(nullptr);

        AdjCellsList single_cell_list{edge_cell_info};
        EdgeData edge_data{std::move(edge_geo_ptr), single_cell_list};
        std::pair<EdgeMap::iterator, bool> insert_status = edge_map.insert(
            std::make_pair(c_edge_vertex_indices, std::move(edge_data)));
        LF_ASSERT_MSG(insert_status.second, "Duplicate not found earlier!");
      } else {
        // Store information about neighboring cell
        AdjCellsList &edge_adj_cellslist((edge_ptr->second).second);
        edge_adj_cellslist.push_back(edge_cell_info);
        // Note that the geometry for this edge is fixed already
      }
    }
    cell_index++;
  }  // end loop over cells

  // Run through the entire associative container for edges
  // and build edge Entities
  // This is the length to be reserved for the edge vector
  size_type no_of_edges = edge_map.size();
  /*
    Initialized vector of Edge entities here
   */
  size_type edge_index = 0;
  for (EdgeMap::value_type &edge : edge_map) {
    // Indices of the two endpoints of the current edge
    // Use this to obtain pointers/references to nodes
    // from the vector of nodes
    const size_type p0(edge.first.first());   // index of first endpoint
    const size_type p1(edge.first.second());  // index of second endpoint
    // geometry of the edge
    GeometryPtr edge_geo_ptr(std::move(edge.second.first));
    /*
      Build edge here adding another element to the edge vector.
      Use the information in p0, p1 to retrieve pointer to node
    */
    edge_index++;
  }  // end loop over all edges

  const size_type no_of_cells = cells.size();

  // Create an auxiliary data structure to store the edges for the cells
  // For every cell and its edges the global edge index is extracted
  // from the auxiliary associative array.
  std::vector<std::array<size_type, 4>> edge_indices(no_of_cells);

  edge_index = 0;
  for (const EdgeMap::value_type &edge : edge_map) {
    // Obtain array of indices of adjacent cells
    AdjCellsList adjacent_cells(edge.second.second);
    for (const auto &adj_cell : adjacent_cells) {
      const size_type adj_cell_index = adj_cell.first;
      // Local index of edge in adjacent cell
      const size_type edge_local_index = adj_cell.second;
      LF_ASSERT_MSG(adj_cell_index < no_of_cells, "adj_cell_idx out of bounds");
      LF_ASSERT_MSG(edge_local_index < cells[adj_cell_index].first.size(),
                    "local edge index out of bounds");
      // Set edge index in auxiliary cell matrix
      edge_indices[adj_cell_index][edge_local_index] = edge_index;
    }
    edge_index++;
  }

  // Now complete information is available for the construction
  // of cells = entities of co-dimension 0
  /*
     Initialize two vectors, one for trilaterals of size `no_of_trilaterals`
     and a second for quadrilaterals of size `no_of_quadrilaterals`
   */
  // Loop over all cells
  cell_index = 0;
  for (const auto &c : cells) {
    // Node indices for the current cell
    const std::vector<size_type> &c_node_indices(c.first);
    const std::array<size_type, 4> &c_edge_indices(edge_indices[cell_index]);
    const GeometryPtr &c_geo_ptr(c.second);
    if (c_node_indices.size() == 3) {
      /*
        Add a trilateral entity to the vector of trilaterals
        Use information in c_node_indices, c_edge_indices to
        obtain pointers to nodes and edges.
        index = cell_index
       */
    } else {
      /*
         Add a quadrilateral entity to the vector of quadrilaterals
         Use information in c_node_indices, c_edge_indices to
         obtain pointers to nodes and edges.
         index = cell_index
       */
    }
    cell_index++;
  }
}

Mesh::Mesh(char dim_world, std::vector<Eigen::VectorXd> nodes,
           std::vector<std::tuple<std::vector<size_type>,
                                  std::unique_ptr<geometry::Geometry>>>
               elements)
    : dim_world_(dim_world) {
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
  std::vector<std::tuple<std::array<size_type, 2>, size_type, char>>
      element_edges{};
  // Reserve auxiliary vector for edges
  element_edges.reserve(4 * elements.size());
  // Run through all cells, given index arrays of vertices and geometry objects
  for (auto &tuple : elements) {
    // Index array of vertices
    auto &element_nodes = std::get<0>(tuple);
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
    auto &end_edge = element_edges[end];
    if (begin == end) {
      // We have found a new edge.
      std::array<std::vector<size_type>, 1> sub_entities;
      // Adjacent element with the lowest index number
      // sub-entity
      auto &element = entities0_[std::get<1>(element_edges[begin])];
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
}  // namespace lf::mesh::hybrid2d
