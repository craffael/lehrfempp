/**
 * @file
 * @brief Implementation from mesh.h
 * @author Raffael Casagrande
 * @date   2018-06-22 04:31:12
 * @copyright MIT License
 */

#include "mesh.h"
#include <iostream>

namespace lf::mesh::hybrid2d {

//CONTROLDECLARECOMMENT(Mesh, output_ctrl_, "hybrid2d_mesh_output_ctrl",
//                      "Diagnostics control for hybrid2d::Mesh");
ADDOPTION(Mesh::output_ctrl_, hybrid2d_mesh_output_ctrl,
                      "Diagnostics control for hybrid2d::Mesh");

base::ForwardRange<const Entity> Mesh::Entities(char codim) const {
  LF_ASSERT_MSG(codim >= 0, "codim negative.");
  LF_ASSERT_MSG(codim <= dim_world_, "codim > dimWorld.");

  auto l = [&](auto i) -> const mesh::Entity & { return **i; };
  switch (codim) {
    case 0: {
      return {base::make_DereferenceLambdaRandomAccessIterator(
                  entity_pointers_[0].begin(), l),
              base::make_DereferenceLambdaRandomAccessIterator(
                  entity_pointers_[0].end(), l)};
    }
    case 1:
      return {segments_.begin(), segments_.end()};
    case 2:
      return {points_.begin(), points_.end()};
    default: {
      LF_VERIFY_MSG(false, "Something is horribly wrong, codim = " +
                               std::to_string(codim) + " is out of bounds.");
      return {
          base::ForwardIterator<const Entity>(static_cast<Entity *>(nullptr)),
          base::ForwardIterator<const Entity>(static_cast<Entity *>(nullptr))};
    }
  }
}

Mesh::size_type Mesh::NumEntities(char codim) const {
  switch (codim) {
    case 0:
      return trias_.size() + quads_.size();
    case 1:
      return segments_.size();
    case 2:
      return points_.size();
    default:
      LF_VERIFY_MSG(false, "codim out of bounds.");
  }
}

Mesh::size_type Mesh::NumEntities(lf::base::RefEl ref_el_type) const {
  switch (ref_el_type) {
    case lf::base::RefEl::kPoint(): {
      return points_.size();
    }
    case lf::base::RefEl::kSegment(): {
      return segments_.size();
    }
    case lf::base::RefEl::kTria(): {
      return trias_.size();
    }
    case lf::base::RefEl::kQuad(): {
      return quads_.size();
    }
    default: {
      LF_ASSERT_MSG(false, "Illegal entity type");
      break;
    }
  }
  return 0;
}

Mesh::size_type Mesh::Index(const Entity &e) const {
  switch (e.Codim()) {
    case 0: {
      if (e.RefEl() == lf::base::RefEl::kTria()) {
        return dynamic_cast<const Triangle &>(e).index();
      }
      if (e.RefEl() == lf::base::RefEl::kQuad()) {
        return dynamic_cast<const Quadrilateral &>(e).index();
      }
      LF_VERIFY_MSG(false, "Illegal cell type");
    }
    case 1:
      return dynamic_cast<const Segment &>(e).index();
    case 2:
      return dynamic_cast<const Point &>(e).index();
    default:
      LF_VERIFY_MSG(false,
                    "Something is horribyl wrong, this entity has codim = " +
                        std::to_string(e.Codim()));
  }
}

const Entity *Mesh::EntityByIndex(dim_t codim, glb_idx_t index) const {
  LF_ASSERT_MSG(codim <= 2, "Illegal codimension " << codim);
  LF_ASSERT_MSG(index < NumEntities(codim),
                "Index " << index << " > " << NumEntities(codim));
  return entity_pointers_[codim][index];
}

bool Mesh::Contains(const Entity &e) const {
  switch (e.Codim()) {
    case 0:
      return (&e >= &trias_.front() && &e <= &trias_.back()) ||
             (&e >= &quads_.front() && &e <= &quads_.back());
    case 1:
      return &e >= &segments_.front() && &e <= &segments_.back();
    case 2:
      return &e >= &points_.front() && &e <= &points_.back();
    default:
      return false;
  }
}

namespace /*anonymous */ {
/** @brief Auxliary class for mesh_from_node_incidence */
class EndpointIndexPair {
 public:
  using size_type = Mesh::size_type;
  // Constructor
  EndpointIndexPair(size_type p0, size_type p1) : p0_(p0), p1_(p1) {
    LF_ASSERT_MSG(p0 != p1, "No loops allowed");
    if (p1 > p0) {
      cmp_p0_ = p0;
      cmp_p1_ = p1;
    } else {
      cmp_p0_ = p1;
      cmp_p1_ = p0;
    }
  }
  // Access operators
  size_type first_node() const { return p0_; }
  size_type second_node() const { return p1_; }
  // The only comparison operator expected from a Map key
  // Edges are considered equal even if they have the opposite orientation
  friend bool operator<(const EndpointIndexPair &e1,
                        const EndpointIndexPair &e2) {
    return ((e1.cmp_p0_ == e2.cmp_p0_) ? (e1.cmp_p1_ < e2.cmp_p1_)
                                       : (e1.cmp_p0_ < e2.cmp_p0_));
  }
  // Reorienting an edge
  void flip() { std::swap(p0_, p1_); }
  // Checking topological equality (taking into account orientation)
  friend bool coincide(const EndpointIndexPair &e1,
                       const EndpointIndexPair &e2) {
    return ((e1.p0_ == e2.p0_) && (e1.p1_ == e2.p1_));
  }

 private:
  size_type p0_, p1_;  // indices of endpoints
  size_type cmp_p0_, cmp_p1_;
};
}  // namespace

// **********************************************************************
// Construction of a 2D hybrid mesh
//
// Data types for arguments
// GeometryPtr = std::unique_ptr<geometry::Geometry>;
// NodeCoordList = std::vector<GeometryPtr>;
// EdgeList = std::vector<std::pair<std::array<size_type, 2>, GeometryPtr>>;
// CellList = std::vector<std::pair<std::array<size_type, 4>, GeometryPtr>>;
// **********************************************************************
Mesh::Mesh(dim_t dim_world, NodeCoordList nodes, EdgeList edges, CellList cells)
    : dim_world_(dim_world) {
  // Auxiliary data type for gathering information about cells adjacent to an
  // edge
  struct AdjCellInfo {
    AdjCellInfo(size_type _cell_idx, size_type _edge_idx)
        : cell_idx(_cell_idx), edge_idx(_edge_idx) {}
    size_type cell_idx;  // index of adjacent cell
    size_type edge_idx;  // local index of edge in adjacent cell
  };
  // Information about cells adjacent to an edge
  using AdjCellsList = std::vector<AdjCellInfo>;
  // Information about edge in auxiliary array
  struct EdgeData {
    EdgeData(GeometryPtr geo_uptr, AdjCellsList _adj_cells_list)
        : geo_uptr(std::move(geo_uptr)),
          adj_cells_list(std::move(_adj_cells_list)),
          edge_global_index(-1),
          reversed(false) {}
    EdgeData(GeometryPtr geo_uptr, AdjCellsList _adj_cells_list,
             glb_idx_t global_index)
        : geo_uptr(std::move(geo_uptr)),
          adj_cells_list(std::move(_adj_cells_list)),
          edge_global_index(global_index),
          reversed(false) {}
    EdgeData(const EdgeData &) = delete;
    EdgeData(EdgeData &&) = default;
    EdgeData &operator=(const EdgeData &) = delete;
    EdgeData &operator=(EdgeData &&) = default;
    ~EdgeData() = default;
    // temporary storage for unique pointer to edge geometry, if provided
    GeometryPtr geo_uptr;
    // Information about neighbors
    AdjCellsList adj_cells_list;
    // Index of the edge
    glb_idx_t edge_global_index;
    // Flag indicating reversed orientation
    bool reversed;
  };

  // Type of associative auxiliary array for edge information
  using EdgeMap = std::map<EndpointIndexPair, EdgeData>;

  // For extracting point coordinates
  const Eigen::MatrixXd zero_point = base::RefEl::kPoint().NodeCoords();

  // ASSUMPTION: The length of the nodes vector gives the number of nodes

  const size_type no_of_nodes(nodes.size());
  if (output_ctrl_ > 0) {
    std::cout << "Constructing mesh: " << no_of_nodes << " nodes" << std::endl;
  }

  // ======================================================================
  // STEP I: Initialize array of edges using pointers to
  //          entries of the array of nodes

  // Register supplied edges in auxiliary map data structure
  if (output_ctrl_ > 0) {
    std::cout << "Initializing edge map" << std::endl;
  }
  EdgeMap edge_map;
  glb_idx_t edge_index = 0;  // position in the array gives index of edge
  for (auto &e : edges) {
    // Node indices of endpoints: the KEY
    std::array<size_type, 2> end_nodes(e.first);
    EndpointIndexPair e_endpoint_idx(end_nodes[0], end_nodes[1]);
    LF_ASSERT_MSG(
        (end_nodes[0] < no_of_nodes) && (end_nodes[1] < no_of_nodes),
        "Illegal edge node numbers " << end_nodes[0] << ", " << end_nodes[1]);
    if (output_ctrl_ > 0) {
      std::cout << "Register edge: " << end_nodes[0] << " <-> " << end_nodes[1]
                << std::endl;
    }
    // Store provided geometry information; information on adjacent cells
    // not yet available.
    LF_ASSERT_MSG(e.second != nullptr,
                  "Edge " << edge_index << ": missing geometry!");
    AdjCellsList empty_cells_list{};
    EdgeData edge_data(std::move(e.second), empty_cells_list, edge_index);
    EdgeMap::value_type edge_info =
        std::make_pair(e_endpoint_idx, std::move(edge_data));
    std::pair<EdgeMap::iterator, bool> insert_status =
        edge_map.insert(std::move(edge_info));
    LF_ASSERT_MSG(insert_status.second,
                  "Duplicate edge " << end_nodes[0] << " <-> " << end_nodes[1]);

    // If one of the endpoints of a edge does not have a geometry, supply it
    // with one inherited from the edge.
    for (int j = 0; j < 2; ++j) {
      if (nodes[end_nodes[j]] == nullptr) {
        // if no geometry for node exists request geomtry for an endpoint of the
        // edge Note: endpoints are entities of relative co-dimension 1
        nodes[end_nodes[j]] = std::move(e.second->SubGeometry(1, j));
      }
    }

    edge_index++;
  }  // end loop over predefined edges
  // ======================================================================

  // ======================================================================
  // Step II: Building edge map from cell information
  //
  // At this point all predefined edges have been stored in the auxiliary
  // associative array, though without information about adjacent cells.
  // The variable edge_index contains the number of edges with externally
  // supplied geometry. The indexing of all extra edges created below must start
  // from this offset.

  if (output_ctrl_ > 10) {
    std::cout << "Edge map after edge registration" << std::endl;
    for (auto &edge_info : edge_map) {
      const EndpointIndexPair &eip(edge_info.first);
      const EdgeData &edat(edge_info.second);
      std::cout << "Edge " << eip.first_node() << " <-> " << eip.second_node()
                << ": ";
      const AdjCellsList &acl(edat.adj_cells_list);
      const GeometryPtr &gptr(edat.geo_uptr);
      for (auto &i : acl) {
        std::cout << "[" << i.cell_idx << "," << i.edge_idx << "] ";
      }
      std::cout << std::endl;
    }
  }

  // Run through cells in order to
  // (i) build edges missing in the list of predefined edges
  // (ii) determine cells adjacent to edges
  size_type cell_index = 0;
  size_type no_of_trilaterals = 0;
  size_type no_of_quadrilaterals = 0;
  // Diagnostics
  if (output_ctrl_ > 0) {
    std::cout << "Scanning list of cells" << std::endl;
  }
  for (const auto &c : cells) {
    // node indices of corners of cell c
    const std::array<size_type, 4> &cell_node_list(c.first);
    // Geometry of current cell
    const GeometryPtr &cell_geometry(c.second);
    // Can be either a trilateral or a quadrilateral
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // A triangle is marked by an invalid node number
    // in the last position
    size_type no_of_vertices;
    if (cell_node_list[3] == idx_nil) {
      no_of_vertices = 3;  // triangle
      no_of_trilaterals++;
    } else {
      no_of_vertices = 4;  // quadrilateral
      no_of_quadrilaterals++;
    }
    // Fix the type of the cell
    base::RefEl ref_el =
        (no_of_vertices == 3) ? base::RefEl::kTria() : base::RefEl::kQuad();
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (output_ctrl_ > 10) {
      std::cout << "Cell " << cell_index;
      if (no_of_vertices == 3) {
        std::cout << ", tria " << no_of_trilaterals << ": ";
      } else {
        std::cout << ", quad " << no_of_quadrilaterals << ": ";
      }
    }

    // Verify validity of vertex indices
    for (int l = 0; l < no_of_vertices; l++) {
      LF_VERIFY_MSG(cell_node_list[l] < no_of_nodes,
                    "Node " << l << " of cell " << cell_index
                            << ": invalid index " << cell_node_list[l]);
    }

    // If the cell has a geometry, it can be used to generate the
    // geometry for its vertices through the SubGeometry() method of Geometry
    // objects
    if (cell_geometry != nullptr) {
      for (int j = 0; j < no_of_vertices; ++j) {
        if (nodes[cell_node_list[j]] == nullptr) {
          // if no geometry for node exists request geomtry for an vertex from
          // the cell Note: vertices are entities of relative co-dimension 2
          nodes[cell_node_list[j]] =
              std::move(cell_geometry->SubGeometry(2, j));
        }
      }
    }

    //  A variant:
    //    There may be cells without a specified geometry.
    //    In case an edge is not equipped with a geometry and not
    //    adjacent to a cell with a specific geometry, assume a straight
    //    edge connecting the two endpoints.
    //

    // Visit all edges of the current cell
    for (unsigned int j = 0; j < ref_el.NumSubEntities(1); j++) {
      // Fetch local indices of endpoints of edge j
      const size_type p0_local_index =
          ref_el.SubSubEntity2SubEntity(1, j, 1, 0);
      const size_type p1_local_index =
          ref_el.SubSubEntity2SubEntity(1, j, 1, 1);
      // Fetch global indices of the endnodes of edge j
      EndpointIndexPair c_edge_vertex_indices(cell_node_list[p0_local_index],
                                              cell_node_list[p1_local_index]);

      if (output_ctrl_ > 10) {
        std::cout << "e(" << j << ") = local " << p0_local_index << " <-> "
                  << p1_local_index << ", global "
                  << c_edge_vertex_indices.first_node() << " <-> "
                  << c_edge_vertex_indices.second_node() << " # ";
      }
      // Store number of cell and the local index j of the edge
      AdjCellInfo edge_cell_info(cell_index, j);
      // Check whether edge exists already
      auto edge_ptr = edge_map.find(c_edge_vertex_indices);
      if (edge_ptr == edge_map.end()) {
        // Inherit geometry of edge from adjacent cell
        GeometryPtr edge_geo_ptr;
        if (cell_geometry) {
          edge_geo_ptr = std::move(cell_geometry->SubGeometry(1, j));
        }
        // Beginning of list of adjacent elements
        AdjCellsList single_cell_list{edge_cell_info};
        EdgeData edge_data(std::move(edge_geo_ptr), single_cell_list);
        std::pair<EdgeMap::iterator, bool> insert_status = edge_map.insert(
            std::make_pair(c_edge_vertex_indices, std::move(edge_data)));
        LF_ASSERT_MSG(insert_status.second, "Duplicate not found earlier!");
        edge_ptr = insert_status.first;  // pointer to newly inserted edge
      } else {
        // Here *edge_ptr is a valid std::pair<EndpointIndexPair,EdgeData>

        // Store information about neighboring cell
        AdjCellsList &edge_adj_cellslist((edge_ptr->second).adj_cells_list);
        edge_adj_cellslist.push_back(edge_cell_info);
        // Check, if the edge possesses geometry information
        GeometryPtr &edge_supplied_geo_uptr((*edge_ptr).second.geo_uptr);
        if (edge_supplied_geo_uptr == nullptr) {
          // Edge does not know its geometry yet. Try to obtain it from the
          // cell.
          if (cell_geometry) {
            edge_supplied_geo_uptr =
                std::move(cell_geometry->SubGeometry(1, j));
            // NOTE: the local orientation of the edge of the cell and that
            // pointed to by edge_ptr can differ. In this case the
            // endpoints of the edge have to be swapped. To detect orientation
            // mismatch, compare c_edge_vertex_indices (cell infomation) and
            // edge_ptr->first (edge information)
            if (!coincide(c_edge_vertex_indices, edge_ptr->first)) {
              // Reverse the direction of the edge
              (edge_ptr->second).reversed = true;
            }
          }
        }
      }
      // At this point edge_ptr points to the map entry for the current edge
    }  // end of loop over edges
    cell_index++;

    // Diagnostics
    if (output_ctrl_ > 10) {
      std::cout << std::endl;
    }
  }  // end loop over cells

  // DIAGNOSTICS
  {
    if (output_ctrl_ > 0) {
      std::cout << "=============================================" << std::endl;
    }
    if (output_ctrl_ > 10) {
      std::cout << "Edge map after cell scan" << std::endl;
      size_type edge_cnt = 0;
      for (auto &edge_info : edge_map) {
        const EndpointIndexPair &eip(edge_info.first);
        const EdgeData &edat(edge_info.second);
        std::cout << "Edge " << edge_cnt << ": " << eip.first_node() << " <-> "
                  << eip.second_node();
        if (edat.edge_global_index == -1) {
          std::cout << " no index : ";
        } else {
          std::cout << ": index = " << edat.edge_global_index << ": ";
        }
        const AdjCellsList &acl(edat.adj_cells_list);
        const GeometryPtr &gptr(edat.geo_uptr);
        for (auto &i : acl) {
          std::cout << "[" << i.cell_idx << "," << i.edge_idx << "] ";
        }
        std::cout << " geo = " << std::endl;
        if (gptr) {
          Eigen::MatrixXd edp_c(
              gptr->Global(base::RefEl::kSegment().NodeCoords()));
          std::cout << edp_c << std::endl;
        } else {
          std::cout << "NO GEOMETRY" << std::endl;
        }
        edge_cnt++;
      }
      std::cout << "=============================================" << std::endl;
    }  // end if(output_ctrl > 10)
  }

  // ======================================================================
  // NEXT STEP : Set up and fill array of nodes: points_
  // In the beginning initialize vector of vertices and do not touch it anymore
  // Initialize vector for Node entities of size `no_of_nodes`
  points_.reserve(no_of_nodes);

  size_type node_index = 0;
  for (hybrid2d::Mesh::GeometryPtr &pt_geo_ptr : nodes) {
    // OLD VERSION. In the new version the geomtry of a point is passed in
    // 'nodes' const Eigen::VectorXd &node_coordinates(v); GeometryPtr point_geo
    // = std::make_unique<geometry::Point>(node_coordinates);
    LF_VERIFY_MSG(pt_geo_ptr != nullptr,
                  "Missing geometry for node " << node_index);
    if (output_ctrl_ > 10) {
      std::cout
          << "-> Adding node " << node_index << " at "
          << (pt_geo_ptr->Global(Eigen::Matrix<double, 0, 1>())).transpose()
          << std::endl;
    }
    points_.emplace_back(node_index, std::move(pt_geo_ptr));
    node_index++;
  }

  // Run through the entire associative container for edges
  // and build edge Entities
  // This is the length to be reserved for the edge vector
  size_type no_of_edges = edge_map.size();

  // Initialized vector of Edge entities here
  segments_.reserve(no_of_edges);

  // Note: the variable edge_index contains the number externally supplied
  // edges, whose index must agree with their position in the
  // 'edges' array.

  // Note: iteration variable cannot be declared const, because
  // moving the geometry pointer changes it!
  for (EdgeMap::value_type &edge : edge_map) {
    // Indices of the two endpoints of the current edge
    // Use this to obtain pointers/references to nodes
    // from the vector of nodes
    size_type p0(edge.first.first_node());   // index of first endpoint
    size_type p1(edge.first.second_node());  // index of second endpoint

    // Swap the indices of the endpoints in case of a geometry mismatch
    // that requires a reversal of orientation.
    if (edge.second.reversed) {
      std::swap(p0, p1);
    }

    const Point *p0_ptr = &points_[p0];  // pointer to first endpoint
    const Point *p1_ptr = &points_[p1];  // pointer to second endpoint

    // obtain geometry of the edge
    GeometryPtr edge_geo_ptr(std::move(edge.second.geo_uptr));
    if (!edge_geo_ptr) {
      // If the edge does not have a geometry build a straight edge
      Eigen::Matrix<double, 2, 2> straight_edge_coords;
      straight_edge_coords.block<2, 1>(0, 0) =
          p0_ptr->Geometry()->Global(zero_point);
      straight_edge_coords.block<2, 1>(0, 1) =
          p1_ptr->Geometry()->Global(zero_point);
      edge_geo_ptr =
          std::make_unique<geometry::SegmentO1>(straight_edge_coords);
    }
    // Determine index of current edge. Two cases have to be distinguished:
    // (i) the edge geometry was specified in the `edges` argument. In this case
    // the edge index must agree with its possition in that array. This position
    /// is stored in the 'edge_global_index' field of the EdgeData structure.
    // (ii) the edge has to be created internally. In this case assign an index
    // larger than the index of any supplied edge.
    if (edge.second.edge_global_index == idx_nil) {
      // Internally created edge needs new index.
      edge.second.edge_global_index = edge_index;
      // Increment 'edge_index', which will give the index of the next
      // internally created edge
      edge_index++;
    }

    // Diagnostics
    if (output_ctrl_ > 10) {
      std::cout << "Registering edge " << edge.second.edge_global_index << ": "
                << p0 << " <-> " << p1 << std::endl;
    }
    // Building edge by adding another element to the edge vector.
    segments_.emplace_back(edge.second.edge_global_index,
                           std::move(edge_geo_ptr), p0_ptr, p1_ptr);
  }  // end loop over all edges
  LF_ASSERT_MSG(edge_index == no_of_edges, "Edge index mismatch");

  // ======================================================================
  // NEXT STEP: Create cells

  const size_type no_of_cells = cells.size();

  // Create an auxiliary data structure to store the edges for the cells
  // For every cell and its edges the global edge index is extracted
  // from the auxiliary associative array.
  std::vector<std::array<size_type, 4>> edge_indices(no_of_cells);

  size_type edge_array_position = 0;  // actual position in edge array
  for (const EdgeMap::value_type &edge : edge_map) {
    // Obtain array of indices of adjacent cells
    AdjCellsList adjacent_cells(edge.second.adj_cells_list);
    for (const auto &adj_cell : adjacent_cells) {
      const size_type adj_cell_index = adj_cell.cell_idx;
      // Local index of edge in adjacent cell
      const size_type edge_local_index = adj_cell.edge_idx;
      LF_ASSERT_MSG(adj_cell_index < no_of_cells, "adj_cell_idx out of bounds");
      LF_ASSERT_MSG(edge_local_index < cells[adj_cell_index].first.size(),
                    "local edge index out of bounds");
      // Set edge index in auxiliary cell matrix
      edge_indices[adj_cell_index][edge_local_index] = edge_array_position;
    }
    edge_array_position++;
  }

  // Diagnostics
  if (output_ctrl_ > 10) {
    std::cout << "########################################" << std::endl;
    std::cout << "Edge array indices for cells " << std::endl;
  }

  // Now complete information is available for the construction
  // of cells = entities of co-dimension 0
  //     Initialize two vectors, one for trilaterals of size `no_of_trilaterals`
  //   and a second for quadrilaterals of size `no_of_quadrilaterals`
  trias_.reserve(no_of_trilaterals);
  quads_.reserve(no_of_quadrilaterals);
  // Loop over all cells
  cell_index = 0;
  for (auto &c : cells) {
    // Node indices for the current cell
    const std::array<size_type, 4> &c_node_indices(c.first);
    const std::array<size_type, 4> &c_edge_indices(edge_indices[cell_index]);
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // A triangle is marked by an invalid node number
    // in the last position (also see above)
    size_type no_of_vertices;
    if (c_node_indices[3] == size_type(-1)) {
      no_of_vertices = 3;  // triangle
    } else {
      no_of_vertices = 4;  // quadrilateral
    }

    // Verify validity of node/ede indices
    for (int l = 0; l < no_of_vertices; l++) {
      LF_VERIFY_MSG(c_node_indices[l] < no_of_nodes,
                    "Node " << l << " of cell " << cell_index
                            << ": invalid index " << c_node_indices[l]);
      LF_VERIFY_MSG(c_edge_indices[l] < no_of_edges,
                    "Edge " << l << " of cell " << cell_index
                            << ": invalid index " << c_edge_indices[l]);
    }

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    GeometryPtr c_geo_ptr(std::move(c.second));
    if (no_of_vertices == 3) {
      // Case of a trilateral

      // Diagnostics
      if (output_ctrl_ > 10) {
        std::cout << "Triangular cell " << cell_index << ": nodes "
                  << c_node_indices[0] << ", " << c_node_indices[1] << ", "
                  << c_node_indices[2] << ", edges " << c_edge_indices[0]
                  << ", " << c_edge_indices[1] << ", " << c_edge_indices[2]
                  << std::endl;
      }
      /*
        Add a trilateral entity to the vector of trilaterals
        Use information in c_node_indices, c_edge_indices to
        obtain pointers to nodes and edges.
        index = cell_index
      */
      const Point *corner0 = &points_[c_node_indices[0]];
      const Point *corner1 = &points_[c_node_indices[1]];
      const Point *corner2 = &points_[c_node_indices[2]];
      const Segment *edge0 = &segments_[c_edge_indices[0]];
      const Segment *edge1 = &segments_[c_edge_indices[1]];
      const Segment *edge2 = &segments_[c_edge_indices[2]];
      if (!c_geo_ptr) {
        // Cell is lacking a geometry and its shape has to
        // be determined from the shape of the edges or
        // location of the vertices
        // At this point only the latter policy is implemented
        // and we build an affine triangle.
        // First assemble corner coordinates into matrix
        Eigen::Matrix<double, 2, 3> triag_corner_coords;
        triag_corner_coords.block<2, 1>(0, 0) =
            corner0->Geometry()->Global(zero_point);
        triag_corner_coords.block<2, 1>(0, 1) =
            corner1->Geometry()->Global(zero_point);
        triag_corner_coords.block<2, 1>(0, 2) =
            corner2->Geometry()->Global(zero_point);

        // Diagnostics
        if (output_ctrl_ > 10) {
          std::cout << "Creating triangle with geometry " << std::endl
                    << triag_corner_coords << std::endl;
        }

        // Then create geometry of an affine triangle
        c_geo_ptr = std::make_unique<geometry::TriaO1>(triag_corner_coords);
        // For later:
        // If blended geometries are available, a cell could also
        // inherit its geometry from the edges
      }
      trias_.emplace_back(cell_index, std::move(c_geo_ptr), corner0, corner1,
                          corner2, edge0, edge1, edge2);
    } else {
      // Case of a quadrilateral

      // Diagnostics
      if (output_ctrl_ > 10) {
        std::cout << "Quadrilateral cell " << cell_index << ": nodes "
                  << c_node_indices[0] << ", " << c_node_indices[1] << ", "
                  << c_node_indices[2] << ", " << c_node_indices[3]
                  << ", edges " << c_edge_indices[0] << ", "
                  << c_edge_indices[1] << ", " << c_edge_indices[2] << ", "
                  << c_edge_indices[3] << std::endl;
      }
      /*
        Add a quadrilateral entity to the vector of quadrilaterals
        Use information in c_node_indices, c_edge_indices to
        obtain pointers to nodes and edges.
        index = cell_index
      */
      const Point *corner0 = &points_[c_node_indices[0]];
      const Point *corner1 = &points_[c_node_indices[1]];
      const Point *corner2 = &points_[c_node_indices[2]];
      const Point *corner3 = &points_[c_node_indices[3]];
      const Segment *edge0 = &segments_[c_edge_indices[0]];
      const Segment *edge1 = &segments_[c_edge_indices[1]];
      const Segment *edge2 = &segments_[c_edge_indices[2]];
      const Segment *edge3 = &segments_[c_edge_indices[3]];
      if (!c_geo_ptr) {
        // Cell is lacking a geometry and its shape has to
        // be determined from the shape of the edges or
        // location of the vertices
        // At this point we can only build a quadrilateral
        // with straight edges ("bilinear quadrilateral")
        // First assemble corner coordinates into matrix

        Eigen::Matrix<double, 2, 4> quad_corner_coords;
        quad_corner_coords.block<2, 1>(0, 0) =
            corner0->Geometry()->Global(zero_point);
        quad_corner_coords.block<2, 1>(0, 1) =
            corner1->Geometry()->Global(zero_point);
        quad_corner_coords.block<2, 1>(0, 2) =
            corner2->Geometry()->Global(zero_point);
        quad_corner_coords.block<2, 1>(0, 3) =
            corner3->Geometry()->Global(zero_point);
        // Then create geometry of an affine triangle

        // Diagnostics
        if (output_ctrl_ > 0) {
          std::cout << "Creating quadrilateral with geometry " << std::endl
                    << quad_corner_coords << std::endl;
        }

        c_geo_ptr = std::make_unique<geometry::QuadO1>(quad_corner_coords);
      }
      quads_.emplace_back(cell_index, std::move(c_geo_ptr), corner0, corner1,
                          corner2, corner3, edge0, edge1, edge2, edge3);
    }
    cell_index++;
  }
  // ======================================================================
  // Initialization of auxiliary entity pointer arrays

  // Order the cells according to their indices!
  // First  fill the array with NIL pointers
  entity_pointers_[0] = std::move(std::vector<const mesh::Entity *>(
      trias_.size() + quads_.size(), nullptr));
  size_type cell_ptr_cnt = 0;

  // First retrieve pointers to triangular cells
  for (int j = 0; j < trias_.size(); j++, cell_ptr_cnt++) {
    // Fetch index of a triangle  by a non-virtual call to Index()
    const glb_idx_t cell_index = lf::mesh::hybrid2d::Mesh::Index(trias_[j]);
    LF_ASSERT_MSG(cell_index < entity_pointers_[0].size(),
                  "Cell index out of range");
    // This index must be unique !
    LF_ASSERT_MSG(entity_pointers_[0][cell_index] == nullptr,
                  "Cell index " << cell_index << " occurs twice!");
    entity_pointers_[0][cell_index] = &(trias_[j]);
  }  // end loop over triangles

  // Second deal with the quadrilateral cells
  for (int j = 0; j < quads_.size(); j++, cell_ptr_cnt++) {
    // Fetch index of a quadrilateral by a non-virtual call to Index()
    const glb_idx_t cell_index = lf::mesh::hybrid2d::Mesh::Index(quads_[j]);
    LF_ASSERT_MSG(cell_index < entity_pointers_[0].size(),
                  "Cell index out of range");
    // This index must be unique !
    LF_ASSERT_MSG(entity_pointers_[0][cell_index] == nullptr,
                  "Cell index " << cell_index << " occurs twice!");
    entity_pointers_[0][cell_index] = &(quads_[j]);
  }
  LF_VERIFY_MSG(cell_ptr_cnt == entity_pointers_[0].size(),
                "Size mismatch for cell counters array");

  // Next the edges and nodes

  entity_pointers_[1] =
      std::move(std::vector<const mesh::Entity *>(segments_.size(), nullptr));
  entity_pointers_[2] =
      std::move(std::vector<const mesh::Entity *>(points_.size(), nullptr));

  for (dim_t codim = 1; codim <= 2; codim++) {
    // Loop over all entities of codimenson codim
    for (const Entity &e : lf::mesh::hybrid2d::Mesh::Entities(codim)) {
      const glb_idx_t entity_index = lf::mesh::hybrid2d::Mesh::Index(e);
      LF_ASSERT_MSG(entity_index < entity_pointers_[codim].size(),
                    "Entity(" << codim << ") index out of range");
      // This index must be unique !
      LF_ASSERT_MSG(
          entity_pointers_[codim][entity_index] == nullptr,
          "Entity(" << codim << ") index " << entity_index << " occurs twice!");
      entity_pointers_[codim][entity_index] =
          static_cast<const mesh::Entity *>(&e);
    }  // loop over entities
  }    // loop over codim's

}  // end of constructor

}  // namespace lf::mesh::hybrid2d
