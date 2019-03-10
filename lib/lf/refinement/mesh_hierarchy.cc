/**
 * @file refinement_hierarchy.cc
 * @brief implementation of global/local refinement methods
 */

#include "mesh_hierarchy.h"
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/utils/utils.h"

namespace lf::refinement {

// Debugging function: check validity of index vectors
bool checkValidIndex(const std::vector<glb_idx_t> &idx_vec) {
  const size_type n_idx = idx_vec.size();
  for (int n = 0; n < n_idx; n++) {
    if (idx_vec[n] == idx_nil) {
      return false;
    }
  }
  return true;
}

CONTROLDECLARECOMMENT(MeshHierarchy, output_ctrl_, "MeshHierarchy_output_ctrl",
                      "Diagnostics control for MeshHierarchy");

CONTROLDECLARECOMMENT(MeshHierarchy, ctrl_, "MeshHierarchy_ctrl",
                      "Output control for MeshHierarchy");

// Implementation of MeshHierarchy
MeshHierarchy::MeshHierarchy(std::shared_ptr<mesh::Mesh> base_mesh,  // NOLINT
                             std::shared_ptr<mesh::MeshFactory> mesh_factory)
    : mesh_factory_(std::move(mesh_factory)) {
  LF_VERIFY_MSG(base_mesh, "No valid mesh supplied");
  LF_VERIFY_MSG(base_mesh->DimMesh() == 2, "Implemented only for 2D meshes");
  // Set coarsest mesh
  meshes_.push_back(base_mesh);
  {
    // Refinement edge has to be set for edges
    refinement_edges_.emplace_back(base_mesh->NumEntities(0), idx_nil);
    for (const mesh::Entity &cell : base_mesh->Entities(0)) {
      lf::base::glb_idx_t cell_index = base_mesh->Index(cell);
      if (cell.RefEl() == lf::base::RefEl::kTria()) {
        (refinement_edges_.back())[cell_index] = LongestEdge(cell);
      }
    }  // end loop over cells

    // Setting up child information
    std::vector<CellChildInfo> cell_child_info(base_mesh->NumEntities(0));
    std::vector<EdgeChildInfo> edge_child_info(base_mesh->NumEntities(1));
    std::vector<PointChildInfo> point_child_info(base_mesh->NumEntities(2));
    cell_child_infos_.push_back(std::move(cell_child_info));
    edge_child_infos_.push_back(std::move(edge_child_info));
    point_child_infos_.push_back(std::move(point_child_info));
  }
  {
    // No parents for entities on the coarsest level
    std::vector<ParentInfo> cell_parent_info(base_mesh->NumEntities(0));
    std::vector<ParentInfo> edge_parent_info(base_mesh->NumEntities(1));
    std::vector<ParentInfo> point_parent_info(base_mesh->NumEntities(2));
    parent_infos_.push_back({std::move(cell_parent_info),
                             std::move(edge_parent_info),
                             std::move(point_parent_info)});
  }
  edge_marked_.push_back(
      std::move(std::vector<bool>(base_mesh->NumEntities(1), false)));
}

void MeshHierarchy::RefineRegular(RefPat ref_pat) {
  LF_VERIFY_MSG(
      ref_pat == RefPat::rp_regular || ref_pat == RefPat::rp_barycentric,
      "Only regular or barycentric uniform refinement possible");
  // Retrieve the finest mesh in the hierarchy
  const mesh::Mesh &finest_mesh(*meshes_.back());
  // Arrays containing refinement information for finest mesh
  std::vector<PointChildInfo> &finest_point_ci(point_child_infos_.back());
  std::vector<EdgeChildInfo> &finest_edge_ci(edge_child_infos_.back());
  std::vector<CellChildInfo> &finest_cell_ci(cell_child_infos_.back());

  LF_VERIFY_MSG(finest_point_ci.size() == finest_mesh.NumEntities(2),
                "length mismatch PointChildInfo vector");
  LF_VERIFY_MSG(finest_edge_ci.size() == finest_mesh.NumEntities(1),
                "length mismatch EdgeChildInfo vector");
  LF_VERIFY_MSG(finest_cell_ci.size() == finest_mesh.NumEntities(0),
                "length mismatch CellChildInfo vector");

  // Flag all points as to be copied
  for (const mesh::Entity &point : finest_mesh.Entities(2)) {
    const lf::base::glb_idx_t point_index = finest_mesh.Index(point);
    PointChildInfo &pt_child_info(finest_point_ci[point_index]);
    // Set the information about a point's children except the child pointer
    pt_child_info.ref_pat = RefPat::rp_copy;
  }
  // Flag all edges as to be split
  for (const mesh::Entity &edge : finest_mesh.Entities(1)) {
    const lf::base::glb_idx_t edge_index = finest_mesh.Index(edge);
    EdgeChildInfo &ed_child_info(finest_edge_ci[edge_index]);
    ed_child_info.ref_pat_ = RefPat::rp_split;
  }
  // All cells are to be refined regularly
  for (const mesh::Entity &cell : finest_mesh.Entities(0)) {
    const lf::base::glb_idx_t cell_index = finest_mesh.Index(cell);
    CellChildInfo &cell_child_info(finest_cell_ci[cell_index]);
    cell_child_info.ref_pat_ = ref_pat;
  }
  // With all refinement patterns set, generate the new mesh
  PerformRefinement();
}

void MeshHierarchy::RefineMarked() {
  // Target the finest mesh
  const lf::mesh::Mesh &finest_mesh(*meshes_.back());
  // Arrays containing refinement information for finest mesh
  std::vector<PointChildInfo> &finest_point_ci(point_child_infos_.back());
  std::vector<EdgeChildInfo> &finest_edge_ci(edge_child_infos_.back());
  std::vector<CellChildInfo> &finest_cell_ci(cell_child_infos_.back());

  LF_VERIFY_MSG(finest_point_ci.size() == finest_mesh.NumEntities(2),
                "length mismatch PointChildInfo vector");
  LF_VERIFY_MSG(finest_edge_ci.size() == finest_mesh.NumEntities(1),
                "length mismatch EdgeChildInfo vector");
  LF_VERIFY_MSG(finest_cell_ci.size() == finest_mesh.NumEntities(0),
                "length mismatch CellChildInfo vector");

  // Flag all points as to be copied
  for (const mesh::Entity &point : finest_mesh.Entities(2)) {
    const glb_idx_t point_index = finest_mesh.Index(point);
    PointChildInfo &pt_child_info(finest_point_ci[point_index]);
    // Set the information about a points children except the child pointer
    pt_child_info.ref_pat = RefPat::rp_copy;
  }

  // Set the split refinement patterrn for all marked edges
  for (const lf::mesh::Entity &edge : finest_mesh.Entities(1)) {
    LF_VERIFY_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
                  "Wrong type for an edge");
    const glb_idx_t edge_index = finest_mesh.Index(edge);
    EdgeChildInfo &ed_ci(finest_edge_ci[edge_index]);
    LF_VERIFY_MSG(ed_ci.ref_pat_ == RefPat::rp_nil,
                  "Edge " << edge_index << " already refined!");
    if ((edge_marked_.back())[edge_index]) {
      // Edge is to be refined
      ed_ci.ref_pat_ = RefPat::rp_split;
    } else {
      // Just copy edge
      ed_ci.ref_pat_ = RefPat::rp_copy;
    }
  }
  // Now all edges are initially marked to be split or copied

  // To keep the mesh conforming refinement might have to propagate
  // This is achieved in the following REPEAT ... UNTIL loop.
  bool refinement_complete;
  do {
    refinement_complete = true;
    // Visit all cells and update their refinement patterns
    for (const lf::mesh::Entity &cell : finest_mesh.Entities(0)) {
      const glb_idx_t cell_index = finest_mesh.Index(cell);
      CellChildInfo &cell_child_info(finest_cell_ci[cell_index]);

      // Global indices of edges
      std::array<glb_idx_t, 4> cell_edge_indices{};

      // Find edges which are marked as split
      std::array<bool, 4> edge_split{{false, false, false, false}};
      // Local indices of edges marked as split
      std::array<sub_idx_t, 4> split_edge_idx{};
      // Array of references to edge sub-entities of current cell
      base::RandomAccessRange<const lf::mesh::Entity> sub_edges(
          cell.SubEntities(1));
      const size_type num_edges = cell.RefEl().NumSubEntities(1);
      LF_VERIFY_MSG(num_edges <= 4, "Too many edges = " << num_edges);
      // Obtain information about current splitting pattern of
      // the edges of the cell
      size_type split_edge_cnt = 0;
      for (int k = 0; k < num_edges; k++) {
        const glb_idx_t edge_index = finest_mesh.Index(sub_edges[k]);
        cell_edge_indices[k] = edge_index;
        edge_split[k] =
            (finest_edge_ci[edge_index].ref_pat_ == RefPat::rp_split);
        if (edge_split[k]) {
          split_edge_idx[split_edge_cnt] = k;
          split_edge_cnt++;
        }
      }
      switch (cell.RefEl()) {
        case lf::base::RefEl::kTria(): {
          // Case of a triangular cell: In this case bisection refinement
          // is performed starting with the refinement edge.
          // Local index of refinement edge for the current triangle, also
          // called the "anchor edge" in the case of repeated  bisection
          const sub_idx_t anchor = refinement_edges_.back()[cell_index];
          LF_VERIFY_MSG(anchor < 3, "Illegal anchor = " << anchor);

          // Refinement edge will always be the anchor edge
          finest_cell_ci[cell_index].anchor_ = anchor;
          const sub_idx_t mod_0 = anchor;
          const sub_idx_t mod_1 = (anchor + 1) % 3;
          const sub_idx_t mod_2 = (anchor + 2) % 3;
          // Flag tuple indicating splitting status of an edge: true <-> split
          std::tuple<bool, bool, bool> split_status(
              {edge_split[mod_0], edge_split[mod_1], edge_split[mod_2]});

          // Determine updated refinement pattern for triangle depending on the
          // splitting status of its edges. If the triangle is subdivided, the
          // refinement must always be split in a first bisection step, even if
          // it may not have been marked as split. In this case refinement may
          // spread to neighboring cells and we have to iterative once more.
          // Therefore the refinement_complete flag has to be set to 'false'.

          if (split_status ==
              std::tuple<bool, bool, bool>({false, false, false})) {
            // No edge to be split: just copy triangle
            LF_VERIFY_MSG(split_edge_cnt == 0, "Wrong number of split edges");
            finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_copy;
          } else if (split_status ==
                     std::tuple<bool, bool, bool>({true, false, false})) {
            // Only refinement edge has to be split by a single bisection
            // No additional edge will be split
            LF_VERIFY_MSG(split_edge_cnt == 1, "Wrong number of split edges");
            finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_bisect;
          } else if (split_status ==
                     std::tuple<bool, bool, bool>({true, true, false})) {
            // Trisection refinement, no extra splitting of edges
            LF_VERIFY_MSG(split_edge_cnt == 2, "Wrong number of split edges");
            finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_trisect;
          } else if (split_status ==
                     std::tuple<bool, bool, bool>({false, true, false})) {
            // Trisection refinement, triggering splitting of refinement edge
            LF_VERIFY_MSG(split_edge_cnt == 1, "Wrong number of split edges");
            finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_trisect;
            finest_edge_ci[cell_edge_indices[anchor]].ref_pat_ =
                RefPat::rp_split;
            refinement_complete = false;
          } else if (split_status ==
                     std::tuple<bool, bool, bool>({true, false, true})) {
            // Trisection refinement (other side), no extra splitting of edges
            LF_VERIFY_MSG(split_edge_cnt == 2, "Wrong number of split edges");
            finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_trisect_left;
          } else if (split_status ==
                     std::tuple<bool, bool, bool>({false, false, true})) {
            // Trisection refinement (other side), triggering splitting of
            // refinement edge
            LF_VERIFY_MSG(split_edge_cnt == 1, "Wrong number of split edges");
            finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_trisect_left;
            finest_edge_ci[cell_edge_indices[anchor]].ref_pat_ =
                RefPat::rp_split;
            refinement_complete = false;
          } else if (split_status ==
                     std::tuple<bool, bool, bool>({true, true, true})) {
            // Quadsection refinement, no extra splitting of edges
            LF_VERIFY_MSG(split_edge_cnt == 3, "Wrong number of split edges");
            finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_quadsect;
          } else if (split_status ==
                     std::tuple<bool, bool, bool>({false, true, true})) {
            // Quadsection refinement requiring splitting of refinement edge
            LF_VERIFY_MSG(split_edge_cnt == 2, "Wrong number of split edges");
            finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_quadsect;
            finest_edge_ci[cell_edge_indices[anchor]].ref_pat_ =
                RefPat::rp_split;
            refinement_complete = false;
          } else {
            LF_VERIFY_MSG(false, "Impossible case");
          }
          break;
        }  // end case of a triangle
        case lf::base::RefEl::kQuad(): {
          // There is no refinement edge for quadrilaterals and so no extra edge
          // splitting will be necessary. The refinement pattern for a
          // quadrilateral will be determined from the number of edges split and
          // their location to each other.
          switch (split_edge_cnt) {
            case 0: {
              // No edge split: quadrilateral has to be copied
              finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_copy;
              break;
            }
            case 1: {
              // One edge split: trisection refinement of the quadrilateral
              // Anchor edge is the split edge
              finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_trisect;
              finest_cell_ci[cell_index].anchor_ = split_edge_idx[0];
              break;
            }
            case 2: {
              if ((split_edge_idx[1] - split_edge_idx[0]) == 2) {
                // If the two split edges are opposite to each other, then
                // bisection of the quadrilateral is the right refinement
                // pattern.
                finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_bisect;
                finest_cell_ci[cell_index].anchor_ = split_edge_idx[0];
              } else {
                // Tthe two split edges are adjacent, this case can be
                // accommodated by quadsection refinement. Anchor is the split
                // edge with the lower index (modulo 4).
                finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_quadsect;
                if (((split_edge_idx[0] + 1) % 4) == split_edge_idx[1]) {
                  finest_cell_ci[cell_index].anchor_ = split_edge_idx[0];
                } else if (((split_edge_idx[1] + 1) % 4) == split_edge_idx[0]) {
                  finest_cell_ci[cell_index].anchor_ = split_edge_idx[1];
                } else {
                  LF_VERIFY_MSG(false,
                                "Quad: impossible situation for 2 split edges");
                }
              }
              break;
            }
            case 3: {
              // Three edges of the quadrilateral are split, which can be
              // accommodated only by the rp_threeedge refinement pattern
              // anchor is the edge with the middle index
              finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_threeedge;
              if (!edge_split[0]) {  // Split edges 1,2,3, middle edge 2
                finest_cell_ci[cell_index].anchor_ = 2;
              } else if (!edge_split[1]) {  // Split edges 0,2,3, middle edge 3
                finest_cell_ci[cell_index].anchor_ = 3;
              } else if (!edge_split[2]) {  // Split edges 0,1,3, middle edge 0
                finest_cell_ci[cell_index].anchor_ = 0;
              } else if (!edge_split[3]) {  // Split edges 0,1,2, middle edge 1
                finest_cell_ci[cell_index].anchor_ = 1;
              } else {
                LF_VERIFY_MSG(false, "Inconsistent split pattern");
              }
              break;
            }
            case 4: {
              // All edges are split => regular refinement
              finest_cell_ci[cell_index].ref_pat_ = RefPat::rp_regular;
              break;
            }
            default: {
              LF_VERIFY_MSG(false, "Illegal number " << split_edge_cnt
                                                     << " of split edges");
              break;
            }
          }  // end switch split_edge_cnt
          break;
        }  // end case of a quadrilateral
        default: {
          LF_VERIFY_MSG(false, "Illegal cell type");
          break;
        }
      }  // end switch cell type
    }    // end loop over cells
  } while (!refinement_complete);

  PerformRefinement();
}  // end RefineMarked

// NOLINTNEXTLINE(google-readability-function-size, hicpp-function-size, readability-function-size)
void MeshHierarchy::PerformRefinement() {
  CONTROLLEDSTATEMENT(output_ctrl_, 10,
                      std::cout << "Entering MeshHierarchy::PerformRefinement: "
                                << meshes_.size() << " levels" << std::endl;)
  // This function expects that the refinement patterns  stored in the
  // vectors point_child_infos_, edge_child_infos_ and cell_child_infos_
  // have been initialized consistently for the finest mesh.

  // This functions will augement these vectors with information about
  // the indices of child entities on a newly created finest mesh.

  // Retrieve the finest mesh in the hierarchy = parent mesh
  const mesh::Mesh &parent_mesh(*meshes_.back());
  // First run through the vertices, create child vertices and register
  // them with the mesh factory
  // Store child indices in an auxiliary array
  {
    std::vector<PointChildInfo> &pt_child_info(point_child_infos_.back());
    size_type new_node_cnt = 0;
    const Hybrid2DRefinementPattern rp_copy_node(lf::base::RefEl::kPoint(),
                                                 RefPat::rp_copy);
    for (const mesh::Entity &node : parent_mesh.Entities(2)) {
      // Obtain index of node in coarse mesh
      const lf::base::glb_idx_t node_index = parent_mesh.Index(node);
      LF_VERIFY_MSG(node_index < pt_child_info.size(),
                    "Node index " << node_index << " out of range");
      // Find position of node in physical coordinates
      const lf::geometry::Geometry &pt_geo(*node.Geometry());
      if (pt_child_info[node_index].ref_pat != RefPat::rp_nil) {
        // Generate a node for the fine mesh at the same position
        std::vector<std::unique_ptr<geometry::Geometry>> pt_child_geo_ptrs(
            pt_geo.ChildGeometry(rp_copy_node, 0));
        LF_VERIFY_MSG(pt_child_geo_ptrs.size() == 1,
                      "A point can only have one chile");
        pt_child_info[node_index].child_point_idx =
            mesh_factory_->AddPoint(std::move(pt_child_geo_ptrs[0]));
        new_node_cnt++;
      }
    }  // end loop over nodes
    CONTROLLEDSTATEMENT(
        output_ctrl_, 10,
        std::cout << new_node_cnt << " new nodes added" << std::endl;)

    // Now traverse the edges. Depending on the refinement pattern,
    // either copy them or split them.
    // Supplement the refinement information for edges accordingly.
    std::vector<EdgeChildInfo> &ed_child_info(edge_child_infos_.back());
    size_type new_edge_cnt = 0;
    for (const mesh::Entity &edge : parent_mesh.Entities(1)) {
      // Fetch global index of edge
      lf::base::glb_idx_t edge_index = parent_mesh.Index(edge);

      // Get indices of endpoints in parent mesh
      auto ed_nodes(edge.SubEntities(1));
      const lf::base::glb_idx_t ed_p0_coarse_idx =
          parent_mesh.Index(ed_nodes[0]);
      const lf::base::glb_idx_t ed_p1_coarse_idx =
          parent_mesh.Index(ed_nodes[1]);

      // Obtain indices of the nodes at the same position in the fine mesh
      const lf::base::glb_idx_t ed_p0_fine_idx =
          pt_child_info[ed_p0_coarse_idx].child_point_idx;
      const lf::base::glb_idx_t ed_p1_fine_idx =
          pt_child_info[ed_p1_coarse_idx].child_point_idx;
      // Prepare request of geometry after refinement
      EdgeChildInfo &edge_ci(ed_child_info[edge_index]);
      const RefPat edge_refpat(edge_ci.ref_pat_);
      Hybrid2DRefinementPattern rp(edge.RefEl(), edge_refpat);

      // Distinguish between different local refinement patterns
      switch (edge_refpat) {
        case RefPat::rp_copy: {
          // Edge has to be duplicated
          std::vector<std::unique_ptr<lf::geometry::Geometry>> ed_copy(
              edge.Geometry()->ChildGeometry(rp, 0));
          LF_VERIFY_MSG(ed_copy.size() == 1,
                        "Copy may create only a single child!");
          // Register the new edge
          CONTROLLEDSTATEMENT(output_ctrl_, 50,
                              std::cout << "Copy edge " << edge_index
                                        << " new edge [" << ed_p0_fine_idx
                                        << "," << ed_p1_fine_idx << "] "
                                        << std::endl;)

          edge_ci.child_edge_idx.push_back(mesh_factory_->AddEntity(
              edge.RefEl(),
              lf::base::ForwardRange<const lf::base::glb_idx_t>(
                  {ed_p0_fine_idx, ed_p1_fine_idx}),
              std::move(ed_copy[0])));
          break;
        }  // end rp_copy
        case rp_split: {
          // Edge has to be split into two halves with a new node created at the
          // midpoint position. First obtain geometry information about new node
          // (sub-entity with relative co-dimension 1)
          std::vector<std::unique_ptr<geometry::Geometry>> edge_nodes_geo_ptrs(
              edge.Geometry()->ChildGeometry(rp, 1));
          LF_VERIFY_MSG(edge_nodes_geo_ptrs.size() == 1,
                        "Split edge with " << edge_nodes_geo_ptrs.size()
                                           << " child nodes!");
          // Register midpoint as new node
          const lf::base::glb_idx_t midpoint_fine_idx =
              mesh_factory_->AddPoint(std::move(edge_nodes_geo_ptrs[0]));
          edge_ci.child_point_idx.push_back(midpoint_fine_idx);
          // Next get the geometry objects for the two child edges (co-dim == 0)
          std::vector<std::unique_ptr<geometry::Geometry>> edge_child_geo_ptrs(
              edge.Geometry()->ChildGeometry(rp, 0));
          LF_VERIFY_MSG(
              edge_child_geo_ptrs.size() == 2,
              "Split edge with " << edge_child_geo_ptrs.size() << " parts!");
          // Register the two new edges
          // CAREFUL: Assignment of endpoints has to match implementation in
          // refinement.cc
          CONTROLLEDSTATEMENT(output_ctrl_, 50,
                              std::cout << "Split Edge " << edge_index
                                        << " new edges [" << ed_p0_fine_idx
                                        << "," << midpoint_fine_idx << "], ["
                                        << midpoint_fine_idx << ","
                                        << ed_p1_fine_idx << "] " << std::endl;)

          edge_ci.child_edge_idx.push_back(mesh_factory_->AddEntity(
              edge.RefEl(),
              lf::base::ForwardRange<const lf::base::glb_idx_t>(
                  {ed_p0_fine_idx, midpoint_fine_idx}),
              std::move(edge_child_geo_ptrs[0])));
          edge_ci.child_edge_idx.push_back(mesh_factory_->AddEntity(
              edge.RefEl(),
              lf::base::ForwardRange<const lf::base::glb_idx_t>(
                  {midpoint_fine_idx, ed_p1_fine_idx}),
              std::move(edge_child_geo_ptrs[1])));
          break;
        }  // end rp_split
        default: {
          LF_VERIFY_MSG(false, "Refinement pattern "
                                   << static_cast<int>(edge_refpat)
                                   << " illegal for edge");
          break;
        }
      }  // end switch refpat
      new_edge_cnt++;
    }  // end edge loop

    CONTROLLEDSTATEMENT(
        output_ctrl_, 50,
        std::cout << new_edge_cnt << " edges added " << std::endl;)

    // Visit all cells, examine their refinement patterns, retrieve indices of
    // their sub-entities, and those of the children.
    std::vector<CellChildInfo> &cell_child_info(cell_child_infos_.back());
    LF_VERIFY_MSG(cell_child_info.size() == parent_mesh.NumEntities(0),
                  "Size mismatch for CellChildInfos");
    for (const mesh::Entity &cell : parent_mesh.Entities(0)) {
      // type of cell
      const lf::base::RefEl ref_el(cell.RefEl());
      const lf::base::size_type num_edges = ref_el.NumSubEntities(1);
      const lf::base::size_type num_vertices = ref_el.NumSubEntities(2);
      // fetch index of current cell
      const lf::base::glb_idx_t cell_index(parent_mesh.Index(cell));

      // Set up refinement object -> variable rp
      CellChildInfo &cell_ci(cell_child_info[cell_index]);
      const RefPat cell_refpat(cell_ci.ref_pat_);
      const sub_idx_t anchor = cell_ci.anchor_;  // anchor edge index

      CONTROLLEDSTATEMENT(
          output_ctrl_, 50,
          std::cout << "Cell " << cell_index << " = " << ref_el.ToString()
                    << ", refpat = " << static_cast<int>(cell_refpat)
                    << ", anchor = " << anchor << std::endl;)

      Hybrid2DRefinementPattern rp(cell.RefEl(), cell_refpat, anchor);

      // Index offsets for refinement patterns requiring an ancchor edge
      std::array<sub_idx_t, 4> mod{};
      if (anchor != idx_nil) {
        for (int k = 0; k < num_vertices; k++) {
          mod[k] = (k + anchor) % num_vertices;
        }
      }

      // Obtain indices of subentities (co-dimension = outer array index)
      std::array<std::vector<lf::base::glb_idx_t>, 3> cell_subent_idx;
      cell_subent_idx[0].push_back(cell_index);
      for (int codim = 1; codim <= 2; codim++) {
        base::RandomAccessRange<const mesh::Entity> subentities(
            cell.SubEntities(codim));
        for (const mesh::Entity &sub_ent : subentities) {
          cell_subent_idx[codim].push_back(parent_mesh.Index(sub_ent));
        }
        LF_VERIFY_MSG(
            cell_subent_idx[codim].size() == ref_el.NumSubEntities(codim),
            ref_el.ToString() << ": only " << cell_subent_idx[codim].size()
                              << " subents of codim = " << codim);
        CONTROLLEDSTATEMENT(
            output_ctrl_, 50,
            std::cout << " Subent(" << codim << ") = [" << std::flush;
            for (unsigned int j
                 : cell_subent_idx[codim]) {
              std::cout << j << "," << std::flush;
            } std::cout
            << "], " << std::flush;
            std::cout << std::endl;)
      }
      // Index information for sub-entities with respect  to fine mesh
      // Retrieve indices of vertices of cell on the fine mesh
      std::array<lf::base::glb_idx_t, 4> vertex_child_idx{
          {idx_nil, idx_nil, idx_nil, idx_nil}};
      for (lf::base::sub_idx_t vt_lidx = 0; vt_lidx < num_vertices; vt_lidx++) {
        LF_VERIFY_MSG(pt_child_info[cell_subent_idx[2][vt_lidx]].ref_pat ==
                          RefPat::rp_copy,
                      "Vertex must have been copied!");
        vertex_child_idx[vt_lidx] =
            pt_child_info[cell_subent_idx[2][vt_lidx]].child_point_idx;
      }

      CONTROLLEDSTATEMENT(
          output_ctrl_, 50, std::cout << ", vt_child_idx = [" << std::flush;
          for (int j = 0; j < num_vertices; j++) {
            std::cout << vertex_child_idx[j] << "," << std::flush;
          } std::cout
          << "], " << std::flush;)

      // Retrieve indices of midpoints of edges, if they exist
      std::array<lf::base::glb_idx_t, 4> edge_midpoint_idx{
          {idx_nil, idx_nil, idx_nil, idx_nil}};
      for (lf::base::sub_idx_t ed_lidx = 0; ed_lidx < num_edges; ed_lidx++) {
        const EdgeChildInfo &ed_ci(ed_child_info[cell_subent_idx[1][ed_lidx]]);
        if (!ed_ci.child_point_idx.empty()) {
          LF_VERIFY_MSG(ed_child_info[cell_subent_idx[1][ed_lidx]].ref_pat_ ==
                            RefPat::rp_split,
                        "Edge with a midpoint must have been split");
          edge_midpoint_idx[ed_lidx] = ed_ci.child_point_idx[0];
        }
      }  // end loop over local edges

      CONTROLLEDSTATEMENT(
          output_ctrl_, 50, std::cout << ", ed_mp_idx = [" << std::flush;
          for (int j = 0; j < num_edges; j++) {
            std::cout << edge_midpoint_idx[j] << "," << std::flush;
          } std::cout
          << "], " << std::endl;)

      // Array of node indices (w.r.t. fine mesh) for sub-cells (triangles or
      // quadrilaterals)
      std::vector<std::vector<glb_idx_t>> child_cell_nodes;
      std::vector<glb_idx_t> tria_ccn_tmp(3);
      std::vector<glb_idx_t> quad_ccn_tmp(4);
      // Array of node indices for interior child edges
      std::vector<std::array<glb_idx_t, 2>> child_edge_nodes;
      std::array<glb_idx_t, 2> cen_tmp{};

      // Local linking of child entities is very different depending on cell
      // type and refinement patterns.
      if (ref_el == lf::base::RefEl::kTria()) {
        // Case of a triangle
        switch (cell_refpat) {
          case RefPat::rp_nil: {
            LF_VERIFY_MSG(false, "Every triangle must be refined");
            break;
          }
          case RefPat::rp_copy: {
            // Just copy the parent triangle
            child_cell_nodes.push_back({vertex_child_idx[0],
                                        vertex_child_idx[1],
                                        vertex_child_idx[2]});
            break;
          }
          case RefPat::rp_bisect: {
            LF_VERIFY_MSG(
                anchor != idx_nil,
                "Anchor must be set for bisection refinement of triangle");
            // Splitting a triangle in two by bisecting anchor edge
            tria_ccn_tmp[0] = vertex_child_idx[mod[0]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = vertex_child_idx[mod[2]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = vertex_child_idx[mod[2]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            // edges
            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = vertex_child_idx[mod[2]];
            child_edge_nodes.push_back(cen_tmp);
            break;
          }
          case RefPat::rp_trisect: {
            LF_VERIFY_MSG(
                anchor != idx_nil,
                "Anchor must be set for trisection refinement of triangle");
            // Bisect through anchor edge first and then bisect through
            // edge with the next larger index (mod 3); creates three
            // child triangles.
            tria_ccn_tmp[0] = vertex_child_idx[mod[0]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = vertex_child_idx[mod[2]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[2]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            // edges
            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = vertex_child_idx[mod[2]];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = edge_midpoint_idx[mod[1]];
            child_edge_nodes.push_back(cen_tmp);
            break;
          }
          case RefPat::rp_trisect_left: {
            LF_VERIFY_MSG(
                anchor != idx_nil,
                "Anchor must be set for trisection refinement of triangle");

            // Bisect through anchor edge first and then bisect through
            // edge with the next smaller index (mod 3); creates three
            // child triangles.
            tria_ccn_tmp[0] = vertex_child_idx[mod[0]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[2]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = vertex_child_idx[mod[2]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[2]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[2]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            // edges
            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = vertex_child_idx[mod[2]];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = edge_midpoint_idx[mod[2]];
            child_edge_nodes.push_back(cen_tmp);
            break;
          }
          case RefPat::rp_quadsect: {
            LF_VERIFY_MSG(
                anchor != idx_nil,
                "Anchor must be set for quadsection refinement of triangle");
            // Bisect through the anchor edge first and then
            // through the two remaining edges; creates four child
            // triangles; every edge is split.
            tria_ccn_tmp[0] = vertex_child_idx[mod[0]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[2]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[2]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[2]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[2]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            // edges
            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = vertex_child_idx[mod[2]];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = edge_midpoint_idx[mod[2]];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = edge_midpoint_idx[mod[1]];
            child_edge_nodes.push_back(cen_tmp);
            break;
          }
          case RefPat::rp_barycentric: {
            //  Split triangle into 6 smaller triangles by connecting
            // the center of gravity with the vertices and the midpoints
            // of the edges.
            // Create a new interior vertex
            std::vector<std::unique_ptr<geometry::Geometry>>
                cell_center_geo_ptrs(cell.Geometry()->ChildGeometry(
                    rp, 2));  // point: co-dim == 2
            LF_VERIFY_MSG(cell_center_geo_ptrs.size() == 1,
                          "Barycentrically refined triangle with "
                              << cell_center_geo_ptrs.size()
                              << " interior child nodes ??");
            // Register midpoint as new node
            const glb_idx_t center_fine_idx =
                mesh_factory_->AddPoint(std::move(cell_center_geo_ptrs[0]));
            cell_ci.child_point_idx.push_back(center_fine_idx);

            tria_ccn_tmp[0] = vertex_child_idx[0];
            tria_ccn_tmp[1] = edge_midpoint_idx[0];
            tria_ccn_tmp[2] = center_fine_idx;
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[1];
            tria_ccn_tmp[1] = edge_midpoint_idx[0];
            tria_ccn_tmp[2] = center_fine_idx;
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[1];
            tria_ccn_tmp[1] = edge_midpoint_idx[1];
            tria_ccn_tmp[2] = center_fine_idx;
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[2];
            tria_ccn_tmp[1] = edge_midpoint_idx[1];
            tria_ccn_tmp[2] = center_fine_idx;
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[2];
            tria_ccn_tmp[1] = edge_midpoint_idx[2];
            tria_ccn_tmp[2] = center_fine_idx;
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[0];
            tria_ccn_tmp[1] = edge_midpoint_idx[2];
            tria_ccn_tmp[2] = center_fine_idx;
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            // edges
            cen_tmp[0] = vertex_child_idx[0];
            cen_tmp[1] = center_fine_idx;
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = vertex_child_idx[1];
            cen_tmp[1] = center_fine_idx;
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = vertex_child_idx[2];
            cen_tmp[1] = center_fine_idx;
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[0];
            cen_tmp[1] = center_fine_idx;
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[1];
            cen_tmp[1] = center_fine_idx;
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[2];
            cen_tmp[1] = center_fine_idx;
            child_edge_nodes.push_back(cen_tmp);
            break;
          }
          case RefPat::rp_regular: {
            // regular refinement into four congruent triangles
            // Child cells
            tria_ccn_tmp[0] = vertex_child_idx[0];
            tria_ccn_tmp[1] = edge_midpoint_idx[0];
            tria_ccn_tmp[2] = edge_midpoint_idx[2];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[1];
            tria_ccn_tmp[1] = edge_midpoint_idx[0];
            tria_ccn_tmp[2] = edge_midpoint_idx[1];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[2];
            tria_ccn_tmp[1] = edge_midpoint_idx[2];
            tria_ccn_tmp[2] = edge_midpoint_idx[1];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = edge_midpoint_idx[0];
            tria_ccn_tmp[1] = edge_midpoint_idx[1];
            tria_ccn_tmp[2] = edge_midpoint_idx[2];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            // Child edges (interior edges)
            cen_tmp[0] = edge_midpoint_idx[0];
            cen_tmp[1] = edge_midpoint_idx[2];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[0];
            cen_tmp[1] = edge_midpoint_idx[1];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[2];
            cen_tmp[1] = edge_midpoint_idx[1];
            child_edge_nodes.push_back(cen_tmp);
            // No new interior node
            break;
          }
          default: {
            LF_VERIFY_MSG(
                false, "Refinement pattern not (yet) implemented for triangle");
            break;
          }
        }  // end switch refpat
      } else if (ref_el == lf::base::RefEl::kQuad()) {
        // Case of a quadrilateral
        switch (cell_refpat) {
          case RefPat::rp_nil: {
            LF_VERIFY_MSG(false, "Every quadrilateral must be refined");
            break;
          }
          case RefPat::rp_copy: {
            // Just copy the parent triangle
            child_cell_nodes.push_back(
                {vertex_child_idx[0], vertex_child_idx[1], vertex_child_idx[2],
                 vertex_child_idx[3]});
            break;
          }
          case RefPat::rp_trisect: {
            LF_VERIFY_MSG(
                anchor != idx_nil,
                "Anchor must be set for trisection refinement of quad");

            // Partition a quad into three triangle, the anchor edge
            // being split in the process
            tria_ccn_tmp[0] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[1] = vertex_child_idx[mod[2]];
            tria_ccn_tmp[2] = vertex_child_idx[mod[3]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[1] = vertex_child_idx[mod[0]];
            tria_ccn_tmp[2] = vertex_child_idx[mod[3]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[1] = vertex_child_idx[mod[1]];
            tria_ccn_tmp[2] = vertex_child_idx[mod[2]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            // edges
            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = vertex_child_idx[mod[2]];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = vertex_child_idx[mod[3]];
            child_edge_nodes.push_back(cen_tmp);
            break;
          }
          case RefPat::rp_quadsect: {
            LF_VERIFY_MSG(
                anchor != idx_nil,
                "Anchor must be set for quadsection refinement of triangle");
            // Partition a quad into four triangle, thus
            // splitting two edges. The one with the smaller sub index is the
            // anchor edge
            tria_ccn_tmp[0] = vertex_child_idx[mod[0]];
            tria_ccn_tmp[1] = vertex_child_idx[mod[3]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[0]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[1]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[0]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[2]];
            tria_ccn_tmp[1] = vertex_child_idx[mod[3]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[1]];
            tria_ccn_tmp[2] = vertex_child_idx[mod[3]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            // edges
            cen_tmp[0] = vertex_child_idx[mod[3]];
            cen_tmp[1] = edge_midpoint_idx[mod[0]];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[mod[1]];
            cen_tmp[1] = edge_midpoint_idx[mod[0]];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = vertex_child_idx[mod[3]];
            cen_tmp[1] = edge_midpoint_idx[mod[1]];
            child_edge_nodes.push_back(cen_tmp);
            break;
          }
          case RefPat::rp_bisect:
          case RefPat::rp_split: {
            LF_VERIFY_MSG(anchor != idx_nil,
                          "Anchor must be set for splitting of quad");

            // Cut a quadrilateral into two
            quad_ccn_tmp[0] = vertex_child_idx[mod[0]];
            quad_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            quad_ccn_tmp[2] = edge_midpoint_idx[mod[2]];
            quad_ccn_tmp[3] = vertex_child_idx[mod[3]];
            LF_ASSERT_MSG(checkValidIndex(quad_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(quad_ccn_tmp);
            ;

            quad_ccn_tmp[0] = vertex_child_idx[mod[1]];
            quad_ccn_tmp[1] = vertex_child_idx[mod[2]];
            quad_ccn_tmp[2] = edge_midpoint_idx[mod[2]];
            quad_ccn_tmp[3] = edge_midpoint_idx[mod[0]];
            LF_ASSERT_MSG(checkValidIndex(quad_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(quad_ccn_tmp);
            ;

            // edges
            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = edge_midpoint_idx[mod[2]];
            child_edge_nodes.push_back(cen_tmp);
            break;
          }
          case RefPat::rp_threeedge: {
            LF_VERIFY_MSG(
                anchor != idx_nil,
                "Anchor must be set for three edge refinement of a quad");

            quad_ccn_tmp[0] = vertex_child_idx[mod[2]];
            quad_ccn_tmp[1] = vertex_child_idx[mod[3]];
            quad_ccn_tmp[2] = edge_midpoint_idx[mod[3]];
            quad_ccn_tmp[3] = edge_midpoint_idx[mod[1]];
            LF_ASSERT_MSG(checkValidIndex(quad_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(quad_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[0]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[3]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            tria_ccn_tmp[0] = edge_midpoint_idx[mod[0]];
            tria_ccn_tmp[1] = edge_midpoint_idx[mod[1]];
            tria_ccn_tmp[2] = edge_midpoint_idx[mod[3]];
            LF_ASSERT_MSG(checkValidIndex(tria_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(tria_ccn_tmp);
            ;

            // edges
            cen_tmp[0] = edge_midpoint_idx[mod[3]];
            cen_tmp[1] = edge_midpoint_idx[mod[1]];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = edge_midpoint_idx[mod[3]];
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[mod[0]];
            cen_tmp[1] = edge_midpoint_idx[mod[1]];
            child_edge_nodes.push_back(cen_tmp);
            break;
          }
          case RefPat::rp_barycentric:
          case RefPat::rp_regular: {
            // Create a new interior vertex
            std::vector<std::unique_ptr<geometry::Geometry>>
                cell_center_geo_ptrs(cell.Geometry()->ChildGeometry(
                    rp, 2));  // point: co-dim == 2
            LF_VERIFY_MSG(cell_center_geo_ptrs.size() == 1,
                          "Regularly refined quadrilateral with "
                              << cell_center_geo_ptrs.size()
                              << " interior child nodes!");
            // Register midpoint as new node
            const glb_idx_t center_fine_idx =
                mesh_factory_->AddPoint(std::move(cell_center_geo_ptrs[0]));
            cell_ci.child_point_idx.push_back(center_fine_idx);

            // Set the node indices (w.r.t. fine mesh) of the four sub-quads
            quad_ccn_tmp[0] = vertex_child_idx[0];
            quad_ccn_tmp[1] = edge_midpoint_idx[0];
            quad_ccn_tmp[2] = center_fine_idx;
            quad_ccn_tmp[3] = edge_midpoint_idx[3];
            LF_ASSERT_MSG(checkValidIndex(quad_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(quad_ccn_tmp);
            ;

            quad_ccn_tmp[0] = vertex_child_idx[1];
            quad_ccn_tmp[1] = edge_midpoint_idx[1];
            quad_ccn_tmp[2] = center_fine_idx;
            quad_ccn_tmp[3] = edge_midpoint_idx[0];
            LF_ASSERT_MSG(checkValidIndex(quad_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(quad_ccn_tmp);
            ;

            quad_ccn_tmp[0] = vertex_child_idx[2];
            quad_ccn_tmp[1] = edge_midpoint_idx[1];
            quad_ccn_tmp[2] = center_fine_idx;
            quad_ccn_tmp[3] = edge_midpoint_idx[2];
            LF_ASSERT_MSG(checkValidIndex(quad_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(quad_ccn_tmp);
            ;

            quad_ccn_tmp[0] = vertex_child_idx[3];
            quad_ccn_tmp[1] = edge_midpoint_idx[2];
            quad_ccn_tmp[2] = center_fine_idx;
            quad_ccn_tmp[3] = edge_midpoint_idx[3];
            LF_ASSERT_MSG(checkValidIndex(quad_ccn_tmp),
                          "refpat = " << (int)cell_refpat << ": illegal index");
            child_cell_nodes.push_back(quad_ccn_tmp);
            ;

            // set node indices of the four new interior edges
            cen_tmp[0] = edge_midpoint_idx[0];
            cen_tmp[1] = center_fine_idx;
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[1];
            cen_tmp[1] = center_fine_idx;
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[2];
            cen_tmp[1] = center_fine_idx;
            child_edge_nodes.push_back(cen_tmp);

            cen_tmp[0] = edge_midpoint_idx[3];
            cen_tmp[1] = center_fine_idx;
            child_edge_nodes.push_back(cen_tmp);
            break;
          }
          default: {
            LF_VERIFY_MSG(
                false,
                "Refinement pattern not (yet) implemented for quadrilateral");
            break;
          }
        }  // end switch cell_refpat
      }    // end quadrilateral case
      else {
        LF_VERIFY_MSG(false, "Unknown cell type" << ref_el.ToString());
      }
      // Register new edges
      {
        std::vector<std::unique_ptr<geometry::Geometry>> cell_edge_geo_ptrs(
            cell.Geometry()->ChildGeometry(rp, 1));  // child edge: co-dim == 2
        const size_type num_new_edges = child_edge_nodes.size();
        LF_VERIFY_MSG(
            num_new_edges == rp.noChildren(1),
            "num_new_edges = " << num_new_edges << " <-> " << rp.noChildren(1));
        for (int k = 0; k < num_new_edges; k++) {
          const std::array<glb_idx_t, 2> &cen(child_edge_nodes[k]);
          CONTROLLEDSTATEMENT(output_ctrl_, 50,
                              std::cout
                                  << ref_el.ToString() << "(" << cell_index
                                  << "), ref_pat = " << (int)cell_refpat
                                  << ": new edge " << k << "[" << cen[0] << ","
                                  << cen[1] << "]" << std::endl;)

          const glb_idx_t new_edge_index = mesh_factory_->AddEntity(
              lf::base::RefEl::kSegment(), {cen[0], cen[1]},
              std::move(cell_edge_geo_ptrs[k]));
          cell_ci.child_edge_idx.push_back(new_edge_index);
        }  // end loop over new edges
      }    // end register new edges
      // Register new cells
      {
        std::vector<std::unique_ptr<geometry::Geometry>> childcell_geo_ptrs(
            cell.Geometry()->ChildGeometry(rp, 0));  // child cell: co-dim == 0
        const size_type num_new_cells = child_cell_nodes.size();
        LF_VERIFY_MSG(
            num_new_cells == rp.noChildren(0),
            "num_new_cells = " << num_new_cells << " <-> " << rp.noChildren(0));
        for (int k = 0; k < num_new_cells; k++) {
          const std::vector<glb_idx_t> &ccn(child_cell_nodes[k]);
          glb_idx_t new_cell_index;
          if (ccn.size() == 3) {
            // New cell is a triangle
            CONTROLLEDSTATEMENT(
                output_ctrl_, 50,
                std::cout << ref_el.ToString() << "(" << cell_index
                          << "), ref_pat = " << (int)cell_refpat
                          << ": new triangle " << k << " [" << ccn[0] << ","
                          << ccn[1] << "," << ccn[2] << "]" << std::endl;)

            new_cell_index = mesh_factory_->AddEntity(
                lf::base::RefEl::kTria(), {ccn[0], ccn[1], ccn[2]},
                std::move(childcell_geo_ptrs[k]));
          } else if (ccn.size() == 4) {
            // New cell is a quadrilateral
            CONTROLLEDSTATEMENT(output_ctrl_, 50,
                                std::cout
                                    << ref_el.ToString() << "(" << cell_index
                                    << "), ref_pat = " << (int)cell_refpat
                                    << ": new quad " << k << " [" << ccn[0]
                                    << "," << ccn[1] << "," << ccn[2] << ","
                                    << ccn[3] << "]" << std::endl;)

            new_cell_index = mesh_factory_->AddEntity(
                lf::base::RefEl::kQuad(), {ccn[0], ccn[1], ccn[2], ccn[3]},
                std::move(childcell_geo_ptrs[k]));
          } else {
            LF_VERIFY_MSG(false,
                          "Child cells must be either triangles or quads");
          }
          cell_ci.child_cell_idx.push_back(new_cell_index);
        }  // end loop over new cells
      }    // end register new cells
    }      // end loop over cells
  }
  // At this point the MeshFactory has complete information to generate the new
  // finest mesh
  meshes_.push_back(mesh_factory_->Build());  // MESH CONSTRUCTION
  mesh::Mesh &child_mesh(*meshes_.back());

  CONTROLLEDSTATEMENT(output_ctrl_, 10,
                      std::cout << "Child mesh" << child_mesh.NumEntities(2)
                                << " nodes, " << child_mesh.NumEntities(1)
                                << " edges, " << child_mesh.NumEntities(0)
                                << " cells." << std::endl;)

  // Create space for data pertaining to the new mesh
  // Note that references to vectors may become invalid
  {
    // Arrays containing parent information
    parent_infos_.push_back(
        {std::move(std::vector<ParentInfo>(child_mesh.NumEntities(0))),
         std::move(std::vector<ParentInfo>(child_mesh.NumEntities(1))),
         std::move(std::vector<ParentInfo>(child_mesh.NumEntities(2)))});

    // Initialize (empty) refinement information for newly created fine mesh
    // Same code as in the constructor of MeshHierarchy
    std::vector<CellChildInfo> fine_cell_child_info(child_mesh.NumEntities(0));
    std::vector<EdgeChildInfo> fine_edge_child_info(child_mesh.NumEntities(1));
    std::vector<PointChildInfo> fine_point_child_info(
        child_mesh.NumEntities(2));
    cell_child_infos_.push_back(std::move(fine_cell_child_info));
    edge_child_infos_.push_back(std::move(fine_edge_child_info));
    point_child_infos_.push_back(std::move(fine_point_child_info));

    // Array containing information about refinement edges
    refinement_edges_.push_back(
        std::move(std::vector<sub_idx_t>(child_mesh.NumEntities(0), idx_nil)));

    // Finally set up vector for edge flags
    edge_marked_.emplace_back(child_mesh.NumEntities(1), false);
  }

  // Finally, we have to initialize the parent pointers for the entities of the
  // newly created finest mesh.
  {
    const size_type n_levels = meshes_.size();
    LF_VERIFY_MSG(n_levels > 1, "A least two levels after refinement");
    // Access parent information for finest level
    std::vector<ParentInfo> &fine_node_parent_info(parent_infos_.back()[2]);
    std::vector<ParentInfo> &fine_edge_parent_info(parent_infos_.back()[1]);
    std::vector<ParentInfo> &fine_cell_parent_info(parent_infos_.back()[0]);

    LF_VERIFY_MSG(fine_node_parent_info.size() == child_mesh.NumEntities(2),
                  "fine_node_parent_info size mismatch");
    LF_VERIFY_MSG(fine_edge_parent_info.size() == child_mesh.NumEntities(1),
                  "fine_edge_parent_info size mismatch");
    LF_VERIFY_MSG(fine_cell_parent_info.size() == child_mesh.NumEntities(0),
                  "fine_edge_parent_info size mismatch");

    // Visit all nodes of the parent mesh and retrieve their children by index
    std::vector<PointChildInfo> &pt_child_info(
        point_child_infos_.at(n_levels - 2));
    for (const mesh::Entity &node : parent_mesh.Entities(2)) {
      // Obtain index of node in coarse mesh
      const glb_idx_t node_index = parent_mesh.Index(node);
      const glb_idx_t child_index = pt_child_info[node_index].child_point_idx;
      if (child_index != idx_nil) {
        LF_VERIFY_MSG(child_index < fine_node_parent_info.size(),
                      "index " << child_index << " illegal for child vertex");
        fine_node_parent_info[child_index].child_number = 0;
        fine_node_parent_info[child_index].parent_ptr = &node;
        fine_node_parent_info[child_index].parent_index = node_index;
      }
    }  // end loop over nodes

    // Traverse edges and set the parent for their interior nodes and child
    // edges
    std::vector<EdgeChildInfo> &ed_child_info(
        edge_child_infos_.at(n_levels - 2));
    for (const mesh::Entity &edge : parent_mesh.Entities(1)) {
      // Obtain index of edge in coarse mesh
      const glb_idx_t edge_index = parent_mesh.Index(edge);
      const EdgeChildInfo &ci_edge(ed_child_info[edge_index]);

      // Visit child edges
      const size_type num_child_edges = ci_edge.child_edge_idx.size();
      for (int l = 0; l < num_child_edges; l++) {
        const glb_idx_t edge_child_idx = ci_edge.child_edge_idx[l];
        LF_VERIFY_MSG(edge_child_idx < fine_edge_parent_info.size(),
                      "index " << edge_child_idx << " illegal for child edge");
        fine_edge_parent_info[edge_child_idx].child_number = l;
        fine_edge_parent_info[edge_child_idx].parent_ptr = &edge;
        fine_edge_parent_info[edge_child_idx].parent_index = edge_index;
      }  // end loop over child edges

      // Visit child nodes
      const size_type num_child_nodes = ci_edge.child_point_idx.size();
      for (int l = 0; l < num_child_nodes; l++) {
        const glb_idx_t node_child_idx = ci_edge.child_point_idx[l];
        LF_VERIFY_MSG(node_child_idx < fine_node_parent_info.size(),
                      "index " << node_child_idx << " illegal for child point");
        fine_node_parent_info[node_child_idx].child_number = l;
        fine_node_parent_info[node_child_idx].parent_ptr = &edge;
        fine_node_parent_info[node_child_idx].parent_index = edge_index;
      }  // end loop over child points
    }    // end loop over edges

    // Loop over cells
    std::vector<CellChildInfo> &cell_child_info(
        cell_child_infos_.at(n_levels - 2));
    for (const mesh::Entity &cell : parent_mesh.Entities(0)) {
      // fetch index of current cell
      const glb_idx_t cell_index(parent_mesh.Index(cell));
      // Get infformation about children
      CellChildInfo &cell_ci(cell_child_info[cell_index]);

      // Visit child cells
      const size_type num_child_cells = cell_ci.child_cell_idx.size();
      for (int l = 0; l < num_child_cells; l++) {
        const glb_idx_t child_cell_idx = cell_ci.child_cell_idx[l];
        LF_VERIFY_MSG(child_cell_idx < fine_cell_parent_info.size(),
                      "index " << child_cell_idx << " illegal for child cell");
        fine_cell_parent_info[child_cell_idx].child_number = l;
        fine_cell_parent_info[child_cell_idx].parent_ptr = &cell;
        fine_cell_parent_info[child_cell_idx].parent_index = cell_index;
      }

      // Visit child edges
      const size_type num_child_edges = cell_ci.child_edge_idx.size();
      for (int l = 0; l < num_child_edges; l++) {
        const glb_idx_t edge_child_idx = cell_ci.child_edge_idx[l];
        LF_VERIFY_MSG(edge_child_idx < fine_edge_parent_info.size(),
                      "index " << edge_child_idx << " illegal for child edge");
        fine_edge_parent_info[edge_child_idx].child_number = l;
        fine_edge_parent_info[edge_child_idx].parent_ptr = &cell;
        fine_edge_parent_info[edge_child_idx].parent_index = cell_index;
      }  // end loop over child edges

      // Visit child nodes
      const size_type num_child_nodes = cell_ci.child_point_idx.size();
      for (int l = 0; l < num_child_nodes; l++) {
        const glb_idx_t node_child_idx = cell_ci.child_point_idx[l];
        LF_VERIFY_MSG(node_child_idx < fine_node_parent_info.size(),
                      "index " << node_child_idx << " illegal for child point");
        fine_node_parent_info[node_child_idx].child_number = l;
        fine_node_parent_info[node_child_idx].parent_ptr = &cell;
        fine_node_parent_info[node_child_idx].parent_index = cell_index;
      }  // end loop over child points
    }    // end loop over cells of parent mesh
  }

  // Finally set refinement edges for fine mesh
  {
    size_type n_levels = meshes_.size();
    LF_VERIFY_MSG(n_levels > 1, "At least two levels after refinement!");
    std::vector<sub_idx_t> &child_ref_edges(refinement_edges_.at(n_levels - 1));
    std::vector<ParentInfo> &cell_parent_info(
        parent_infos_.at(n_levels - 1)[0]);
    std::vector<CellChildInfo> &parent_cell_ci(
        cell_child_infos_.at(n_levels - 2));
    std::vector<sub_idx_t> &parent_ref_edges(
        refinement_edges_.at(n_levels - 2));

    // Traverse the cells of the fine mesh
    CONTROLLEDSTATEMENT(output_ctrl_, 10,
                        std::cout << "Setting refinement edges" << std::endl;)

    for (const mesh::Entity &fine_cell : child_mesh.Entities(0)) {
      const glb_idx_t cell_index = child_mesh.Index(fine_cell);
      child_ref_edges[cell_index] = idx_nil;
      // Refinement edge relevant for triangles onyl
      if (fine_cell.RefEl() == lf::base::RefEl::kTria()) {
        // pointer to cell whose refinement has created the current one
        const mesh::Entity *parent_ptr =
            cell_parent_info[cell_index].parent_ptr;
        LF_VERIFY_MSG(parent_ptr != nullptr,
                      "Cell " << cell_index << " has no parent, paren_index = "
                              << cell_parent_info[cell_index].parent_index);

        if (parent_ptr->RefEl() == lf::base::RefEl::kTria()) {
          // Current cell was created  by splitting a triangle
          // Inheritance rules for refinement edges apply
          // Together with refinement pattern allows identification of child
          // cell
          const sub_idx_t fine_cell_child_number =
              (parent_infos_.back())[0][cell_index].child_number;
          const glb_idx_t parent_index = parent_mesh.Index(*parent_ptr);
          LF_VERIFY_MSG(parent_index < parent_mesh.NumEntities(0),
                        "parent_index = " << parent_index << " out of range");

          CONTROLLEDSTATEMENT(
              output_ctrl_, 100,
              std::cout << "Cell " << cell_index << ": " << std::flush;
              std::cout << "triangle child " << fine_cell_child_number
                        << " of parent " << parent_index << std::endl;)

          const CellChildInfo &parent_ci(parent_cell_ci[parent_index]);
          LF_VERIFY_MSG(
              parent_ci.child_cell_idx[fine_cell_child_number] == cell_index,
              "Parent child index mismatch!");
          const RefPat parent_ref_pat = parent_ci.ref_pat_;

          CONTROLLEDSTATEMENT(
              output_ctrl_, 100,
              std::cout << "ref_pat = " << (int)parent_ref_pat << std::flush;)

          switch (parent_ref_pat) {
            case RefPat::rp_nil: {
              LF_VERIFY_MSG(false,
                            "Parent cannot carry nil refinement pattern");
              break;
            }
            case RefPat::rp_copy: {
              CONTROLLEDSTATEMENT(output_ctrl_, 100,
                                  std::cout << "COPY" << std::flush;)
              // Inherit refinement edge from parent triangle
              child_ref_edges[cell_index] = parent_ref_edges[parent_index];
              break;
            }
            case RefPat::rp_bisect: {
              CONTROLLEDSTATEMENT(output_ctrl_, 100,
                                  std::cout << "BISECT" << std::flush;)
              // Both children have refinement edge 2
              LF_VERIFY_MSG(fine_cell_child_number < 2,
                            "Only 2 children for rp_bisect");
              child_ref_edges[cell_index] = 2;
              break;
            }
            case RefPat::rp_trisect: {
              CONTROLLEDSTATEMENT(output_ctrl_, 100,
                                  std::cout << "TRISECT" << std::flush;)
              // Refinement edges: 0 -> 2, 1 -> 1, 2 -> 0
              LF_VERIFY_MSG(fine_cell_child_number < 3,
                            "Only 3 children for rp_trisect");
              switch (fine_cell_child_number) {
                case 0: {
                  child_ref_edges[cell_index] = 2;
                  break;
                }
                case 1: {
                  child_ref_edges[cell_index] = 1;
                  break;
                }
                case 2: {
                  child_ref_edges[cell_index] = 0;
                  break;
                }
              }
              break;
            }
            case RefPat::rp_trisect_left: {
              CONTROLLEDSTATEMENT(output_ctrl_, 100,
                                  std::cout << "TRISECT_LEFT" << std::flush;)
              // Refinement edges: 0 -> 2, 1 -> 1, 2 -> 0
              LF_VERIFY_MSG(fine_cell_child_number < 4,
                            "Only 3 children for rp_quadsect");
              switch (fine_cell_child_number) {
                case 0: {
                  child_ref_edges[cell_index] = 2;
                  break;
                }
                case 1: {
                  child_ref_edges[cell_index] = 2;
                  break;
                }
                case 2: {
                  child_ref_edges[cell_index] = 0;
                  break;
                }
              }
              break;
            }
            case RefPat::rp_quadsect: {
              CONTROLLEDSTATEMENT(output_ctrl_, 100,
                                  std::cout << "QUADSECT" << std::flush;)
              // Refinement edges: 0 -> 2, 1 -> 0, 2 -> 0, 3-> 0
              switch (fine_cell_child_number) {
                case 0: {
                  child_ref_edges[cell_index] = 2;
                  break;
                }
                case 1: {
                  child_ref_edges[cell_index] = 0;
                  break;
                }
                case 2: {
                  child_ref_edges[cell_index] = 0;
                  break;
                }
                case 3: {
                  child_ref_edges[cell_index] = 0;
                  break;
                }
                default: {
                  LF_VERIFY_MSG(false, "Illegal child number");
                  break;
                }
              }
              break;
            }
            case rp_regular: {
              CONTROLLEDSTATEMENT(output_ctrl_, 100,
                                  std::cout << "REGULAR" << std::flush;)
              // Inherit the refinement edge of the parent triangle
              const sub_idx_t parent_ref_edge_idx =
                  parent_ref_edges[parent_index];
              switch (parent_ref_edge_idx) {
                case 0: {
                  switch (fine_cell_child_number) {
                    case 0: {
                      child_ref_edges[cell_index] = 0;
                      break;
                    }
                    case 1: {
                      child_ref_edges[cell_index] = 0;
                      break;
                    }
                    case 2: {
                      child_ref_edges[cell_index] = 1;
                      break;
                    }
                    case 3: {
                      child_ref_edges[cell_index] = 1;
                      break;
                    }
                    default: {
                      LF_VERIFY_MSG(false, "Illegal child number");
                      break;
                    }
                  }
                  break;
                }
                case 1: {
                  switch (fine_cell_child_number) {
                    case 0: {
                      child_ref_edges[cell_index] = 1;
                      break;
                    }
                    case 1: {
                      child_ref_edges[cell_index] = 2;
                      break;
                    }
                    case 2: {
                      child_ref_edges[cell_index] = 2;
                      break;
                    }
                    case 3: {
                      child_ref_edges[cell_index] = 2;
                      break;
                    }
                    default: {
                      LF_VERIFY_MSG(false, "Illegal child number");
                      break;
                    }
                  }
                  break;
                }
                case 2: {
                  switch (fine_cell_child_number) {
                    case 0: {
                      child_ref_edges[cell_index] = 2;
                      break;
                    }
                    case 1: {
                      child_ref_edges[cell_index] = 1;
                      break;
                    }
                    case 2: {
                      child_ref_edges[cell_index] = 0;
                      break;
                    }
                    case 3: {
                      child_ref_edges[cell_index] = 0;
                      break;
                    }
                    default: {
                      LF_VERIFY_MSG(false, "Illegal child number");
                      break;
                    }
                  }
                  break;
                }
              }  // end switch parent ref_edge_idx
              break;
            }
            case rp_barycentric: {
              CONTROLLEDSTATEMENT(output_ctrl_, 100,
                                  std::cout << "BARYCENTRIC" << std::flush;)
              // In the case of barycentric refinement choose the longest edge
              // as refinement edge for every child triangle
              child_ref_edges[cell_index] = LongestEdge(fine_cell);
              break;
            }
            default: {
              LF_VERIFY_MSG(false, "Illegal refinement type for a triangle");
              break;
            }
          }  // end switch parent_ref_pat
          CONTROLLEDSTATEMENT(
              output_ctrl_, 100,
              std::cout << " ref edge = " << child_ref_edges[cell_index]
                        << std::endl;)

        }  // end treatment of triangular child cell
        else if (parent_ptr->RefEl() == lf::base::RefEl::kQuad()) {
          // Parent is a quadrilateral:
          // refinement edge will be set to the longest edge
          child_ref_edges[cell_index] = LongestEdge(fine_cell);
        } else {
          LF_VERIFY_MSG(false, "Unknown parent cell type");
        }
      }
    }
  }
}

sub_idx_t MeshHierarchy::LongestEdge(const lf::mesh::Entity &T) const {
  LF_VERIFY_MSG(T.Codim() == 0, "Entity must be a call");
  // Obtain iterator over the edges
  const size_type num_edges = T.RefEl().NumSubEntities(1);
  base::RandomAccessRange<const lf::mesh::Entity> sub_edges(T.SubEntities(1));
  double max_len = 0.0;
  sub_idx_t idx_longest_edge = 0;
  Eigen::MatrixXd mp_refc(1, 1);
  mp_refc(0, 0) = 0.5;
  for (int k = 0; k < num_edges; k++) {
    // Approximate length by "1-point quadrature"
    const double approx_length =
        (sub_edges[k].Geometry()->IntegrationElement(mp_refc))[0];
    if (max_len < approx_length) {
      idx_longest_edge = k;
      max_len = approx_length;
    }
  }
  return idx_longest_edge;
}

std::ostream &MeshHierarchy::PrintInfo(std::ostream &o) const {
  o << "MeshHierarchy, " << NumLevels() << " levels: " << std::endl;
  for (int level = 0; level < NumLevels(); ++level) {
    const lf::mesh::Mesh &mesh{*getMesh(level)};
    o << "l=" << level << ": ";
    if ((ctrl_ & kout_meshinfo) != 0) {
      LF_ASSERT_MSG(false, "Not yet implemented");
      // TODO(raffael), when output for lf::mesh::Mesh has been fixed
      o << mesh << std::endl;
    } else {
      o << static_cast<int>(mesh.DimMesh()) << "D -> "
        << static_cast<int>(mesh.DimWorld()) << "D, " << mesh.NumEntities(0)
        << " cells [" << mesh.NumEntities(lf::base::RefEl::kTria()) << " tria, "
        << mesh.NumEntities(lf::base::RefEl::kQuad()) << " quad], "
        << mesh.NumEntities(1) << " edges, " << mesh.NumEntities(2) << " nodes"
        << std::endl;
    }
  }
  return o;
}

// Utility function for generating a hierarchy of meshes
/* SAM_LISTING_BEGIN_1 */
std::shared_ptr<MeshHierarchy> GenerateMeshHierarchyByUniformRefinemnt(
    std::shared_ptr<lf::mesh::Mesh> mesh_p, lf::base::size_type ref_lev,
    RefPat ref_pat) {
  LF_ASSERT_MSG(mesh_p != nullptr, "No valid mesh supplied!");
  // Set up the builder object for mesh entities, here suitable for a 2D hybrid
  // mesh comprising triangles and quadrilaterals
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  // Create a mesh hierarchy with a single level
  std::shared_ptr<MeshHierarchy> multi_mesh_p =
      std::make_shared<MeshHierarchy>(std::move(mesh_p), mesh_factory_ptr);
  // Perform the desired number of steps of uniform refinement
  for (int refstep = 0; refstep < ref_lev; ++refstep) {
    // Conduct regular refinement of all cells of the currently finest mesh.
    // This adds another mesh to the sequence of meshes. 
    multi_mesh_p->RefineRegular(ref_pat);
  }
  return multi_mesh_p;
}
/* SAM_LISTING_END_1 */

}  // namespace lf::refinement
