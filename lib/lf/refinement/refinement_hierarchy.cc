/**
 * @file refinement_hierarchy.cc
 * @brief implementation of global/local refinement methods
 */

#include "refinement_hierarchy.h"

namespace lf::refinement {

  // Implementation of MeshHierarchy
  MeshHierarchy::MeshHierarchy(std::shared_ptr<mesh::Mesh> base_mesh,
			       mesh::MeshFactory &mesh_factory)
    : mesh_factory_(mesh_factory) {
    LF_VERIFY_MSG(base_mesh, "No valid mesh supplied");
    LF_VERIFY_MSG(base_mesh->DimMesh() == 2, "Implemented only for 2D meshes");
    // Set coarsest mesh
    meshes_.push_back(base_mesh);
    {
      // Refinement edge has to be set for edges
      refinement_edges_.push_back
	(std::vector<lf::base::sub_idx_t>(base_mesh->Size(0), idx_nil));
      for (const mesh::Entity &cell : base_mesh->Entities(0)) {
	lf::base::glb_idx_t cell_index = base_mesh->Index(cell);
	if (cell.RefEl() == lf::base::RefEl::kTria()) {
	  (refinement_edges_.back())[cell_index] = LongestEdge(cell);
	}
      } // end loop over cells
      
      // Setting up child information
      std::vector<CellChildInfo> cell_child_info(base_mesh->Size(0));
      std::vector<EdgeChildInfo> edge_child_info(base_mesh->Size(1));
      std::vector<PointChildInfo> point_child_info(base_mesh->Size(2));
      cell_child_infos_.push_back(std::move(cell_child_info));
      edge_child_infos_.push_back(std::move(edge_child_info));
      point_child_infos_.push_back(std::move(point_child_info));
    }
    {
      // No parents for entities on the coarsest level
      std::vector<ParentInfo> cell_parent_info(base_mesh->Size(0));
      std::vector<ParentInfo> edge_parent_info(base_mesh->Size(1));
      std::vector<ParentInfo> point_parent_info(base_mesh->Size(2));
      parent_infos_.push_back({std::move(cell_parent_info),
	    std::move(edge_parent_info),
	    std::move(point_parent_info)});
    }
  }

  void MeshHierarchy::RefineRegular(RefPat ref_pat) {
    LF_VERIFY_MSG(ref_pat == RefPat::rp_regular ||
		  ref_pat == RefPat::rp_barycentric,
		  "Only regular or barycentric uniform refinement possible");
    // Retrieve the finest mesh in the hierarchy
    const mesh::Mesh &finest_mesh(*meshes_.back());
    // Flag all points as to be copied
    for (const mesh::Entity &point : finest_mesh.Entities(2)) {
      const lf::base::glb_idx_t point_index = finest_mesh.Index(point);
      PointChildInfo &pt_child_info((point_child_infos_.back())[point_index]);
      // Set the information about a points children except the child pointer
      pt_child_info.ref_pat_ = RefPat::rp_copy;
    }
    // Flag all edges as to be split
    for (const mesh::Entity &edge : finest_mesh.Entities(1)) {
      const lf::base::glb_idx_t edge_index = finest_mesh.Index(edge);
      EdgeChildInfo &ed_child_info(((edge_child_infos_.back()))[edge_index]);
      ed_child_info.ref_pat_ = RefPat::rp_split;
    }
    // All cells are to be refined regularly
    for (const mesh::Entity &cell : finest_mesh.Entities(0)) {
      const lf::base::glb_idx_t cell_index = finest_mesh.Index(cell);
      CellChildInfo &cell_child_info(((cell_child_infos_.back()))[cell_index]);
      cell_child_info.ref_pat_ = ref_pat;
    }
    // With all refinement patterns set, generate the new mesh
    PerformRefinement();
  }

  template <typename Marker>
  void MeshHierarchy::MarkEdges(Marker &&marker) {
    // Retrieve the finest mesh in the hierarchy
    const mesh::Mesh &finest_mesh(*meshes_.back());
    // Run through the edges = entities of co-dimension 1
    for (const mesh::Entity &edge : finest_mesh.Entities(1)) {
      lf::base::glb_idx_t edge_index = finest_mesh.Index(edge);
      (edge_marked_.back())[edge_index] = marker(finest_mesh, edge);
    }
  }

  void MeshHierarchy::RefineMarked(void) {
    // Target the finest mesh
    const lf::mesh::Mesh &finest_mesh(*meshes_.back());

    // Flag all points as to be copied
    for (const mesh::Entity &point : finest_mesh.Entities(2)) {
      const glb_idx_t point_index = finest_mesh.Index(point);
      PointChildInfo &pt_child_info((point_child_infos_.back())[point_index]);
      // Set the information about a points children except the child pointer
      pt_child_info.ref_pat_ = RefPat::rp_copy;
    }
    
    // Set the split refinement patterrn for all marked edges
    for (const lf::mesh::Entity &edge : finest_mesh.Entities(1)) {
      LF_VERIFY_MSG(edge.RefEl() == lf::base::RefEl::kSegment(),
		    "Wrong type for an edge");
      const glb_idx_t edge_index = finest_mesh.Index(edge);
      EdgeChildInfo &ed_ci(edge_child_infos_.back()[edge_index]);
      LF_VERIFY_MSG(ed_ci.ref_pat_ == RefPat::rp_nil,
		    "Edge " << edge_index << " already refined!");
      if ((edge_marked_.back())[edge_index]) {
	// Edge is to be refined
	ed_ci.ref_pat_ = RefPat::rp_split;
      }
      else {
	// Just copy edge
	ed_ci.ref_pat_ = RefPat::rp_copy;
      }
    }
    // Now all edges are initially marked to be split or copied

    // To keep the mesh conforming refinement might have to propagate
    bool refinement_complete;
    do {
      refinement_complete = true;
      // Visit all cells and update their refinement patterns
      for (const lf::mesh::Entity &cell : finest_mesh.Entities(0)) {
	const glb_idx_t cell_index = finest_mesh.Index(cell);
	CellChildInfo &cell_child_info(cell_child_infos_.back()[cell_index]);
	// Find edges which are marked as split
	std::array<bool,4> edge_split;
	// Array of references to edge sub-entities of current cell
	base::RandomAccessRange<const lf::mesh::Entity> sub_edges(cell.SubEntities(1));
	const size_type num_edges = cell.RefEl().NumSubEntities(1);
	LF_VERIFY_MSG(num_edges < 4,"Too many edges = " << num_edges);
	size_type split_edge_cnt = 0;
	for (int k=0; k < num_edges; k++) {
	  const glb_idx_t edge_index = finest_mesh.Index(sub_edges[k]);
	  edge_split[k] =
	    (edge_child_infos_.back()[edge_index].ref_pat_ == RefPat::rp_split);
	  if (edge_split[k]) split_edge_cnt++;
	}
	const sub_idx_t ref_edge_idx = refinement_edges_.back()[cell_index];
	switch (cell.RefEl()) {
	case lf::base::RefEl::kTria(): {
	  
	  break;
	}
	case lf::base::RefEl::kQuad(): {
	  
	  break;
	}
	default: {
	  LF_VERIFY_MSG(false,"Illegal cell type");
	  break;
	}
	} // end switch cell type
      } // end loop over cells
    }
    while (!refinement_complete);

    PerformRefinement();
  } // end RefineMarked

  void MeshHierarchy::PerformRefinement(void) {
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
    std::vector<PointChildInfo> &pt_child_info(point_child_infos_.back());
    const Hybrid2DRefinementPattern rp_copy_node(lf::base::RefEl::kPoint(),RefPat::rp_copy);
    for (const mesh::Entity &node : parent_mesh.Entities(2)) {
      // Obtain index of node in coarse mesh
      const lf::base::glb_idx_t node_index = parent_mesh.Index(node);
      LF_VERIFY_MSG(node_index < pt_child_info.size(),
		    "Node index " << node_index << " out of range");
      // Find position of node in physical coordinates
      const lf::geometry::Geometry &pt_geo(*node.Geometry());
      if (pt_child_info[node_index].ref_pat_ != RefPat::rp_nil) {
	// Generate a node for the fine mesh at the same position
	std::vector<std::unique_ptr<geometry::Geometry>> pt_child_geo_ptrs(pt_geo.ChildGeometry(rp_copy_node,0));
	LF_VERIFY_MSG(pt_child_geo_ptrs.size() == 1,
		      "A point can only have one chile");
	pt_child_info[node_index].child_point_idx_
	  = mesh_factory_.AddPoint(std::move(pt_child_geo_ptrs[0]));
      }
    } // end loop over nodes 

    // Now traverse the edges. Depending on the refinement pattern,
    // either copy them or split them.
    // Supplement the refinement information for edges accordingly.
    std::vector<EdgeChildInfo> &ed_child_info(edge_child_infos_.back());
    for (const mesh::Entity &edge : parent_mesh.Entities(1)) {
      // Fetch global index of edge
      lf::base::glb_idx_t edge_index = parent_mesh.Index(edge);
    
      // Get indices of endpoints in parent mesh
      auto ed_nodes(edge.SubEntities(1));
      const lf::base::glb_idx_t ed_p0_coarse_idx = parent_mesh.Index(ed_nodes[0]);
      const lf::base::glb_idx_t ed_p1_coarse_idx = parent_mesh.Index(ed_nodes[1]);
    
      // Obtain indices of the nodes at the same position in the fine mesh
      const lf::base::glb_idx_t ed_p0_fine_idx =
	pt_child_info[ed_p0_coarse_idx].child_point_idx_;
      const lf::base::glb_idx_t ed_p1_fine_idx =
	pt_child_info[ed_p1_coarse_idx].child_point_idx_;
      // Prepare request of geometry after refinement
      EdgeChildInfo &edge_ci(ed_child_info[edge_index]);
      const RefPat edge_refpat(edge_ci.ref_pat_);
      Hybrid2DRefinementPattern rp(edge.RefEl(),edge_refpat);
    
      // Distinguish between different local refinement patterns
      switch (edge_refpat) {
      case RefPat::rp_copy: {
        // Edge has to be duplicated
        std::vector<std::unique_ptr<lf::geometry::Geometry>>
	  ed_copy(edge.Geometry()->ChildGeometry(rp, 0));
        LF_VERIFY_MSG(ed_copy.size() == 1, "Copy may create only a single child!");
	// Register the new edge
	edge_ci.child_edge_idx_.push_back
	  (mesh_factory_.AddEntity
	   (edge.RefEl(),lf::base::ForwardRange<const lf::base::glb_idx_t>
	    ({ed_p0_fine_idx,ed_p1_fine_idx}),std::move(ed_copy[0])));
        break;
      }  // end rp_copy
      case rp_split: {
	// Edge has to be split into two halves with a new node created at the midpoint
	// position.
	// First obtain geometry information about new node (sub-entity with relative
	// co-dimension 1)
	std::vector<std::unique_ptr<geometry::Geometry>>
	  edge_nodes_geo_ptrs(edge.Geometry()->ChildGeometry(rp,1));
	LF_VERIFY_MSG(edge_nodes_geo_ptrs.size() == 1,
		      "Split edge with " << edge_nodes_geo_ptrs.size() << " child nodes!");
	// Register midpoint as new node
	const lf::base::glb_idx_t midpoint_fine_idx =
	  mesh_factory_.AddPoint(std::move(edge_nodes_geo_ptrs[0]));
	edge_ci.child_point_idx_.push_back(midpoint_fine_idx);
	// Next get the geometry objects for the two child edges (co-dim == 0)
	std::vector<std::unique_ptr<geometry::Geometry>>
	  edge_child_geo_ptrs(edge.Geometry()->ChildGeometry(rp,0));
	LF_VERIFY_MSG(edge_child_geo_ptrs.size() == 2,
		      "Split edge with " << edge_child_geo_ptrs.size() << " parts!");
	// Register the two new edges
	// CAREFUL: Assignment of endpoints has to match implementation in
	// refinement.cc
	edge_ci.child_edge_idx_.push_back
	  (mesh_factory_.AddEntity
	   (edge.RefEl(),lf::base::ForwardRange<const lf::base::glb_idx_t>
	    ({ed_p0_fine_idx,midpoint_fine_idx}),std::move(edge_child_geo_ptrs[0])));
	edge_ci.child_edge_idx_.push_back
	  (mesh_factory_.AddEntity
	   (edge.RefEl(),lf::base::ForwardRange<const lf::base::glb_idx_t>
	    ({midpoint_fine_idx,ed_p1_fine_idx}),std::move(edge_child_geo_ptrs[1])));
	break;
      } // end rp_split
      default: {
	LF_VERIFY_MSG(false,"Refinement pattern " << (int)edge_refpat << " illegal for edge");
	break;
      }
      }    // end switch refpat
    }      // end edge loop

    // Visit all cells, examine their refinement patterns, retrieve indices of
    // their sub-entities, and those of the children.
    std::vector<CellChildInfo> &cell_child_info(cell_child_infos_.back());
    for(const mesh::Entity &cell : parent_mesh.Entities(0)) {
      // type of cell
      const lf::base::RefEl ref_el(cell.RefEl());
      const lf::base::size_type num_edges = ref_el.NumSubEntities(1);
      const lf::base::size_type num_vertices = ref_el.NumSubEntities(2);
      // fetch index of current cell
      const lf::base::glb_idx_t cell_index(parent_mesh.Index(cell));
    
      // Obtain indices of subentities (co-dimension = outer array index)
      std::array<std::vector<lf::base::glb_idx_t>,3> cell_subent_idx;
      cell_subent_idx[0].push_back(cell_index);
      for (int codim = 1; codim <= 2; codim++) {
	base::RandomAccessRange<const mesh::Entity> subentities(cell.SubEntities(codim));
	for (const mesh::Entity &sub_ent : subentities) {
	  cell_subent_idx[codim].push_back(parent_mesh.Index(sub_ent));
	}
	LF_VERIFY_MSG(cell_subent_idx[codim].size() == ref_el.NumSubEntities(codim),
		      ref_el.ToString() << ": only " << cell_subent_idx[codim].size()
		      << " subents of codim = " << codim);
      }
      // Index information for sub-entities with respect  to fine mesh
      // Retrieve indices of vertices of cell on the fine mesh
      std::array<lf::base::glb_idx_t,4> vertex_child_idx({idx_nil,idx_nil,idx_nil,idx_nil});
      for (lf::base::sub_idx_t vt_lidx = 0; vt_lidx < num_vertices; vt_lidx++) {
	LF_VERIFY_MSG(pt_child_info[cell_subent_idx[2][vt_lidx]].ref_pat_ == RefPat::rp_copy,
		      "Vertex must have been copied!");
	vertex_child_idx[vt_lidx] = pt_child_info[cell_subent_idx[2][vt_lidx]].child_point_idx_;
      }
      // Retrieve indices of midpoints of edges, if they exist
      std::array<lf::base::glb_idx_t,4> edge_midpoint_idx({idx_nil,idx_nil,idx_nil,idx_nil});
      for (lf::base::sub_idx_t ed_lidx = 0; ed_lidx < num_edges; ed_lidx++) {
	const EdgeChildInfo &ed_ci(ed_child_info[cell_subent_idx[1][ed_lidx]]);
	if (ed_ci.child_point_idx_.size() > 0) {
	  LF_VERIFY_MSG(ed_child_info[cell_subent_idx[1][ed_lidx]].ref_pat_ == RefPat::rp_split,
			"Edge with a midpoint must have been split");
	  edge_midpoint_idx[ed_lidx] = ed_ci.child_point_idx_[0];
	}
      } // end loop over local edges

      // Set up refinement object -> v
      CellChildInfo &cell_ci(cell_child_info[cell_index]);
      const RefPat cell_refpat(cell_ci.ref_pat_);
      const sub_idx_t anchor = cell_ci.anchor_; // anchor edge index
      Hybrid2DRefinementPattern rp(cell.RefEl(),cell_refpat,anchor);

      // Index offsets for refinement patterns requiring an ancchor edge
      std::array<sub_idx_t,4> mod;
      if (anchor != idx_nil) {
	for (int k =0; k < num_vertices; k++) mod[k] = (k + anchor) % num_vertices;
      }

      // Array of node indices (w.r.t. fine mesh) for sub-cells (triangles or quadrilaterals)
      std::vector<std::vector<glb_idx_t>> child_cell_nodes;
      std::vector<glb_idx_t> tria_ccn_tmp(3);
      std::vector<glb_idx_t> quad_ccn_tmp(4);
      // Array of node indices for interior child edges
      std::vector<std::array<glb_idx_t,2>> child_edge_nodes;
      std::array<glb_idx_t,2> cen_tmp;
    
      // Local linking of child entities is very different depending on cell type
      // and refinement patterns.
      if (ref_el == lf::base::RefEl::kTria()) {
	// Case of a triangle
	switch (cell_refpat) {
	case RefPat::rp_nil: {
	  LF_VERIFY_MSG(false,"Every triangle must be refined");
	  break;
	}
	case RefPat::rp_copy: {
	  // Just copy the parent triangle
	  child_cell_nodes.push_back
	    ({vertex_child_idx[0],vertex_child_idx[1],vertex_child_idx[2]});
	  break;
	}
	case RefPat::rp_bisect: {
	  LF_VERIFY_MSG(anchor != idx_nil,
			"Anchor must be set for bisection refinement of triangle");
	  // Splitting a triangle in two by bisecting anchor edge
	  tria_ccn_tmp[0] = vertex_child_idx[mod[0]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = vertex_child_idx[mod[2]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = vertex_child_idx[mod[2]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

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
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[2]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[2]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

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
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = vertex_child_idx[mod[2]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[2]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[2]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

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
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[2]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[2]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[2]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

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
	  //  Split triangle into 5 smaller triangles by connecting
	  // the center of gravity with the vertices and the midpoints
	  // of the edges.
	  // Create a new interior vertex
	  std::vector<std::unique_ptr<geometry::Geometry>>
	    cell_center_geo_ptrs(cell.Geometry()->ChildGeometry(rp,2)); // point: co-dim == 2
	  LF_VERIFY_MSG(cell_center_geo_ptrs.size() == 1,
			"Barycentrically refined triangle with "
			<< cell_center_geo_ptrs.size() << " interior child nodes ??");
	  // Register midpoint as new node
	  const glb_idx_t center_fine_idx  =
	    mesh_factory_.AddPoint(std::move(cell_center_geo_ptrs[0]));
	  cell_ci.child_point_idx_.push_back(center_fine_idx);

 	  tria_ccn_tmp[0] = vertex_child_idx[0];
	  tria_ccn_tmp[1] = edge_midpoint_idx[0];
	  tria_ccn_tmp[2] = center_fine_idx;
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[1];
	  tria_ccn_tmp[1] = edge_midpoint_idx[0];
	  tria_ccn_tmp[2] = center_fine_idx;
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[1];
	  tria_ccn_tmp[1] = edge_midpoint_idx[1];
	  tria_ccn_tmp[2] = center_fine_idx;
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[2];
	  tria_ccn_tmp[1] = edge_midpoint_idx[1];
	  tria_ccn_tmp[2] = center_fine_idx;
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[2];
	  tria_ccn_tmp[1] = edge_midpoint_idx[2];
	  tria_ccn_tmp[2] = center_fine_idx;
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[0];
	  tria_ccn_tmp[1] = edge_midpoint_idx[2];
	  tria_ccn_tmp[2] = center_fine_idx;
	  child_cell_nodes.push_back(tria_ccn_tmp);

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
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[1];
	  tria_ccn_tmp[1] = edge_midpoint_idx[0];
	  tria_ccn_tmp[2] = edge_midpoint_idx[1];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[2];
	  tria_ccn_tmp[1] = edge_midpoint_idx[2];
	  tria_ccn_tmp[2] = edge_midpoint_idx[1];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = edge_midpoint_idx[0];
	  tria_ccn_tmp[1] = edge_midpoint_idx[1];
	  tria_ccn_tmp[2] = edge_midpoint_idx[2];
	  child_cell_nodes.push_back(tria_ccn_tmp);

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
	  LF_VERIFY_MSG(false,"Refinement pattern not (yet) implemented for triangle");
	  break;
	}
	} // end switch refpat
      }
      else if (ref_el == lf::base::RefEl::kQuad()) {
	// Case of a quadrilateral
	switch (cell_refpat) {
	case RefPat::rp_nil: {
	  LF_VERIFY_MSG(false,"Every quadrilateral must be refined");
	  break;
	}
	case RefPat::rp_copy: {
	  // Just copy the parent triangle
	  child_cell_nodes.push_back
	    ({vertex_child_idx[0],vertex_child_idx[1],vertex_child_idx[2],vertex_child_idx[3]});
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
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[1] = vertex_child_idx[mod[0]];
	  tria_ccn_tmp[2] = vertex_child_idx[mod[3]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[1] = vertex_child_idx[mod[1]];
	  tria_ccn_tmp[2] = vertex_child_idx[mod[2]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

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
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[1]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[0]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[2]];
	  tria_ccn_tmp[1] = vertex_child_idx[mod[3]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[1]];
	  tria_ccn_tmp[2] = vertex_child_idx[mod[3]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

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
	  child_cell_nodes.push_back(quad_ccn_tmp);

	  quad_ccn_tmp[0] = vertex_child_idx[mod[1]];
	  quad_ccn_tmp[1] = vertex_child_idx[mod[2]];
	  quad_ccn_tmp[2] = edge_midpoint_idx[mod[2]];
	  quad_ccn_tmp[3] = edge_midpoint_idx[mod[0]];
	  child_cell_nodes.push_back(quad_ccn_tmp);

	  // edges
	  cen_tmp[0] = edge_midpoint_idx[mod[0]];
	  cen_tmp[1] = edge_midpoint_idx[mod[2]];
	  child_edge_nodes.push_back(cen_tmp);
	  break;
	}
	case RefPat::rp_threeedge: {
	  LF_VERIFY_MSG(anchor != idx_nil,
			"Anchor must be set for three edge refinement of a quad");

	  quad_ccn_tmp[0] = vertex_child_idx[mod[2]];
	  quad_ccn_tmp[1] = vertex_child_idx[mod[3]];
	  quad_ccn_tmp[2] = edge_midpoint_idx[mod[3]];
	  quad_ccn_tmp[3] = edge_midpoint_idx[mod[1]];
	  child_cell_nodes.push_back(quad_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[0]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[3]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = vertex_child_idx[mod[1]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[1]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

	  tria_ccn_tmp[0] = edge_midpoint_idx[mod[0]];
	  tria_ccn_tmp[1] = edge_midpoint_idx[mod[1]];
	  tria_ccn_tmp[2] = edge_midpoint_idx[mod[3]];
	  child_cell_nodes.push_back(tria_ccn_tmp);

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
	    cell_center_geo_ptrs(cell.Geometry()->ChildGeometry(rp,2)); // point: co-dim == 2
	  LF_VERIFY_MSG(cell_center_geo_ptrs.size() == 1,
			"Regularly refined quadrilateral with "
			<< cell_center_geo_ptrs.size() << " interior child nodes!");
	  // Register midpoint as new node
	  const glb_idx_t center_fine_idx  =
	    mesh_factory_.AddPoint(std::move(cell_center_geo_ptrs[0]));
	  cell_ci.child_point_idx_.push_back(center_fine_idx);

	  // Set the node indices (w.r.t. fine mesh) of the four sub-quads
	  quad_ccn_tmp[0] = vertex_child_idx[0];
	  quad_ccn_tmp[1] = edge_midpoint_idx[0];
	  quad_ccn_tmp[2] = center_fine_idx;
	  quad_ccn_tmp[3] = edge_midpoint_idx[3];
	  child_cell_nodes.push_back(quad_ccn_tmp);

	  quad_ccn_tmp[0] = vertex_child_idx[1];
	  quad_ccn_tmp[1] = edge_midpoint_idx[1];
	  quad_ccn_tmp[2] = center_fine_idx;
	  quad_ccn_tmp[3] = edge_midpoint_idx[0];
	  child_cell_nodes.push_back(quad_ccn_tmp);

	  quad_ccn_tmp[0] = vertex_child_idx[2];
	  quad_ccn_tmp[1] = edge_midpoint_idx[1];
	  quad_ccn_tmp[2] = center_fine_idx;
	  quad_ccn_tmp[3] = edge_midpoint_idx[2];
	  child_cell_nodes.push_back(quad_ccn_tmp);

	  quad_ccn_tmp[0] = vertex_child_idx[3];
	  quad_ccn_tmp[1] = edge_midpoint_idx[2];
	  quad_ccn_tmp[2] = center_fine_idx;
	  quad_ccn_tmp[3] = edge_midpoint_idx[3];
	  child_cell_nodes.push_back(quad_ccn_tmp);

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
	  LF_VERIFY_MSG(false,"Refinement pattern not (yet) implemented for quadrilateral");
	  break;
	}
	}// end switch cell_refpat
      } // end quadrilateral case
      else {
	LF_VERIFY_MSG(false,"Unknown cell type" << ref_el.ToString());
      }
      // Register new edges
      {
	std::vector<std::unique_ptr<geometry::Geometry>>
	  cell_edge_geo_ptrs(cell.Geometry()->ChildGeometry(rp,1)); // child edge: co-dim == 2
	const size_type num_new_edges = child_edge_nodes.size();
	LF_VERIFY_MSG(num_new_edges == rp.noChildren(1),
		      "num_new_edges = " << num_new_edges << " <-> " << rp.noChildren(1));
	for (int l = 0; l < num_new_edges; l++) {
	  const std::array<glb_idx_t,2> &cen(child_edge_nodes[l]);
	  const glb_idx_t new_edge_index =
	    mesh_factory_.AddEntity(lf::base::RefEl::kSegment(),{cen[0],cen[1]},
				    std::move(cell_edge_geo_ptrs[l]));
	  cell_ci.child_edge_idx_.push_back(new_edge_index);
	} // end loop over new edges
      } // end register new edges
      // Register new cells
      {
	std::vector<std::unique_ptr<geometry::Geometry>>
	  childcell_geo_ptrs(cell.Geometry()->ChildGeometry(rp,0)); // child cell: co-dim == 0
	const size_type num_new_cells = child_cell_nodes.size();
	LF_VERIFY_MSG(num_new_cells == rp.noChildren(0),
		      "num_new_cells = " << num_new_cells << " <-> " << rp.noChildren(0));
	for (int l = 0; l < num_new_cells; l++) {
	  const std::vector<glb_idx_t> &ccn(child_cell_nodes[l]);
	  glb_idx_t new_cell_index;
	  if (ccn.size() == 3) {
	    // New cell is a triangle
	    new_cell_index =
	      mesh_factory_.AddEntity(lf::base::RefEl::kTria(),{ccn[0],ccn[1],ccn[2]},
				      std::move(childcell_geo_ptrs[l]));
	  }
	  else if (ccn.size() == 4) {
	    // New cell is a quadrilateral
	    new_cell_index =
	      mesh_factory_.AddEntity(lf::base::RefEl::kQuad(),{ccn[0],ccn[1],ccn[2],ccn[3]},
				      std::move(childcell_geo_ptrs[l]));
	  }
	  else {
	    LF_VERIFY_MSG(false,"Child cells must be either triangles or quads");
	  }
	  cell_ci.child_cell_idx_.push_back(new_cell_index);
	} // end loop over new cells
      } // end register new cells
    }// end loop over cells

    // At this point the MeshFactory has complete information to generate the new finest mesh
    meshes_.push_back(mesh_factory_.Build());
    mesh::Mesh &child_mesh(*meshes_.back());
    
    // Finally, we have to initialize the parent pointers for the entities of the newly
    // created finest mesh.
    {
      parent_infos_.push_back
	({ std::vector<ParentInfo>(child_mesh.Size(0)),
	    std::vector<ParentInfo>(child_mesh.Size(1)),
	    std::vector<ParentInfo>(child_mesh.Size(2))});
    
      // Visit all nodes of the parent mesh and retrieve their children by index
      for (const mesh::Entity &node : parent_mesh.Entities(2)) {
	// Obtain index of node in coarse mesh
	const glb_idx_t node_index = parent_mesh.Index(node);
	const glb_idx_t child_index = pt_child_info[node_index].child_point_idx_;
	if (child_index != idx_nil) {
	  LF_VERIFY_MSG(child_index < (parent_infos_.back())[2].size(),
			"index " << child_index << " illegal for child vertex");
	  ((parent_infos_.back()[2])[child_index]).child_number_ = 0;
	  ((parent_infos_.back()[2])[child_index]).parent_ptr_ = &node;
	}
      } // end loop over nodes 

      // Traverse edges and set the parent for their interior nodes and child edges
      for (const mesh::Entity &edge : parent_mesh.Entities(1)) {
	// Obtain index of edge in coarse mesh
	const glb_idx_t edge_index = parent_mesh.Index(edge);
	const EdgeChildInfo &ci_edge(ed_child_info[edge_index]);

	// Visit child edges
	const size_type num_child_edges = ci_edge.child_edge_idx_.size();
	for (int l = 0; l<num_child_edges; l++) {
	  const glb_idx_t edge_child_idx = ci_edge.child_edge_idx_[l];
	  LF_VERIFY_MSG(edge_child_idx < (parent_infos_.back())[1].size(),
			"index " << edge_child_idx << " illegal for child edge");
	  ((parent_infos_.back()[1])[edge_child_idx]).child_number_ = l;
	  ((parent_infos_.back()[1])[edge_child_idx]).parent_ptr_ = &edge;
	} // end loop over child edges

	// Visit child nodes
	const size_type num_child_nodes = ci_edge.child_point_idx_.size();
	for (int l = 0; l<num_child_nodes; l++) {
	  const glb_idx_t node_child_idx = ci_edge.child_point_idx_[l];
	  LF_VERIFY_MSG(node_child_idx < (parent_infos_.back())[2].size(),
			"index " << node_child_idx << " illegal for child point");
	  ((parent_infos_.back()[2])[node_child_idx]).child_number_ = l;
	  ((parent_infos_.back()[2])[node_child_idx]).parent_ptr_ = &edge;
	} // end loop over child points
      } // end loop over edges


      // Loop over cells
      for(const mesh::Entity &cell : parent_mesh.Entities(0)) {
	// fetch index of current cell
	const glb_idx_t cell_index(parent_mesh.Index(cell));
	// Get infformation about children
	CellChildInfo &cell_ci(cell_child_info[cell_index]);

	// Visit child cells
	const size_type num_child_cells = cell_ci.child_cell_idx_.size();
	for (int l=0; l<num_child_cells; l++) {
	  const glb_idx_t child_cell_idx = cell_ci.child_cell_idx_[l];
	  LF_VERIFY_MSG(child_cell_idx < (parent_infos_.back())[0].size(),
			"index " << child_cell_idx << " illegal for child cell");
	  ((parent_infos_.back()[0])[child_cell_idx]).child_number_ = l;
	  ((parent_infos_.back()[0])[child_cell_idx]).parent_ptr_ = &cell;
	}

	// Visit child edges
	const size_type num_child_edges = cell_ci.child_edge_idx_.size();
	for (int l = 0; l<num_child_edges; l++) {
	  const glb_idx_t edge_child_idx = cell_ci.child_edge_idx_[l];
	  LF_VERIFY_MSG(edge_child_idx < (parent_infos_.back())[1].size(),
			"index " << edge_child_idx << " illegal for child edge");
	  ((parent_infos_.back()[1])[edge_child_idx]).child_number_ = l;
	  ((parent_infos_.back()[1])[edge_child_idx]).parent_ptr_ = &cell;
	} // end loop over child edges

	// Visit child nodes
	const size_type num_child_nodes = cell_ci.child_point_idx_.size();
	for (int l = 0; l<num_child_nodes; l++) {
	  const glb_idx_t node_child_idx = cell_ci.child_point_idx_[l];
	  LF_VERIFY_MSG(node_child_idx < (parent_infos_.back())[2].size(),
			"index " << node_child_idx << " illegal for child point");
	  ((parent_infos_.back()[2])[node_child_idx]).child_number_ = l;
	  ((parent_infos_.back()[2])[node_child_idx]).parent_ptr_ = &cell;
	} // end loop over child points
      } // end loop over cells of parent mesh
    }
    
    // Initialize (empty) refinement information for newly created fine mesh
    // Same code as in the constructor of MeshHierarchy
    {
     std::vector<CellChildInfo> fine_cell_child_info(child_mesh.Size(0));
     std::vector<EdgeChildInfo> fine_edge_child_info(child_mesh.Size(1));
     std::vector<PointChildInfo> fine_point_child_info(child_mesh.Size(2));
     cell_child_infos_.push_back(std::move(fine_cell_child_info));
     edge_child_infos_.push_back(std::move(fine_edge_child_info));
     point_child_infos_.push_back(std::move(fine_point_child_info));
    }
    // Finally set refinement edges for fine mesh
    {
      std::vector<sub_idx_t> &parent_ref_edges(refinement_edges_.back());
      refinement_edges_.push_back(std::vector<sub_idx_t>(child_mesh.Size(0),idx_nil));
      // Traverse the cells of the fine mesh
      for (const mesh::Entity &fine_cell : child_mesh.Entities(0)) {
	// Refinement edge relevant for triangles onyl
	if (fine_cell.RefEl() == lf::base::RefEl::kTria()) {
	  const glb_idx_t cell_index = child_mesh.Index(fine_cell);
	  // pointer to cell whose refinement has created the current one
	  const mesh::Entity *parent_ptr = (parent_infos_.back())[0][cell_index].parent_ptr_;
	  LF_VERIFY_MSG(parent_ptr != nullptr,"Every cell on the fine mesh must have a parent!");
	 
	  if (parent_ptr->RefEl() == lf::base::RefEl::kTria()) {
	    // Current cell was created  by splitting a triangle
	    // Inheritance rules for refinement edges apply
	    // Together with refinement pattern allows identification of child cell
	    const sub_idx_t fine_cell_child_number = (parent_infos_.back())[0][cell_index].child_number_;
	    const glb_idx_t parent_index = parent_mesh.Index(*parent_ptr);
	    LF_VERIFY_MSG(parent_index < parent_mesh.Size(0),
			  "parent_index = " << parent_index << " out of range");
	    const CellChildInfo parent_ci(cell_child_info[parent_index]);
	    LF_VERIFY_MSG(parent_ci.child_cell_idx_[fine_cell_child_number] == cell_index,
			  "Parent child index mismatch!");
	    const RefPat parent_ref_pat = parent_ci.ref_pat_;

	    switch (parent_ref_pat) {
	    case RefPat::rp_nil: {
	      LF_VERIFY_MSG(false,"Parent cannot carry nil refinement pattern");
	      break;
	    }
	    case RefPat::rp_copy: {
	      // Inherit refinement edge from parent triangle
	      refinement_edges_.back()[cell_index] = parent_ref_edges[parent_index];
	      break;
	    }
	    case RefPat::rp_bisect: {
	      // Both children have refinement edge 2
	      LF_VERIFY_MSG(fine_cell_child_number < 2,"Only 2 children for rp_bisect");
	      refinement_edges_.back()[cell_index] = 2;
	      break;
	    }
	    case RefPat::rp_trisect: {
	      // Refinement edges: 0 -> 2, 1 -> 1, 2 -> 0
	      LF_VERIFY_MSG(fine_cell_child_number < 3,"Only 3 children for rp_trisect");
	      switch (fine_cell_child_number) {
	      case 0: { refinement_edges_.back()[cell_index] = 2; break; }
	      case 1: { refinement_edges_.back()[cell_index] = 1; break; }
	      case 2: { refinement_edges_.back()[cell_index] = 0; break; }
	      }
	      break;
	    }
	    case RefPat::rp_trisect_left: {
	      // Refinement edges: 0 -> 2, 1 -> 1, 2 -> 0
	      LF_VERIFY_MSG(fine_cell_child_number < 4,"Only 3 children for rp_quadsect");
	      switch (fine_cell_child_number) {
	      case 0: { refinement_edges_.back()[cell_index] = 2; break; }
	      case 1: { refinement_edges_.back()[cell_index] = 2; break; }
	      case 2: { refinement_edges_.back()[cell_index] = 0; break; }
	      }
	      break;
	    }
	    case RefPat::rp_quadsect: {
	      // Refinement edges: 0 -> 2, 1 -> 0, 2 -> 0, 3-> 0
	      switch (fine_cell_child_number) {
	      case 0: { refinement_edges_.back()[cell_index] = 2; break; }
	      case 1: { refinement_edges_.back()[cell_index] = 0; break; }
	      case 2: { refinement_edges_.back()[cell_index] = 0; break; }
	      case 3: { refinement_edges_.back()[cell_index] = 0; break; }
	      default: { LF_VERIFY_MSG(false,"Illegal child number"); break; }
	      }
	      break;
	    }
	    case rp_regular: {
	      // Inherit the refinement edge of the parent triangle
	      const sub_idx_t parent_ref_edge_idx = parent_ref_edges[parent_index];
	      switch (parent_ref_edge_idx)  {
	      case 0: {
		switch (fine_cell_child_number) {
		case 0: { refinement_edges_.back()[cell_index] = 0; break; }
		case 1: { refinement_edges_.back()[cell_index] = 0; break; }
		case 2: { refinement_edges_.back()[cell_index] = 1; break; }
		case 3: { refinement_edges_.back()[cell_index] = 1; break; }
		default: { LF_VERIFY_MSG(false,"Illegal child number"); break; }
		}
		break;
	      }
	      case 1: {
		switch (fine_cell_child_number) {
		case 0: { refinement_edges_.back()[cell_index] = 1; break; }
		case 1: { refinement_edges_.back()[cell_index] = 2; break; }
		case 2: { refinement_edges_.back()[cell_index] = 2; break; }
		case 3: { refinement_edges_.back()[cell_index] = 2; break; }
		default: { LF_VERIFY_MSG(false,"Illegal child number"); break; }
		}
		break;
	      }
	      case 2: {
		switch (fine_cell_child_number) {
		case 0: { refinement_edges_.back()[cell_index] = 2; break; }
		case 1: { refinement_edges_.back()[cell_index] = 1; break; }
		case 2: { refinement_edges_.back()[cell_index] = 0; break; }
		case 3: { refinement_edges_.back()[cell_index] = 0; break; }
		default: { LF_VERIFY_MSG(false,"Illegal child number"); break; }
		}
		break;
	      }
	      } // end switch parent ref_edge_idx
	      break;
	    }
	    case rp_barycentric: {
	      // In the case of barycentric refinement choose the longest edge as
	      // refinement edge for every child triangle
	      refinement_edges_.back()[cell_index] = LongestEdge(fine_cell);
	      break;
	    }
	    default: {
	      LF_VERIFY_MSG(false,"Illegal refinement type for a triangle");
	      break;
	    }
	    } // end switch parent_ref_pat
	  } // end treatment of triangular child cell
	  else if (parent_ptr->RefEl() == lf::base::RefEl::kQuad()) {
	    // Parent is a quadrilateral:
	    // refinement edge will be set to the longest edge
	    refinement_edges_.back()[cell_index] = LongestEdge(fine_cell);
	  }
	  else {
	    LF_VERIFY_MSG(false,"Unknown parent cell type");
	  }
	}
      }
    }
  }

  sub_idx_t MeshHierarchy::LongestEdge(const lf::mesh::Entity &T) const {
    LF_VERIFY_MSG(T.Codim() == 0,"Entity must be a call");
    // Obtain iterator over the edges
    const size_type num_edges = T.RefEl().NumSubEntities(1);
    base::RandomAccessRange<const lf::mesh::Entity> sub_edges(T.SubEntities(1));
    double max_len = 0.0;
    sub_idx_t idx_longest_edge = 0;
    Eigen::MatrixXd mp_refc(1,1); mp_refc(0,0) = 0.5;
    for (int k=0; k < num_edges; k++) {
      // Approximate length by "1-point quadrature"
      const double approx_length = (sub_edges[k].Geometry()->IntegrationElement(mp_refc))[0];
      if (max_len < approx_length) {
	idx_longest_edge = k;
	max_len = approx_length;
      }
    }
    return idx_longest_edge;
  }
  
}  // namespace lf::refinement
