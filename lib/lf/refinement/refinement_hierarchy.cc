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
      refinement_edges_.push_back(
				  std::vector<lf::base::sub_idx_t>(base_mesh->Size(0), -1));
      for (const mesh::Entity &cell : base_mesh->Entities(0)) {
	lf::base::glb_idx_t cell_index = base_mesh->Index(cell);
	if (cell.RefEl() == lf::base::RefEl::kTria()) {
	  // TODO: set refinement edge to longest edge
	  (refinement_edges_.back())[cell_index] = -1;
	}
      }
      // Setting upe child information
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

  void MeshHierarchy::RefineRegular(void) {
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
      cell_child_info.ref_pat_ = RefPat::rp_regular;
    }
    // With all refinement patterns set, generate the new mesh
    PerformRefinement();
  }

  template <typename Marker>
  void MeshHierarchy::MarkEdges(Marker &&marker) {
    // Retrieve the finest mesh in the hierarchy
    const mesh::Mesh &finest_mesh(*meshes_.back());
    // Run through the edges = entities of co-dimension 1
    for (mesh::Entity &edge : finest_mesh.Entities(1)) {
      lf::base::glb_idx_t edge_index = finest_mesh.Index(edge);
      (edge_marked_.back())[edge_index] = marker(finest_mesh, edge);
    }
  }

  void MeshHierarchy::RefineMarked(void) {}

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
    for (const mesh::Entity &node : parent_mesh.Entities(2)) {
      // Obtain index of node in coarse mesh
      const lf::base::glb_idx_t node_index = parent_mesh.Index(node);
      // Find position of node in physical coordinates
      const lf::geometry::Geometry &pt_geo(*node.Geometry());
      Eigen::MatrixXd pt_coords(pt_geo.Global(Eigen::Matrix<double, 0, 1>()));
      if (pt_child_info[node_index].ref_pat_ != RefPat::rp_nil) {
	// Generate a node for the fine mesh at the same position
	// ADJUST to new version of AddPoint() using ChildGeometry() !!!!!!!
	pt_child_info[node_index].child_point_idx_ = mesh_factory_.AddPoint(pt_coords);
      }} // end loop over nodes 

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
	LF_VERIFY_MSG(edge_nodes_geo_ptrs.size() != 1,
		      "Split edge with " << edge_nodes_geo_ptrs.size() << " child nodes!");
	// Register midpoint as new node
	const lf::base::glb_idx_t midpoint_fine_idx =
	  mesh_factory_.AddPoint(std::move(edge_nodes_geo_ptrs[0]));
	edge_ci.child_point_idx_.push_back(midpoint_fine_idx);
	// Next get the geometry objects for the two child edges (co-dim == 0)
	std::vector<std::unique_ptr<geometry::Geometry>>
	  edge_child_geo_ptrs(edge.Geometry()->ChildGeometry(rp,0));
	LF_VERIFY_MSG(edge_child_geo_ptrs.size() != 2,
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
      // std::array<std::vector<const mesh::Entity *>,3> cell_subent_ptrs; NOT NEEDED
      std::array<std::vector<lf::base::glb_idx_t>,3> cell_subent_idx;
      cell_subent_idx[0].push_back(cell_index);
      // cell_subent_ptrs[0].push_back(&cell); NOT NEEDED
      for (int codim = 1; codim <= 2; codim++) {
	base::RandomAccessRange<const mesh::Entity> subentities(cell.SubEntities(codim));
	for (const mesh::Entity &sub_ent : subentities) {
	  // cell_subent_ptrs[codim].push_back(&sub_ent); NOT NEEDED
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

      // Index offsets
      std::array<sub_idx_t,4> mod;
      if (anchor != idx_nil) {
	for (int l =0; l < num_vertices; l++) mod[l] = (l + anchor) % num_vertices;
      }

      // Array of node indices (w.r.t. fine mesh) for sub-cells (triangles or quadrilaterals)
      std::vector<std::vector<glb_idx_t>> child_cell_nodes;
      std::vector<glb_idx_t> ccn_tmp {}; 
      // Array of node indices for interior child edges
      std::vector<std::array<glb_idx_t,2>> child_edge_nodes;
      std::array<glb_idx_t,2> cen_tmp;
    
      // Local linking of child entities is very different depending on cell type
      // and refinement patterns.
      if (ref_el == lf::base::RefEl::kTria()) {
	ccn_tmp.resize(3); // Only triangular children
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
	case RefPat::rp_regular: {
	  // regular refinement into four congruent triangles
	  // Child cells
	  ccn_tmp[0] = vertex_child_idx[0];
	  ccn_tmp[1] = edge_midpoint_idx[0];
	  ccn_tmp[2] = edge_midpoint_idx[2];
	  child_cell_nodes.push_back(ccn_tmp);

	  ccn_tmp[0] = vertex_child_idx[1];
	  ccn_tmp[1] = edge_midpoint_idx[0];
	  ccn_tmp[2] = edge_midpoint_idx[1];
	  child_cell_nodes.push_back(ccn_tmp);

	  ccn_tmp[0] = vertex_child_idx[2];
	  ccn_tmp[1] = edge_midpoint_idx[2];
	  ccn_tmp[2] = edge_midpoint_idx[1];
	  child_cell_nodes.push_back(ccn_tmp);

	  ccn_tmp[0] = edge_midpoint_idx[0];
	  ccn_tmp[1] = edge_midpoint_idx[1];
	  ccn_tmp[2] = edge_midpoint_idx[2];
	  child_cell_nodes.push_back(ccn_tmp);

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
	case RefPat::rp_regular: {
	  // Create a new interior vertex
	  std::vector<std::unique_ptr<geometry::Geometry>>
	    cell_center_geo_ptrs(cell.Geometry()->ChildGeometry(rp,2)); // point: co-dim == 2
	  LF_VERIFY_MSG(cell_center_geo_ptrs.size() != 1,
			"Regularly refined quadrilateral with "
			<< cell_center_geo_ptrs.size() << " interior child nodes!");
	  // Register midpoint as new node
	  const glb_idx_t center_fine_idx  =
	    mesh_factory_.AddPoint(std::move(cell_center_geo_ptrs[0]));
	  cell_ci.child_point_idx_.push_back(center_fine_idx);

	  // Set the node indices (w.r.t. fine mesh) of the four sub-quads
	  ccn_tmp.resize(4);
	
	  ccn_tmp[0] = vertex_child_idx[0];
	  ccn_tmp[1] = edge_midpoint_idx[0];
	  ccn_tmp[2] = center_fine_idx;
	  ccn_tmp[3] = edge_midpoint_idx[3];
	  child_cell_nodes.push_back(ccn_tmp);

	  ccn_tmp[0] = vertex_child_idx[1];
	  ccn_tmp[1] = edge_midpoint_idx[1];
	  ccn_tmp[2] = center_fine_idx;
	  ccn_tmp[3] = edge_midpoint_idx[0];
	  child_cell_nodes.push_back(ccn_tmp);

	  ccn_tmp[0] = vertex_child_idx[2];
	  ccn_tmp[1] = edge_midpoint_idx[1];
	  ccn_tmp[2] = center_fine_idx;
	  ccn_tmp[3] = edge_midpoint_idx[2];
	  child_cell_nodes.push_back(ccn_tmp);

	  ccn_tmp[0] = vertex_child_idx[3];
	  ccn_tmp[1] = edge_midpoint_idx[2];
	  ccn_tmp[2] = center_fine_idx;
	  ccn_tmp[3] = edge_midpoint_idx[3];
	  child_cell_nodes.push_back(ccn_tmp);

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

	  cen_tmp[0] = edge_midpoint_idx[2];
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

      }
    }
  }

}  // namespace lf::refinement
