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
    // Initialize child information
    std::vector<ChildInfo> cell_child_info(base_mesh->Size(0));
    // Nothing special has to be done for edges and points
    std::vector<ChildInfo> edge_child_info(base_mesh->Size(1));
    std::vector<ChildInfo> point_child_info(base_mesh->Size(2));
    child_infos_.push_back({std::move(cell_child_info),
                            std::move(edge_child_info),
                            std::move(point_child_info)});
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
    ChildInfo &pt_child_info((child_infos_.back())[2][point_index]);
    // Set the information about a points children except the child pointer
    pt_child_info.num_children_ = 1;
    pt_child_info.ref_pat_ = RefPat::rp_copy;
  }
  // Flag all edges as to be split
  for (const mesh::Entity &edge : finest_mesh.Entities(1)) {
    const lf::base::glb_idx_t edge_index = finest_mesh.Index(edge);
    ChildInfo &ed_child_info(((child_infos_.back())[1])[edge_index]);
    ed_child_info.num_children_ = 2;
    ed_child_info.ref_pat_ = RefPat::rp_split;
  }
  // All cells are to be refined regularly
  for (const mesh::Entity &cell : finest_mesh.Entities(0)) {
    const lf::base::glb_idx_t cell_index = finest_mesh.Index(cell);
    ChildInfo &cell_child_info(((child_infos_.back())[0])[cell_index]);
    cell_child_info.num_children_ = 4;
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
  // Retrieve the finest mesh in the hierarchy = parent mesh
  const mesh::Mesh &parent_mesh(*meshes_.back());
  // First run through the vertices, create child vertices and register
  // them with the mesh factory
  // Store child indices in an auxiliary array
  std::vector<ChildInfo> &pt_child_info((child_infos_.back())[2]);
  std::vector<lf::base::glb_idx_t> pt_child_idx(parent_mesh.Size(2), -1);
  for (const mesh::Entity &node : parent_mesh.Entities(2)) {
    const lf::base::glb_idx_t node_index = parent_mesh.Index(node);
    const lf::geometry::Geometry &pt_geo(*node.Geometry());
    Eigen::MatrixXd pt_coords(pt_geo.Global(Eigen::Matrix<double, 0, 1>()));
    if (pt_child_info[node_index].ref_pat_ != RefPat::rp_nil) {
      pt_child_idx[node_index] = mesh_factory_.AddPoint(pt_coords);
    }
  }

  // Now traverse the edges. Depending on the refinement pattern,
  // either copy them or split them
  // In auxiliary arrays store the indices of the endpoints
  // of the new edges and those of midpoints in the case of split edges
  std::vector<lf::base::glb_idx_t> edge_part00_idx(parent_mesh.Size(1), -1);
  std::vector<lf::base::glb_idx_t> edge_part1_idx(parent_mesh.Size(1), -1);
  std::vector<lf::base::glb_idx_t> edge_mp_idx(parent_mesh.Size(1), -1);
  for (const mesh::Entity &edge : parent_mesh.Entities(1)) {
    lf::base::glb_idx_t edge_index = parent_mesh.Index(edge);
    const ChildInfo &edge_ci((child_infos_.back())[1][edge_index]);
    const RefPat edge_refpat(edge_ci.ref_pat_);
    switch (edge_refpat) {
      case RefPat::rp_copy: {
        // Edge has to be duplicated
        Hybrid2DRefinementPattern rp_copy(edge.RefEl(), edge_refpat);
        std::vector<std::unique_ptr<lf::geometry::Geometry>> ed_copy(
            edge.Geometry()->ChildGeometry(rp_copy, 0));
        LF_VERIFY_MSG(ed_copy.size() == 1, "Copy creates only a single child!");
        // Get indices of endpoints in parent mesh

        // Fetch indices of endpoints in
        break;
      }  // end rp_copy
    }    // end switch refpat
  }      // end edge loop
}

}  // namespace lf::refinement
