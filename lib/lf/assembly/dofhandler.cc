/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Rudimentary implementation of a general DOF handler interface
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#include "dofhandler.h"

namespace lf::assemble {
UniformFEDofHandler::UniformFEDofHandler(std::shared_ptr<lf::mesh::Mesh> mesh,
                                         const LocalStaticDOFs &locdof)
    : mesh_(std::move(mesh)) {
  LF_VERIFY_MSG((mesh_->DimMesh() == 2),
                "UniformFEDofHandler implemented for 2D only");

  // First fetch number of local shape functions covering each entity type
  no_dofs_[kNodeOrd] = locdof.TotalNoLocDofs(lf::base::RefEl::kPoint());
  no_dofs_[kEdgeOrd] = locdof.TotalNoLocDofs(lf::base::RefEl::kSegment());
  num_dofs_tria_ = locdof.TotalNoLocDofs(lf::base::RefEl::kTria());
  num_dofs_quad_ = locdof.TotalNoLocDofs(lf::base::RefEl::kQuad());
  no_dofs_[kCellOrd] = std::max(num_dofs_tria_, num_dofs_quad_);

  gdof_idx_t dof_idx = 0;
  // Step I: Set indices for shape functions on nodes
  // Total number of degrees of freedom on nodes (entities of co-dim = 2)
  const size_type num_dofs_nodes = mesh_->Size(2) * no_dofs_[kNodeOrd];
  dofs_[kNodeOrd].resize(num_dofs_nodes);
  // Run through nodes
  for (const lf::mesh::Entity &node : mesh_->Entities(2)) {
    const glb_idx_t node_idx = mesh_->Index(node);
    glb_idx_t node_dof_offset = node_idx * no_dofs_[kNodeOrd];
    for (int j = 0; j < locdof.NoLocDofs(lf::base::RefEl::kPoint()); j++) {
      dofs_[kNodeOrd][node_dof_offset++] = dof_idx;
      dof_entities_.push_back(&node);  // Store entity for current dof
      dof_idx++;                       // Move on to next index
    }
  }

  // Step II: Set indices for shape functions on edges
  // Total number of degrees of freedom belonging to edges (entities of co-dim =
  // 1)
  const size_type num_dofs_edges = mesh_->Size(1) * no_dofs_[kEdgeOrd];
  dofs_[kEdgeOrd].resize(num_dofs_edges);
  // Visit all edges
  for (const lf::mesh::Entity &edge : mesh_->Entities(1)) {
    const glb_idx_t edge_idx = mesh_->Index(edge);
    glb_idx_t edge_dof_offset = edge_idx * no_dofs_[kEdgeOrd];

    // Obtain indices for basis functions sitting at endpoints
    for (const lf::mesh::Entity &endpoint : edge.SubEntities(1)) {
      const glb_idx_t ep_idx(mesh_->Index(endpoint));
      glb_idx_t ep_dof_offset = ep_idx * no_dofs_[kNodeOrd];
      // Copy indices of shape functions from nodes to edge
      for (int j = 0; j < no_dofs_[kNodeOrd]; j++) {
        dofs_[kEdgeOrd][edge_dof_offset++] = dofs_[kNodeOrd][ep_dof_offset++];
      }
    }
    // Set indices for interior edge degrees of freedom
    for (int j = 0; j < locdof.NoLocDofs(lf::base::RefEl::kSegment()); j++) {
      dofs_[kEdgeOrd][edge_dof_offset++] = dof_idx;
      dof_entities_.push_back(&edge);
      dof_idx++;
    }
  }

  // Step III: Set indices for shape functions on cells
  const size_type max_num_dof_cells = mesh_->Size(0) * no_dofs_[kCellOrd];
  dofs_[kCellOrd].resize(max_num_dof_cells);
  
  // Number of (non-)interior shape functins for edges
  const size_type no_int_dof_edge = locdof.NoLocDofs(lf::base::RefEl::kSegment());
  const size_type num_ext_dof_edge = no_dofs_[kEdgeOrd] - no_int_dof_edge;
  
  // Visit all cells
  for (const lf::mesh::Entity &cell : mesh_->Entities(0)) {
    const glb_idx_t cell_idx = mesh_->Index(cell);
    // Offset for cell dof indices in large dof index vector
    glb_idx_t cell_dof_offset = cell_idx * no_dofs_[kCellOrd];
    
    // Obtain indices for basis functions in vertices
    for (const lf::mesh::Entity &vertex : cell.SubEntities(2)) {
      const glb_idx_t vt_idx(mesh_->Index(vertex));
      glb_idx_t vt_dof_offset = vt_idx * no_dofs_[kNodeOrd];
      // Copy indices of shape functions from nodes to cell
      for (int j = 0; j < no_dofs_[kNodeOrd]; j++) {
        dofs_[kCellOrd][cell_dof_offset++] = dofs_[kNodeOrd][vt_dof_offset++];
      }
    }
    
    // Collect indices of interior shape functions of edges
    // Internal ordering will depend on the orientation of the edge
    lf::base::RandomAccessRange<const lf::mesh::Orientation>
      edge_orientations(cell.RelativeOrientations());
    lf::base::RandomAccessRange<const lf::mesh::Entity>
      edges(cell.SubEntities(1));
    // Loop over edges
    for (int ed_sub_idx = 0; ed_sub_idx < cell.RefEl().NumSubEntities(1); ed_sub_idx++) {
      const glb_idx_t edge_idx = mesh_->Index(edges[ed_sub_idx]);
      glb_idx_t edge_int_dof_offset = edge_idx * no_dofs_[kEdgeOrd] + num_ext_dof_edge;
      // Copy indices of shape functions from edges to cell
      // The order, in which they are copied depends on the relative orientation of
      // the edge w.r.t. the cell
      switch (edge_orientations[ed_sub_idx]) {
      case lf::mesh::Orientation::positive: { 
	for (int j = 0; j < no_int_dof_edge; j++) {
	  dofs_[kCellOrd][cell_dof_offset++] = dofs_[kEdgeOrd][edge_int_dof_offset+j];
	}
	break;
      }
      case lf::mesh::Orientation::negative: { 
	for (int j = no_int_dof_edge-1; j>=0; j--) {
	  dofs_[kCellOrd][cell_dof_offset++] = dofs_[kEdgeOrd][edge_int_dof_offset+j];
	}
	break;
      }
      } // end switch
    }
    
    // Set indices for interior edge degrees of freedom depending on the type of cell
    // Here we add new degrees of freedom
    size_type num_int_dofs_cell;
    if (cell.RefEl() == lf::base::RefEl::kTria())
      num_int_dofs_cell = locdof.NoLocDofs(lf::base::RefEl::kTria());
    else if (cell.RefEl() == lf::base::RefEl::kQuad())
      num_int_dofs_cell = locdof.NoLocDofs(lf::base::RefEl::kQuad());
    else LF_ASSERT_MSG(false,"Illegal entity type");
      
    for (int j = 0; j < num_int_dofs_cell; j++) {
      dofs_[kCellOrd][cell_dof_offset++] = dof_idx;
      dof_entities_.push_back(&cell);
      dof_idx++;
    }
  }
  // Finally store total number of shape functions on the mesh.
  num_dof_ = dof_idx;
} // end contructore

lf::base::RandomAccessRange<const gdof_idx_t>
UniformFEDofHandler::GetGlobalDofs(lf::base::RefEl ref_el_type,
                                   glb_idx_t entity_index) const {
  dim_t codim;  // Co-dimension of entity in a 2D mesh
  size_type no_loc_dofs;
  switch (ref_el_type) {
    case lf::base::RefEl::kPoint(): {
      codim = 0;
      no_loc_dofs = no_dofs_[kNodeOrd];
      break;
    }
    case lf::base::RefEl::kSegment(): {
      codim = 1;
      no_loc_dofs = no_dofs_[kEdgeOrd];
      break;
    }
    case lf::base::RefEl::kTria(): {
      codim = 2;
      no_loc_dofs = num_dofs_tria_;
      break;
    }
    case lf::base::RefEl::kQuad(): {
      codim = 2;
      no_loc_dofs = num_dofs_quad_;
      break;
    }
    default: {
      LF_VERIFY_MSG(false, "Illegal entity type");
      break;
    }
  }
  LF_ASSERT_MSG((mesh_->Size(codim) >= entity_index),
                "Index " << entity_index << " out of range");
  // Pointers to range of dof indices
  const gdof_idx_t *begin =
      dofs_[codim].data() + (no_dofs_[codim] * entity_index);
  const gdof_idx_t *end = begin + no_loc_dofs;
  return {begin, end};
}

lf::base::RandomAccessRange<const gdof_idx_t>
UniformFEDofHandler::GetGlobalDofs(const lf::mesh::Entity &entity) const {
  return GetGlobalDofs(entity.RefEl(), mesh_->Index(entity));
}

size_type UniformFEDofHandler::GetNoDofs(const lf::mesh::Entity &entity) const {
  switch (entity.RefEl()) {
    case lf::base::RefEl::kPoint(): {
      return no_dofs_[kNodeOrd];
    }
    case lf::base::RefEl::kSegment(): {
      return no_dofs_[kEdgeOrd];
    }
    case lf::base::RefEl::kTria(): {
      return num_dofs_tria_;
    }
    case lf::base::RefEl::kQuad(): {
      return num_dofs_quad_;
    }
    default: {
      LF_VERIFY_MSG(false, "Illegal entity type");
      break;
    }
  }
  return (size_type)0;
}

}  // namespace lf::assemble
