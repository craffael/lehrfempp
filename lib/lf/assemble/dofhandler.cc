/***************************************************************************
 * LehrFEM++ - A simple C++ finite element libray for teaching
 * Developed from 2018 at the Seminar of Applied Mathematics of ETH Zurich,
 * lead developers Dr. R. Casagrande and Prof. R. Hiptmair
 ***************************************************************************/

/**
 * @file
 * @brief Implementation of a general DOF handler interface
 * @author Ralf Hiptmair
 * @date August 2018
 * @copyright MIT License
 */

#include "dofhandler.h"

namespace lf::assemble {

// Default output flag
unsigned int DofHandler::output_ctrl_ = 0;

// Implementation of output operator for interface class
std::ostream &operator<<(std::ostream &o, const DofHandler &dof_handler) {
  // Obtain pointer to the underlying mesh
  auto mesh = dof_handler.Mesh();
  // Number of degrees of freedom managed by the DofHandler object
  const lf::assemble::size_type N_dofs(dof_handler.NumDofs());

  o << "DofHandler(" << dof_handler.NumDofs() << " dofs)";
  if (DofHandler::output_ctrl_ > 0) {
    // More detailed output
    o << std::endl;
    if (DofHandler::output_ctrl_ % 2 == 0) {
      for (lf::base::dim_t codim = 0; codim <= mesh->DimMesh(); codim++) {
        // Visit all entities of a specific codimension
        for (const lf::mesh::Entity *e : mesh->Entities(codim)) {
          // Fetch unique index of current entity supplied by mesh object
          const lf::base::glb_idx_t e_idx = mesh->Index(*e);
          // Number of shape functions covering current entity
          const lf::assemble::size_type no_dofs(dof_handler.NumLocalDofs(*e));
          // Obtain global indices of those shape functions ...
          nonstd::span<const gdof_idx_t> doflist(
              dof_handler.GlobalDofIndices(*e));
          // and print them
          o << *e << ' ' << e_idx << ": " << no_dofs << " dofs = [";
          for (const lf::assemble::gdof_idx_t &dof : doflist) {
            o << dof << ' ';
          }
          o << ']';
          if (DofHandler::output_ctrl_ % 5 == 0) {
            // Also output indices of interior shape functions
            nonstd::span<const gdof_idx_t> intdoflist(
                dof_handler.InteriorGlobalDofIndices(*e));
            o << " int = [";
            for (lf::assemble::gdof_idx_t int_dof : intdoflist) {
              o << int_dof << ' ';
            }
            o << ']';
          }
          o << std::endl;
        }
      }
    }
    if (DofHandler::output_ctrl_ % 3 == 0) {
      // List entities associated wit the dofs managed by the current
      // DofHandler object
      for (lf::assemble::gdof_idx_t dof_idx = 0; dof_idx < N_dofs; dof_idx++) {
        const lf::mesh::Entity &e(dof_handler.Entity(dof_idx));
        o << "dof " << dof_idx << " -> " << e << " " << mesh->Index(e)
          << std::endl;
      }
    }
  }
  return o;
}

// ----------------------------------------------------------------------
// Implementation UniformFEDofHandler
// ----------------------------------------------------------------------

UniformFEDofHandler::UniformFEDofHandler(
    std::shared_ptr<const lf::mesh::Mesh> mesh, dof_map_t dofmap)
    : mesh_(std::move(mesh)), num_dofs_() {
  LF_ASSERT_MSG((mesh_->DimMesh() == 2), "Can handle 2D meshes only");

  // For checking whether a key was found
  auto map_end = dofmap.end();

  // Get no of interior dof specified for nodes
  auto map_el_pt = dofmap.find(lf::base::RefEl::kPoint());
  if (map_el_pt != map_end) {
    num_loc_dof_point_ = map_el_pt->second;
  }

  // Get no of interior dof specified for edges
  auto map_el_ed = dofmap.find(lf::base::RefEl::kSegment());
  if (map_el_ed != map_end) {
    num_loc_dof_segment_ = map_el_ed->second;
  }

  // Get no of interior dof specified for triangles
  auto map_el_tr = dofmap.find(lf::base::RefEl::kTria());
  if (map_el_tr != map_end) {
    num_loc_dof_tria_ = map_el_tr->second;
  }

  // Get no of interior dof specified for quads
  auto map_el_qd = dofmap.find(lf::base::RefEl::kQuad());
  if (map_el_qd != map_end) {
    num_loc_dof_quad_ = map_el_qd->second;
  }

  // If an entity type is not represented in the map, we assume that
  // no shape functions are attached to those entities

  // Initialize total number of shape functions covering an entity.
  InitTotalNumDofs();

  // Initializatin of dof index arrays
  initIndexArrays();
}

void UniformFEDofHandler::InitTotalNumDofs() {
  num_dofs_[kNodeOrd] = num_loc_dof_point_;
  num_dofs_[kEdgeOrd] = 2 * num_loc_dof_point_ + num_loc_dof_segment_;
  num_dofs_tria_ =
      3 * num_loc_dof_point_ + 3 * num_loc_dof_segment_ + num_loc_dof_tria_;
  num_dofs_quad_ =
      4 * num_loc_dof_point_ + 4 * num_loc_dof_segment_ + num_loc_dof_quad_;
  num_dofs_[kCellOrd] = std::max(num_dofs_tria_, num_dofs_quad_);
}

void UniformFEDofHandler::initIndexArrays() {
  // This method assumes a proper initialization
  // of the data in no_loc_dof_* and num_dofs_, num_dof_tria, num_dofs_quad_
  gdof_idx_t dof_idx = 0;

  // Step I: Set indices for shape functions on nodes
  // Total number of degrees of freedom on nodes (entities of co-dim = 2)
  const size_type no_nodes = mesh_->NumEntities(2);
  const size_type num_dofs_nodes = no_nodes * num_dofs_[kNodeOrd];
  dofs_[kNodeOrd].resize(num_dofs_nodes);
  // Run through nodes in the order given by their numbering
  for (glb_idx_t node_idx = 0; node_idx < no_nodes; node_idx++) {
    const mesh::Entity *node_p{mesh_->EntityByIndex(2, node_idx)};
    LF_ASSERT_MSG(mesh_->Index(*node_p) == node_idx, "Node index mismatch");
    // Beginning of section for concrete node in the dof index vector
    // for entities of co-dimension 2
    glb_idx_t node_dof_offset = node_idx * num_dofs_[kNodeOrd];
    for (unsigned j = 0; j < num_loc_dof_point_; j++) {
      dofs_[kNodeOrd][node_dof_offset++] = dof_idx;
      dof_entities_.push_back(node_p);  // Store entity for current dof
      dof_idx++;                        // Move on to next index
    }
  }

  // Step II: Set indices for shape functions on edges
  // Total number of degrees of freedom belonging to edges (entities of co-dim =
  // 1)
  const size_type no_edges = mesh_->NumEntities(1);
  const size_type num_dofs_edges = no_edges * num_dofs_[kEdgeOrd];
  dofs_[kEdgeOrd].resize(num_dofs_edges);
  // Visit all edges
  // Old implementation, see remarks above
  // for (const lf::mesh::Entity &edge : mesh_->Entities(1)) {
  for (glb_idx_t edge_idx = 0; edge_idx < no_edges; edge_idx++) {
    // Obtain pointer to edge entity
    const mesh::Entity *edge_p{mesh_->EntityByIndex(1, edge_idx)};
    LF_ASSERT_MSG(mesh_->Index(*edge_p) == edge_idx, "Edge index mismatch");
    // Beginning of section for concrete edge in the dof index vector
    // for entities of co-dimension 1
    glb_idx_t edge_dof_offset = edge_idx * num_dofs_[kEdgeOrd];

    // Obtain indices for basis functions sitting at endpoints
    for (const lf::mesh::Entity *endpoint : edge_p->SubEntities(1)) {
      const glb_idx_t ep_idx(mesh_->Index(*endpoint));
      glb_idx_t ep_dof_offset = ep_idx * num_dofs_[kNodeOrd];
      // Copy indices of shape functions from nodes to edge
      for (unsigned j = 0; j < num_dofs_[kNodeOrd]; j++) {
        dofs_[kEdgeOrd][edge_dof_offset++] = dofs_[kNodeOrd][ep_dof_offset++];
      }
    }
    // Set indices for interior edge degrees of freedom
    for (unsigned j = 0; j < num_loc_dof_segment_; j++) {
      dofs_[kEdgeOrd][edge_dof_offset++] = dof_idx;
      dof_entities_.push_back(edge_p);
      dof_idx++;
    }
  }

  // Step III: Set indices for shape functions on cells
  const size_type no_cells = mesh_->NumEntities(0);
  const size_type max_num_dof_cells = no_cells * num_dofs_[kCellOrd];
  dofs_[kCellOrd].resize(max_num_dof_cells);

  // Number of (non-)interior shape functins for edges
  const size_type no_int_dof_edge = num_loc_dof_segment_;
  const size_type num_ext_dof_edge = num_dofs_[kEdgeOrd] - no_int_dof_edge;

  // Visit all cells
  // Old implementation without strong link between cell
  // indices and ordering of global shape functions
  // for (const lf::mesh::Entity &cell : mesh_->Entities(0)) {
  for (glb_idx_t cell_idx = 0; cell_idx < no_cells; cell_idx++) {
    // Obtain pointer to current ell
    const mesh::Entity *cell_p{mesh_->EntityByIndex(0, cell_idx)};
    LF_ASSERT_MSG(cell_idx == mesh_->Index(*cell_p), "cell index mismatch");
    // Offset for cell dof indices in large dof index vector
    glb_idx_t cell_dof_offset = cell_idx * num_dofs_[kCellOrd];

    // Obtain indices for basis functions in vertices
    for (const lf::mesh::Entity *vertex : cell_p->SubEntities(2)) {
      const glb_idx_t vt_idx(mesh_->Index(*vertex));
      glb_idx_t vt_dof_offset = vt_idx * num_dofs_[kNodeOrd];
      // Copy indices of shape functions from nodes to cell
      for (unsigned j = 0; j < num_dofs_[kNodeOrd]; j++) {
        dofs_[kCellOrd][cell_dof_offset++] = dofs_[kNodeOrd][vt_dof_offset++];
      }
    }

    // Collect indices of interior shape functions of edges
    // Internal ordering will depend on the orientation of the edge
    auto edge_orientations = cell_p->RelativeOrientations();
    auto edges = cell_p->SubEntities(1);
    // Loop over edges
    const size_type no_edges_cell = cell_p->RefEl().NumSubEntities(1);
    for (int ed_sub_idx = 0; ed_sub_idx < no_edges_cell; ed_sub_idx++) {
      const glb_idx_t edge_idx = mesh_->Index(*edges[ed_sub_idx]);
      glb_idx_t edge_int_dof_offset =
          edge_idx * num_dofs_[kEdgeOrd] + num_ext_dof_edge;
      // Copy indices of shape functions from edges to cell
      // The order, in which they are copied depends on the relative orientation
      // of the edge w.r.t. the cell
      switch (edge_orientations[ed_sub_idx]) {
        case lf::mesh::Orientation::positive: {
          for (int j = 0; j < no_int_dof_edge; j++) {
            dofs_[kCellOrd][cell_dof_offset++] =
                dofs_[kEdgeOrd][edge_int_dof_offset + j];
          }
          break;
        }
        case lf::mesh::Orientation::negative: {
          for (int j = static_cast<int>(no_int_dof_edge - 1); j >= 0; j--) {
            dofs_[kCellOrd][cell_dof_offset++] =
                dofs_[kEdgeOrd][edge_int_dof_offset + j];
          }
          break;
        }
      }  // end switch
    }

    // Set indices for interior cell degrees of freedom depending on the type of
    // cell. Here we add new degrees of freedom
    size_type num_int_dofs_cell;
    if (cell_p->RefEl() == lf::base::RefEl::kTria()) {
      num_int_dofs_cell = num_loc_dof_tria_;
    } else if (cell_p->RefEl() == lf::base::RefEl::kQuad()) {
      num_int_dofs_cell = num_loc_dof_quad_;
    } else {
      LF_ASSERT_MSG(
          false, "Illegal cell type; only triangles and quads are supported");
    }

    // enlist new interior cell-associated dofs
    for (unsigned j = 0; j < num_int_dofs_cell; j++) {
      dofs_[kCellOrd][cell_dof_offset++] = dof_idx;
      dof_entities_.push_back(cell_p);
      dof_idx++;
    }
  }
  // Finally store total number of shape functions on the mesh.
  num_dof_ = dof_idx;
}  // end constructor

nonstd::span<const gdof_idx_t> UniformFEDofHandler::GlobalDofIndices(
    lf::base::RefEl ref_el_type, glb_idx_t entity_index) const {
  // Co-dimension of entity in a 2D mesh
  const dim_t codim = 2 - ref_el_type.Dimension();
  const size_type no_covered_dofs = NumCoveredDofs(ref_el_type);

  LF_ASSERT_MSG((mesh_->NumEntities(codim) > entity_index),
                "Index " << entity_index << " out of range");
  // Pointers to range of dof indices
  const gdof_idx_t *begin =
      dofs_[codim].data() + (num_dofs_[codim] * entity_index);
  const gdof_idx_t *end = begin + no_covered_dofs;
  return {begin, end};
}

nonstd::span<const gdof_idx_t> UniformFEDofHandler::GlobalDofIndices(
    const lf::mesh::Entity &entity) const {
  return GlobalDofIndices(entity.RefEl(), mesh_->Index(entity));
}

nonstd::span<const gdof_idx_t> UniformFEDofHandler::InteriorGlobalDofIndices(
    lf::base::RefEl ref_el_type, glb_idx_t entity_index) const {
  // Co-dimension of entity in a 2D mesh
  const dim_t codim = 2 - ref_el_type.Dimension();
  const size_type no_covered_dofs = NumCoveredDofs(ref_el_type);
  const size_type no_loc_dofs = NumInteriorDofs(ref_el_type);

  LF_ASSERT_MSG((mesh_->NumEntities(codim) > entity_index),
                "Index " << entity_index << " out of range");
  // Pointers to range of dof indices
  const gdof_idx_t *begin =
      dofs_[codim].data() + (num_dofs_[codim] * entity_index);
  const gdof_idx_t *end = begin + no_covered_dofs;
  begin += (no_covered_dofs - no_loc_dofs);
  return {begin, end};
}

nonstd::span<const gdof_idx_t> UniformFEDofHandler::InteriorGlobalDofIndices(
    const lf::mesh::Entity &entity) const {
  return InteriorGlobalDofIndices(entity.RefEl(), mesh_->Index(entity));
}

size_type UniformFEDofHandler::NumLocalDofs(
    const lf::mesh::Entity &entity) const {
  return GetNumLocalDofs(entity.RefEl(), 0);
}

size_type UniformFEDofHandler::NumInteriorDofs(
    const lf::mesh::Entity &entity) const {
  return NumInteriorDofs(entity.RefEl());
}

// ----------------------------------------------------------------------
// Implementation DynamicFEDofHandler
// ----------------------------------------------------------------------

nonstd::span<const gdof_idx_t> DynamicFEDofHandler::GlobalDofIndices(
    const lf::mesh::Entity &entity) const {
  // Topological type
  const lf::base::RefEl ref_el_type = entity.RefEl();
  // Co-dimension of entity in a 2D mesh
  const dim_t codim = 2 - ref_el_type.Dimension();
  // Index of current mesh entity
  const glb_idx_t entity_index = mesh_p_->Index(entity);
  // Offset of dof indices for current entity
  const size_type idx_offset = offsets_[codim][entity_index];
  // Number of shape functions covering current entity
  const size_type no_covered_dofs =
      offsets_[codim][entity_index + 1] - idx_offset;
  // Pointers to range of dof indices
  const gdof_idx_t *begin = dofs_[codim].data() + idx_offset;
  const gdof_idx_t *end = begin + no_covered_dofs;
  return {begin, end};
}

nonstd::span<const gdof_idx_t> DynamicFEDofHandler::InteriorGlobalDofIndices(
    const lf::mesh::Entity &entity) const {
  // Topological type
  const lf::base::RefEl ref_el_type = entity.RefEl();
  // Co-dimension of entity in a 2D mesh
  const dim_t codim = 2 - ref_el_type.Dimension();
  // Index of current mesh entity
  const glb_idx_t entity_index = mesh_p_->Index(entity);
  // Offset of indices for current entity
  const size_type idx_offset = offsets_[codim][entity_index];
  // Number of shape functions covering current entity
  const size_type no_covered_dofs =
      offsets_[codim][entity_index + 1] - idx_offset;
  // Number of shape functions associated with the current entity
  const size_type no_loc_dofs = num_int_dofs_[codim][entity_index];
  // Pointers to range of dof indices
  const gdof_idx_t *begin = dofs_[codim].data() + idx_offset;
  const gdof_idx_t *end = begin + no_covered_dofs;
  // Move pointer to location of indices for interior dofs
  // This is the only difference to GlobalDofIndices()
  begin += (no_covered_dofs - no_loc_dofs);
  return {begin, end};
}

size_type DynamicFEDofHandler::NumLocalDofs(
    const lf::mesh::Entity &entity) const {
  const lf::base::RefEl ref_el_type = entity.RefEl();
  const dim_t codim = 2 - ref_el_type.Dimension();
  const glb_idx_t entity_index = mesh_p_->Index(entity);
  return (offsets_[codim][entity_index + 1] - offsets_[codim][entity_index]);
}

size_type DynamicFEDofHandler::NumInteriorDofs(
    const lf::mesh::Entity &entity) const {
  const lf::base::RefEl ref_el_type = entity.RefEl();
  const dim_t codim = 2 - ref_el_type.Dimension();
  const glb_idx_t entity_index = mesh_p_->Index(entity);
  const size_type no_loc_dofs = num_int_dofs_[codim][entity_index];
  return no_loc_dofs;
}

}  // namespace lf::assemble
