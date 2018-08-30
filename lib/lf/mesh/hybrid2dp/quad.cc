/**
 * @file
 * @brief Implementations from quad.h
 * @author Raffael Casagrande
 * @date   2018-06-22 04:05:20
 * @copyright MIT License
 */

#include "quad.h"
#include "point.h"
#include "segment.h"

namespace lf::mesh::hybrid2dp {
Quadrilateral::Quadrilateral(size_type index,
                             std::unique_ptr<geometry::Geometry>&& geometry,
                             const Point* corner0, const Point* corner1,
                             const Point* corner2, const Point* corner3,
                             const Segment* edge0, const Segment* edge1,
                             const Segment* edge2, const Segment* edge3)
    : index_(index),
      geometry_(std::move(geometry)),
      nodes_({corner0, corner1, corner2, corner3}),
      edges_({edge0, edge1, edge2, edge3}) {
  LF_VERIFY_MSG(corner0 != nullptr, "Invalid pointer to corner 0");
  LF_VERIFY_MSG(corner1 != nullptr, "Invalid pointer to corner 1");
  LF_VERIFY_MSG(corner2 != nullptr, "Invalid pointer to corner 2");
  LF_VERIFY_MSG(corner3 != nullptr, "Invalid pointer to corner 3");
  LF_VERIFY_MSG(edge0 != nullptr, "Invalid pointer to edge 0");
  LF_VERIFY_MSG(edge1 != nullptr, "Invalid pointer to edge 1");
  LF_VERIFY_MSG(edge2 != nullptr, "Invalid pointer to edge 2");
  LF_VERIFY_MSG(edge3 != nullptr, "Invalid pointer to edge 3");
  if (geometry_) {
    LF_VERIFY_MSG(geometry_->DimLocal() == 2,
                  "Geometry must describe a 2D cell");
    LF_VERIFY_MSG(geometry_->RefEl() == base::RefEl::kQuad(),
                  "Cell geometry must fit a quad");
  }

  // Consistency check: make sure that sub-entities of sub-entities
  // agree with sub-sub-entities

  const lf::base::RefEl ref_el_quad(lf::base::RefElType::kQuad);
  // Run through all edges
  for (int ed_loc_idx = 0; ed_loc_idx < 4; ed_loc_idx++) {
    // Pointer to current edge
    const Segment* ed_ptr = edges_.at(ed_loc_idx);
    // Check for segement type
    LF_VERIFY_MSG(ed_ptr->RefEl() == lf::base::RefEl::kSegment(),
                  "Wrong type for edge " << ed_loc_idx);
    // Fetch nodes from current edge
    base::RandomAccessRange<const Entity> ed_nodes(ed_ptr->SubEntities(1));

    // Obtain local indices of endpoints of edge
    const size_type loc_idx_p0 =
        ref_el_quad.SubSubEntity2SubEntity(1, ed_loc_idx, 1, 0);
    const size_type loc_idx_p1 =
        ref_el_quad.SubSubEntity2SubEntity(1, ed_loc_idx, 1, 1);
    const Point* p0_ptr = nodes_.at(loc_idx_p0);
    const Point* p1_ptr = nodes_.at(loc_idx_p1);

    // Verify that the nodes of the quadrilateral are nodes of the edge as well
    LF_VERIFY_MSG((ed_nodes[0] == *p0_ptr) || (ed_nodes[0] == *p1_ptr),
                  "Node 0 of edge " << ed_loc_idx << " not a quad node");
    LF_VERIFY_MSG((ed_nodes[1] == *p0_ptr) || (ed_nodes[1] == *p1_ptr),
                  "Node 1 of edge " << ed_loc_idx << " not a quad node");
  }

  // Finally set relative orientations for the edges. Edge i has positive
  // orientation, if its first node agrees with vertex i
  for (int ed_loc_idx = 0; ed_loc_idx < 4; ed_loc_idx++) {
    // Fetch nodes of current edge
    base::RandomAccessRange<const Entity> ed_nodes(
        edges_[ed_loc_idx]->SubEntities(1));
    edge_ori_[ed_loc_idx] = (ed_nodes[0] == *nodes_[ed_loc_idx])
                                ? lf::mesh::Orientation::positive
                                : lf::mesh::Orientation::negative;
  }
  
}  // end constructor

// Access to sub-entities
base::RandomAccessRange<const mesh::Entity> Quadrilateral::SubEntities(
    char rel_codim) const {
  auto l = [&](auto i) -> const mesh::Entity& { return **i; };
  switch (rel_codim) {
    case 2:
      return {
          base::make_DereferenceLambdaRandomAccessIterator(nodes_.begin(), l),
          base::make_DereferenceLambdaRandomAccessIterator(nodes_.end(), l)};
    case 1:
      return {
          base::make_DereferenceLambdaRandomAccessIterator(edges_.begin(), l),
          base::make_DereferenceLambdaRandomAccessIterator(edges_.end(), l)};
    case 0:
      return {this, this + 1};
    default:
      LF_VERIFY_MSG(
          false, "Quadrilateral: rel_codim " << rel_codim << " out of range");
  }
}
}  // namespace lf::mesh::hybrid2dp
