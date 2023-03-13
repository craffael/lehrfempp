/**
 * @file
 * @brief Implementation of Triangle class
 * @author Raffael Casagrande
 * @date   2018-06-22 04:01:10
 * @copyright MIT License
 */

#include "triangle.h"

#include "point.h"
#include "segment.h"

namespace lf::mesh::hybrid2d {

Triangle::Triangle(size_type index,
                   std::unique_ptr<geometry::Geometry>&& geometry,
                   const Point* corner0, const Point* corner1,
                   const Point* corner2, const Segment* edge0,
                   const Segment* edge1, const Segment* edge2)
    : index_(index),
      geometry_(std::move(geometry)),
      nodes_({corner0, corner1, corner2}),
      edges_({edge0, edge1, edge2}),
      edge_ori_(),
      this_(this) {
  LF_VERIFY_MSG(corner0 != nullptr, "Invalid pointer to corner 0");
  LF_VERIFY_MSG(corner1 != nullptr, "Invalid pointer to corner 1");
  LF_VERIFY_MSG(corner2 != nullptr, "Invalid pointer to corner 2");
  LF_VERIFY_MSG(edge0 != nullptr, "Invalid pointer to edge 0");
  LF_VERIFY_MSG(edge1 != nullptr, "Invalid pointer to edge 1");
  LF_VERIFY_MSG(edge2 != nullptr, "Invalid pointer to edge 2");
  if (geometry_) {
    LF_VERIFY_MSG(geometry_->DimLocal() == 2,
                  "Geometry must describe a 2D cell");
    LF_VERIFY_MSG(geometry_->RefEl() == base::RefEl::kTria(),
                  "Cell geometry must fit a triangle");
  }

  // Consistency check: make sure that sub-entities of sub-entities
  // agree with sub-sub-entities

  const lf::base::RefEl ref_el_tria(lf::base::RefElType::kTria);
  // Run through all edges
  for (int ed_loc_idx = 0; ed_loc_idx < 3; ed_loc_idx++) {
    // Pointer to current edge
    const Segment* ed_ptr = edges_.at(ed_loc_idx);
    // Check for segement type
    LF_VERIFY_MSG(ed_ptr->RefEl() == lf::base::RefEl::kSegment(),
                  "Wrong type for edge " << ed_loc_idx);
    // Fetch nodes from current edge
    auto ed_nodes = ed_ptr->SubEntities(1);

    // Obtain local indices of endpoints of edge
    const size_type loc_idx_p0 =
        ref_el_tria.SubSubEntity2SubEntity(1, ed_loc_idx, 1, 0);
    const size_type loc_idx_p1 =
        ref_el_tria.SubSubEntity2SubEntity(1, ed_loc_idx, 1, 1);
    const Point* p0_ptr = nodes_.at(loc_idx_p0);
    const Point* p1_ptr = nodes_.at(loc_idx_p1);

    // Verify that the nodes of the triangle are nodes of the edge as well
    LF_VERIFY_MSG((ed_nodes[0] == p0_ptr) || (ed_nodes[0] == p1_ptr),
                  "Node 0 of edge " << ed_loc_idx << " not a triangle node");
    LF_VERIFY_MSG((ed_nodes[1] == p0_ptr) || (ed_nodes[1] == p1_ptr),
                  "Node 1 of edge " << ed_loc_idx << " not a triangle node");
  }

  // Finally set relative orientations for the edges. Edge i has positive
  // orientation, if its first node agrees with vertex i
  for (int ed_loc_idx = 0; ed_loc_idx < 3; ed_loc_idx++) {
    // Fetch nodes of current edge
    auto ed_nodes = edges_[ed_loc_idx]->SubEntities(1);
    edge_ori_[ed_loc_idx] = (ed_nodes[0] == nodes_[ed_loc_idx])
                                ? lf::mesh::Orientation::positive
                                : lf::mesh::Orientation::negative;
  }
}

// Acessing sub-entities
nonstd::span<const Entity* const> Triangle::SubEntities(
    unsigned rel_codim) const {
  auto l = [&](auto i) -> const mesh::Entity& { return **i; };
  switch (rel_codim) {
    case 2:
      return {reinterpret_cast<const Entity* const*>(nodes_.data()), 3};
    case 1:
      return {reinterpret_cast<const Entity* const*>(edges_.data()), 3};
    case 0:
      return {&this_, 1};
    default:
      LF_VERIFY_MSG(false, "Triangle: rel_codim out of range");
  }
}
}  // namespace lf::mesh::hybrid2d
