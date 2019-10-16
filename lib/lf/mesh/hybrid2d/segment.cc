/**
 * @file
 * @brief Some implementation of the Segment class.
 * @author Raffael Casagrande
 * @date   2018-06-22 03:58:14
 * @copyright MIT License
 */

#include "segment.h"
#include "point.h"

namespace lf::mesh::hybrid2d {
nonstd::span<const Entity* const> Segment::SubEntities(
    unsigned rel_codim) const {
  // An impressive way to do double dereferencing!
  auto l = [&](auto i) -> const mesh::Entity& { return **i; };
  if (rel_codim == 1) {
    return {reinterpret_cast<const Entity* const*>(&nodes_[0]), 2};
  }
  // If the segment itself is requested return reference to itself
  if (rel_codim == 0) {
    return {&this_, 1};
  }
  LF_VERIFY_MSG(false, "Segment: rel_codim out of range");
}
}  // namespace lf::mesh::hybrid2d
