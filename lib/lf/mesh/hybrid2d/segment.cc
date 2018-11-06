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
base::RandomAccessRange<const mesh::Entity> Segment::SubEntities(
    char rel_codim) const {
  // An impressive way to do double dereferencing!
  auto l = [&](auto i) -> const mesh::Entity& { return **i; };
  if (rel_codim == 1) {
    return {base::make_DereferenceLambdaRandomAccessIterator(nodes_.begin(), l),
            base::make_DereferenceLambdaRandomAccessIterator(nodes_.end(), l)};
  }
  // If the segment itself is requested return reference to itself
  if (rel_codim == 0) {
    return {this, this + 1};
  }
  LF_VERIFY_MSG(false, "Segment: rel_codim out of range");
}
}  // namespace lf::mesh::hybrid2d
