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
      LF_VERIFY_MSG(false, "Quadrilateral: rel_codim out of range");
  }
}
}  // namespace lf::mesh::hybrid2dp
