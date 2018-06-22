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

namespace lf::mesh::hybrid2dp {
base::RandomAccessRange<const mesh::Entity> Triangle::SubEntities(
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
      LF_VERIFY_MSG(false, "Triangle: rel_codim out of range");
  }
}
}  // namespace lf::mesh::hybrid2dp
