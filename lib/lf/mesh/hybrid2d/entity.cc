#include "entity.h"

namespace lf::mesh::hybrid2d {
// Implementation of sub-entity access methods

base::RandomAccessRange<const mesh::Entity> Edge::SubEntities(
    char rel_codim) const {
  auto l = [&](auto i) -> const mesh::Entity& { return **i; };
  if (rel_codim == 1) {
    return {base::make_DereferenceLambdaRandomAccessIterator(nodes_.begin(), l),
            base::make_DereferenceLambdaRandomAccessIterator(nodes_.end(), l)};
  }
  if (rel_codim == 1) {
    return {this, this + 1};
  }
  LF_VERIFY_MSG(false, "Edge: rel_codim out of range");
}

base::RandomAccessRange<const mesh::Entity> Trilateral::SubEntities(
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
      LF_VERIFY_MSG(false, "Trilateral: rel_codim out of range");
  }
}

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
      LF_VERIFY_MSG(false, "Trilateral: rel_codim out of range");
  }
}

}  // namespace lf::mesh::hybrid2d
