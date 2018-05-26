#ifndef __7a3f1903d42141a3b1135e8e5ad72c1c
#define __7a3f1903d42141a3b1135e8e5ad72c1c
#include "lf/mesh/entity.h"
#include "lf/mesh/mesh_interface.h"

namespace lf::mesh::hybrid2d {

class Mesh;

template <char CODIM>
class Entity : public mesh::Entity {
  using size_t = mesh::Mesh::size_t;
public:

  // needed by std::vector
  Entity() = default;
  Entity(Entity&&) = default;
  Entity& operator=(Entity&& ) = default;

  char Codim() const override { return CODIM; }

  base::RandomAccessRange<const mesh::Entity> SubEntities(char codim) const override {
    if(codim == 0) {
      // return ourselves as the only element:
      base::RandomAccessIterator<const mesh::Entity> start(this);
      return {start, start+1};
    }

    return { base::DereferenceRandomAccessIterator(sub_entities_[codim-1].begin()),
      base::DereferenceRandomAccessIterator(sub_entities_[codim-1].end())};
  }

  mesh::Geometry* Geometry() const override { return geometry_.get(); }
  base::RefEl RefEl() const override { return geometry_->RefEl(); }

  // constructor, is called from Mesh
  Entity(size_t index, std::unique_ptr<mesh::Geometry>&& geometry,
         std::array<std::vector<mesh::Entity*>, 2 - CODIM> sub_entities)
    : index_(index),
      geometry_(std::move(geometry)),
      sub_entities_(std::move(sub_entities)) {
  }

private:
  size_t index_; // zero-based index of this entity.
  std::unique_ptr<mesh::Geometry> geometry_;
  std::array<std::vector<mesh::Entity*>, 2 - CODIM> sub_entities_;

  friend class Mesh;

};

}

#endif // __7a3f1903d42141a3b1135e8e5ad72c1c
