#ifndef __7a3f1903d42141a3b1135e8e5ad72c1c
#define __7a3f1903d42141a3b1135e8e5ad72c1c

#include <lf/base/base.h>
#include "lf/mesh/entity.h"
#include "lf/mesh/mesh_interface.h"

namespace lf::mesh::hybrid2d {

class Mesh;

template <char CODIM>
class Entity : public mesh::Entity {
  using size_type = mesh::Mesh::size_type;

 public:
  // needed by std::vector
  Entity() = default;
  Entity(Entity&&) = default;
  Entity& operator=(Entity&&) = default;

  char Codim() const override { return CODIM; }

  base::RandomAccessRange<const mesh::Entity> SubEntities(
    char codim) const override;


  mesh::Geometry* Geometry() const override { return geometry_.get(); }
  base::RefEl RefEl() const override { return geometry_->RefEl(); }

  bool operator==(const mesh::Entity& rhs) const override {
    return this == &rhs;
  }

  // constructor, is called from Mesh
  Entity(Mesh* mesh, size_type index,
         std::unique_ptr<mesh::Geometry>&& geometry,
         std::array<std::vector<size_type>, 2 - CODIM> sub_entities)
      : mesh_(mesh),
        index_(index),
        geometry_(std::move(geometry)),
        sub_entities_(std::move(sub_entities)) {}

 private:
  Mesh* mesh_;
  size_type index_;  // zero-based index of this entity.
  std::unique_ptr<mesh::Geometry> geometry_;
  std::array<std::vector<size_type>, 2 - CODIM> sub_entities_;

  friend class Mesh;
};

}

#include "mesh.h"

namespace lf::mesh::hybrid2d {

template <char CODIM>
base::RandomAccessRange<const mesh::Entity> Entity<CODIM>::SubEntities(
  char codim) const {
  switch (2-CODIM - codim) {
    case 2:
      // return ourselves as the only element:
      return {this, this + 1};
    case 1:
      return {
        base::DereferenceLambdaRandomAccessIterator(
                                                    sub_entities_[codim - 1].
                                                    begin(),
                                                    [&](auto i) -> const mesh::Entity& {
                                                      return mesh_->entities1_[*i];
                                                    }),
        base::DereferenceLambdaRandomAccessIterator(
                                                    sub_entities_[codim - 1].
                                                    begin(),
                                                    [&](auto i) -> const mesh::Entity& {
                                                      return mesh_->entities1_[*i];
                                                    })
      };
    case 0:
      return {
        base::DereferenceLambdaRandomAccessIterator(
                                                    sub_entities_[codim - 1].
                                                    begin(),
                                                    [&](auto i) -> const mesh::Entity& {
                                                      return mesh_->entities2_[*i];
                                                    }),
        base::DereferenceLambdaRandomAccessIterator(
                                                    sub_entities_[codim - 1].
                                                    begin(),
                                                    [&](auto i) -> const mesh::Entity& {
                                                      return mesh_->entities2_[*i];
                                                    })
      };
    default:
    LF_VERIFY_MSG(false, "codim is out of bounds.");
  }
}

}  // namespace lf::mesh::hybrid2d

#endif  // __7a3f1903d42141a3b1135e8e5ad72c1c
