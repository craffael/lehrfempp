/**
 * @file
 * @brief Declares the AllCodimMeshDataSet
 * @author Raffael Casagrande
 * @date   2018-06-30 01:48:52
 * @copyright MIT License
 */

#ifndef __df4311acf6554f11919b7b1edfc5b3dd
#define __df4311acf6554f11919b7b1edfc5b3dd

#include "lf/base/lf_assert.h"
#include "mesh_data_set.h"
#include "mesh_interface.h"

namespace lf::mesh {

/**
 * @brief Assigns to every entity(all codims) in a mesh a value of type `T`
 * @tparam T The type of value to store with every entity
 * @sa MeshDataSet
 */
template <class T>
class AllCodimMeshDataSet : public MeshDataSet<T> {
  using size_type = Mesh::size_type;
  using dim_t = base::RefEl::dim_t;

 public:
  AllCodimMeshDataSet(const AllCodimMeshDataSet&) = delete;
  AllCodimMeshDataSet(AllCodimMeshDataSet&&) = delete;
  AllCodimMeshDataSet& operator=(const AllCodimMeshDataSet&) = delete;
  AllCodimMeshDataSet& operator=(AllCodimMeshDataSet&&) = delete;
  ~AllCodimMeshDataSet() = default;

  /**
   * @brief Create a new AllCodimMeshDataSet and [Default
   * initialize](https://en.cppreference.com/w/cpp/language/default_initialization)
   * the data.
   * @param mesh The mesh that contains the entities.
   *
   */
  AllCodimMeshDataSet(std::shared_ptr<Mesh> mesh)
      : MeshDataSet<T>(),
        dim_mesh_(mesh->DimMesh()),
        mesh_(std::move(mesh)),
        data_(dim_mesh_ + 1) {
    for (dim_t codim = 0; codim <= dim_mesh_; ++codim) {
      data_[codim].resize(mesh_->Size(codim));
    }
  }

  /**
   * @brief Create a new AllCodimMeshDataSet and initialize the data of every
   * entity with the given value (`T` must be copyable!)
   * @param mesh The mesh that contains the entities.
   * @param init_value The initial value that should be assigned to every
   * entity.
   */
  template <class = typename std::enable_if<
                std::is_copy_constructible<T>::value>::type>
  AllCodimMeshDataSet(std::shared_ptr<Mesh> mesh, T init_value)
      : MeshDataSet<T>(),
        dim_mesh_(mesh->DimMesh()),
        mesh_(std::move(mesh)),
        data_(dim_mesh_ + 1) {
    for (dim_t codim = 0; codim <= dim_mesh_; ++codim) {
      data_[codim] = std::vector<T>(mesh_->Size(codim), init_value);
    }
  }

  T& data(const Entity& e) override {
    LF_ASSERT_MSG(DefinedOn(e), "MeshDataSet is not defined on this entity.");
    return data_[e.Codim()][mesh_->Index(e)];
  }
  const T& data(const Entity& e) const override {
    LF_ASSERT_MSG(DefinedOn(e), "MeshDataSet is not defined on this entity.");
    return data_[e.Codim()][mesh_->Index(e)];
  }
  bool DefinedOn(const Entity& e) const override { return mesh_->Contains(e); }

 private:
  dim_t dim_mesh_;
  std::shared_ptr<Mesh> mesh_;
  std::vector<std::vector<T>> data_;
};

}  // namespace lf::mesh

#endif  // __df4311acf6554f11919b7b1edfc5b3dd
