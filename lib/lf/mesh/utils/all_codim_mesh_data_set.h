/**
 * @file
 * @brief Declares the AllCodimMeshDataSet
 * @author Raffael Casagrande
 * @date   2018-06-30 01:48:52
 * @copyright MIT License
 */

#ifndef __df4311acf6554f11919b7b1edfc5b3dd
#define __df4311acf6554f11919b7b1edfc5b3dd

#include <lf/mesh/mesh.h>

#include <boost/container/vector.hpp>

#include "lf/base/lf_assert.h"
#include "mesh_data_set.h"

namespace lf::mesh::utils {

/**
 * @brief Assigns to every entity(all codims) in a mesh a value of type `T`
 * @tparam T The type of value to store with every entity
 *
 * @sa MeshDataSet
 */
template <class T>
class AllCodimMeshDataSet : public MeshDataSet<T> {
 public:
  using size_type = Mesh::size_type;
  using dim_t = base::RefEl::dim_t;

  AllCodimMeshDataSet(const AllCodimMeshDataSet&) = delete;
  AllCodimMeshDataSet(AllCodimMeshDataSet&&) noexcept = default;
  AllCodimMeshDataSet& operator=(const AllCodimMeshDataSet&) = delete;
  AllCodimMeshDataSet& operator=(AllCodimMeshDataSet&&) noexcept = default;
  ~AllCodimMeshDataSet() override = default;

  /**
   * @brief Get a (modifiable) reference to the data stored with entity e.
   * @param e The entity whose data should be retrieved/modified
   * @return  A reference to the stored data.
   *
   * @note The behavior of this method is undefined if `DefinedOn(e) == false`!
   */
  [[nodiscard]] T& operator()(const Entity& e) {
    LF_ASSERT_MSG(DefinedOn(e), "MeshDataSet is not defined on this entity.");
    return data_[e.Codim()][mesh_->Index(e)];
  }
  // NOLINTNEXTLINE(readability-const-return-type)
  [[nodiscard]] const T& operator()(const Entity& e) const override {
    LF_ASSERT_MSG(DefinedOn(e), "MeshDataSet is not defined on this entity.");
    return data_[e.Codim()][mesh_->Index(e)];
  }
  [[nodiscard]] bool DefinedOn(const Entity& e) const override {
    return mesh_->Contains(e);
  }

  /** @brief set up default-initialized data arrays
   *
   * @param mesh shared pointer to underlying mesh
   *
   * allocates arrays containing a value of type T for _every_ entity of the
   * mesh
   */
  explicit AllCodimMeshDataSet(
      const std::shared_ptr<const lf::mesh::Mesh>& mesh)
      : MeshDataSet<T>(),
        dim_mesh_(mesh->DimMesh()),
        mesh_(mesh),
        data_(dim_mesh_ + 1) {
    for (dim_t codim = 0; codim <= dim_mesh_; ++codim) {
      data_[codim].resize(mesh_->NumEntities(codim));
    }
  }

  /** @brief Set up data array for entities and initialize it with a given value
   *
   * @param mesh shared pointer to underlying mesh
   * @param init_value value to be copied in every cell of the array
   *
   * Allocates an array for data of type T index by _all_ entities of the
   * mesh.
   */
  template <class X,
            class = typename std::enable_if_t<std::is_convertible_v<X, T> &&
                                              std::is_copy_constructible_v<X>>>
  AllCodimMeshDataSet(const std::shared_ptr<const lf::mesh::Mesh>& mesh,
                      X init_value)
      : MeshDataSet<T>(),
        dim_mesh_(mesh->DimMesh()),
        mesh_(mesh),
        data_(dim_mesh_ + 1) {
    for (dim_t codim = 0; codim <= dim_mesh_; ++codim) {
      data_[codim] =
          boost::container::vector<T>(mesh_->NumEntities(codim), init_value);
    }
  }

 private:
  dim_t dim_mesh_;
  std::shared_ptr<const lf::mesh::Mesh> mesh_;

  // Note(craffael): We use boost container here which doesn't have the ridicule
  // specialization that std::vector<bool> has. In this way we can actually
  // return references in operator() const.
  boost::container::vector<boost::container::vector<T>> data_;
};

/**
 * @brief Create a new AllCodimMeshDataSet and [Default
 * initialize](https://en.cppreference.com/w/cpp/language/default_initialization)
 * the data.
 * @param mesh The mesh that contains the entities.
 *
 */
template <class T>
std::shared_ptr<AllCodimMeshDataSet<T>> make_AllCodimMeshDataSet(
    const std::shared_ptr<const lf::mesh::Mesh>& mesh) {
  using impl_t = AllCodimMeshDataSet<T>;
  auto t = new AllCodimMeshDataSet<T>(mesh);
  return std::shared_ptr<impl_t>(t);
}

/**
 * @brief Create a new AllCodimMeshDataSet and initialize the data of every
 * entity with the given value (`T` must be copyable!)
 * @param mesh The mesh that contains the entities.
 * @param init_value The initial value that should be assigned to every
 * entity.
 */
template <class T, class = typename std::enable_if<
                       std::is_copy_constructible<T>::value>::type>
std::shared_ptr<AllCodimMeshDataSet<T>> make_AllCodimMeshDataSet(
    const std::shared_ptr<const lf::mesh::Mesh>& mesh, T init_value) {
  using impl_t = AllCodimMeshDataSet<T>;
  return std::shared_ptr<impl_t>(new impl_t(mesh, init_value));
}

}  // namespace lf::mesh::utils

#endif  // __df4311acf6554f11919b7b1edfc5b3dd
