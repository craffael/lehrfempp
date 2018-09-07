#ifndef __bae24f2390174bff85f65c2c2e558a9e
#define __bae24f2390174bff85f65c2c2e558a9e

#include <boost/container/vector.hpp>
#include "mesh_data_set.h"

namespace lf::mesh::utils {

/**
 * @brief A MeshDataSet that attaches data of type `T` to every entity of a mesh
 * that has a specified codimension.
 * @tparam T The type of data that should be stored with the entities.
 * @sa MeshDataSet
 */
template <class T>
class CodimMeshDataSet : public MeshDataSet<T> {
 public:
  using dim_t = base::RefEl::dim_t;

  /**
   * @brief Get a (modifiable) reference to the data stored with entity e.
   * @param e The entity whose data should be retrieved/modified
   * @return  A reference to the stored data.
   *
   * @note The behavior of this method is undefined if `DefinedOn(e) == false`!
   */
  T& operator()(const Entity& e) {
    LF_ASSERT_MSG(DefinedOn(e), "MeshDataSet not defined on this entity.");
    return data_[mesh_->Index(e)];
  }
  const T operator()(const Entity& e) const override {
    LF_ASSERT_MSG(DefinedOn(e), "MeshDataSet not defined on this entity.");
    return data_[mesh_->Index(e)];
  }
  bool DefinedOn(const Entity& e) const override {
    return e.Codim() == codim_ && mesh_->Contains(e);
  }

 private:
  // template magic to not use std::vector<bool> but boost::dynamic_bitset
  // instead
  template <class X, class Y = int>
  struct Container {
    // static_assert(false, "gugus");
    using type_t = std::vector<X>;
  };

  template <class Y>
  struct Container<bool, Y> {
    using type_t = boost::container::vector<bool>;
  };

  /** pointer to underlying mesh. Owned by the data set ! */
  std::shared_ptr<const Mesh> mesh_;
  /** data vector */
  typename Container<T, int>::type_t data_;
  /** co-dimension */
  dim_t codim_;

  CodimMeshDataSet(std::shared_ptr<const Mesh> mesh, dim_t codim)
      : MeshDataSet<T>(),
        mesh_(std::move(mesh)),
        data_(mesh_->Size(codim)),
        codim_(codim) {}

  template <class = typename std::enable_if<
                std::is_copy_constructible<T>::value>::type>
  CodimMeshDataSet(std::shared_ptr<const Mesh> mesh, dim_t codim, T init)
      : MeshDataSet<T>(),
        mesh_(std::move(mesh)),
        data_(mesh_->Size(codim), init),
        codim_(codim) {}

  // Friends:
  template <class S>
  friend std::shared_ptr<CodimMeshDataSet<S>> make_CodimMeshDataSet(
      std::shared_ptr<const Mesh> mesh, base::dim_t codim);

  template <class S, class>
  friend std::shared_ptr<CodimMeshDataSet<S>> make_CodimMeshDataSet(
      std::shared_ptr<const Mesh> mesh, base::dim_t codim, S init);
};

/**
 * @brief Create a new CodimMeshDataSet that attaches data of type `T` with
 *        every entity with codimension `codim`. The data is [Default
 * initialized](https://en.cppreference.com/w/cpp/language/default_initialization)
 * @param mesh The mesh that contains the entities.
 * @param codim The codimension of the entities whith which the data is
 * stored.
 *
 */
template <class T>
std::shared_ptr<CodimMeshDataSet<T>> make_CodimMeshDataSet(
    std::shared_ptr<const Mesh> mesh, base::dim_t codim) {
  using impl_t = CodimMeshDataSet<T>;
  return std::shared_ptr<impl_t>(new impl_t(std::move(mesh), codim));
}

/**
 * @brief Create a new CodimMeshDataSet that attached data of type `T` with
 *        every entity with codimension `codim`. The data of every entity
 *        is initialized to the given value (`T` must be copyable!)
 * @param mesh The mesh that contains the entities.
 * @param codim The codimension of the entities with which the data is stored.
 * @param init The initial value that should be assigned to every entity.
 */
template <class T, class = typename std::enable_if<
                       std::is_copy_constructible<T>::value>::type>
std::shared_ptr<CodimMeshDataSet<T>> make_CodimMeshDataSet(
    std::shared_ptr<const Mesh> mesh, base::dim_t codim, T init) {
  using impl_t = CodimMeshDataSet<T>;
  return std::shared_ptr<impl_t>(new impl_t(std::move(mesh), codim, init));
}

}  // namespace lf::mesh::utils

#endif  // __bae24f2390174bff85f65c2c2e558a9e
