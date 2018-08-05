#ifndef __bae24f2390174bff85f65c2c2e558a9e
#define __bae24f2390174bff85f65c2c2e558a9e

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
   * @brief Create a new CodimMeshDataSet that attaches data of type `T` with
   *        every entity with codimension `codim`. The data is [Default
   * initialized](https://en.cppreference.com/w/cpp/language/default_initialization)
   * @param mesh The mesh that contains the entities.
   * @param codim The codimension of the entities whith which the data is
   * stored.
   *
   */
  CodimMeshDataSet(std::shared_ptr<Mesh> mesh, dim_t codim)
      : MeshDataSet<T>(),
        mesh_(std::move(mesh)),
        data_(mesh_->Size(codim)),
        codim_(codim) {}

  /**
   * @brief Create a new CodimMeshDataSet that attached data of type `T` with
   *        every entity with codimension `codim`. The data of every entity
   *        is initialized to the given value (`T` must be copyable!)
   * @param mesh The mesh that contains the entities.
   * @param codim The codimension of the entities with which the data is stored.
   * @param init The initial value that should be assigned to every entity.
   */
  template <class = typename std::enable_if<
                std::is_copy_constructible<T>::value>::type>
  CodimMeshDataSet(std::shared_ptr<Mesh> mesh, dim_t codim, T init)
      : MeshDataSet<T>(),
        mesh_(std::move(mesh)),
        data_(mesh_->Size(codim), init),
        codim_(codim) {}

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
  std::shared_ptr<Mesh> mesh_;
  std::vector<T> data_;
  dim_t codim_;
};

}  // namespace lf::mesh::utils

#endif  // __bae24f2390174bff85f65c2c2e558a9e
