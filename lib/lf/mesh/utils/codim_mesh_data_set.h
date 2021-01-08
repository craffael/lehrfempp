#ifndef __bae24f2390174bff85f65c2c2e558a9e
#define __bae24f2390174bff85f65c2c2e558a9e

#include <boost/container/vector.hpp>
#include <utility>
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
  using size_type = Mesh::size_type;
  using dim_t = base::RefEl::dim_t;

  /**
   * @brief construct data vector indexed by entities of a particular
   * co-dimension
   * @param mesh pointer to underlying mesh
   * @param codim co-dimension of indexing entities for the given mesh
   *
   * @note This constructor does not initialize the data that is attached to the
   * entities. They are initialized arbitrarily!
   */
  CodimMeshDataSet(std::shared_ptr<const Mesh> mesh, dim_t codim)
      : MeshDataSet<T>(),
        mesh_(std::move(mesh)),
        data_(mesh_->NumEntities(codim)),
        codim_(codim) {}

  /**
   * @brief construct and initialize data vector indexed by entities of a
   *        particula co-dimension
   * @param mesh pointer to underlying mesh
   * @param codim co-dimension of indexing entities for the given mesh
   * @param init value to be copied into all data elements
   */
  template <class = typename std::enable_if<
                std::is_copy_constructible<T>::value>::type>
  CodimMeshDataSet(std::shared_ptr<const Mesh> mesh, dim_t codim, T init)
      : MeshDataSet<T>(),
        mesh_(std::move(mesh)),
        data_(mesh_->NumEntities(codim), init),
        codim_(codim) {}

  /**
   * @brief Get a (modifiable) reference to the data stored with entity e.
   * @param e The entity whose data should be retrieved/modified
   * @return  A reference to the stored data.
   *
   * @note The behavior of this method is undefined if `DefinedOn(e) == false`!
   */
  [[nodiscard]] T& operator()(const Entity& e) {
    LF_ASSERT_MSG(DefinedOn(e), "MeshDataSet is not defined on this entity.");
    return data_[mesh_->Index(e)];
  }
  // NOLINTNEXTLINE(readability-const-return-type)
  [[nodiscard]] const T& operator()(const Entity& e) const override {
    LF_ASSERT_MSG(DefinedOn(e), "MeshDataSet not defined on this entity.");
    return data_[mesh_->Index(e)];
  }
  [[nodiscard]] bool DefinedOn(const Entity& e) const override {
    return e.Codim() == codim_ && mesh_->Contains(e);
  }

 private:
  /** pointer to underlying mesh. Owned by the data set ! */
  std::shared_ptr<const Mesh> mesh_;
  /** data vector */
  boost::container::vector<T> data_;
  /** co-dimension */
  dim_t codim_;
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
    const std::shared_ptr<const lf::mesh::Mesh>& mesh, base::dim_t codim) {
  using impl_t = CodimMeshDataSet<T>;
  return std::shared_ptr<impl_t>(new impl_t(mesh, codim));
}

/**
 * @brief Create a new CodimMeshDataSet that attaches data of type `T` with
 *        every entity with codimension `codim`. The data of every entity
 *        is initialized to the given value (`T` must be copyable!)
 * @param mesh The mesh that contains the entities.
 * @param codim The codimension of the entities with which the data is stored.
 * @param init The initial value that should be assigned to every entity.
 */
template <class T, class = typename std::enable_if<
                       std::is_copy_constructible<T>::value>::type>
std::shared_ptr<CodimMeshDataSet<T>> make_CodimMeshDataSet(
    const std::shared_ptr<const Mesh>& mesh, base::dim_t codim, T init) {
  using impl_t = CodimMeshDataSet<T>;
  return std::shared_ptr<impl_t>(new impl_t(mesh, codim, init));
}

}  // namespace lf::mesh::utils

#endif  // __bae24f2390174bff85f65c2c2e558a9e
