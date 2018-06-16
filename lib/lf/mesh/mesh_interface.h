
#ifndef __b86b0c9cb0fd48da931a1f24421b8842
#define __b86b0c9cb0fd48da931a1f24421b8842

#include <lf/base/base.h>

#include "entity.h"

namespace lf::mesh {
/**
 * @brief Abstract interface for objects representing a single mesh
 *
 * This abstract base class desccribes the basic functionality of objects
 * that manage single-level conforming finite element meshes. These objects
 * essentially boil down to containers for mesh entities of different
 * co-dimensions. Thus they allow sequential traversal of these entities.
 *
 * Another important functionality concerns the management of entity indices,
 * which have to provide a consecutive numbering of entities of a specific
 * co-dimension starting from zero.
 */
class Mesh {
 protected:
  Mesh() = default;
  Mesh(const Mesh&) = default;
  Mesh(Mesh&&) = default;
  Mesh& operator=(const Mesh&) = default;
  Mesh& operator=(Mesh&&) = default;

 public:
  using size_type = unsigned int;
  /**
   * @brief The dimension of the manifold described by the mesh, or
   *        equivalently the maximum dimension of the reference elements
   *        present in the mesh.
   */
  virtual char DimMesh() const = 0;

  /**
   * @brief The dimension of the Euclidean space in which the mesh is
   *        embedded.
   */
  virtual char DimWorld() const = 0;

  /**
   * @brief All entities of a given codimension.
   * @param codim The codimension of the entities that should be returned.
   * @return A base::ForwardRange that can be used to traverse the entities.
   *
   * @sa Entity
   *
   * Principal access method for entities distinguished only by their
   * co-dimension Hence, all cells of a mesh are covered by the range returned
   * when giving co-dimension 0, regardless of their concrete shape.
   */
  virtual base::ForwardRange<const Entity> Entities(char codim) const = 0;

  /**
   * @brief The number of Entities that have the given codimension.
   * @param codim The codimension of the entities that should be counted.
   * @return That number of entities that have the given codimension.
   */
  virtual size_type Size(char codim) const = 0;

  /**
   * @brief Acess to the index of a mesh entity of any co-dimension
   * @param e Entity whose index is requested
   * @return index ranging from 0 to no. of entities of the same co-dimension-1
   *
   * It is a strict convention in LehrFEM++ that all entities of the same
   * co-dimension belonging to a mesh are endowed with an integer index. These
   * indices are guaranteed to be contiguous and to range from 0 to
   * `Size(codim)-1`.
   * @note The index of a mesh entity is NOT related to its position in the
   * range returned by the Entities() method.
   */
  virtual size_type Index(const Entity& e) const = 0;

  /**
   * @brief virtual destructor
   */
  virtual ~Mesh() = default;
};
}  // namespace lf::mesh

#endif  // __b86b0c9cb0fd48da931a1f24421b8842
