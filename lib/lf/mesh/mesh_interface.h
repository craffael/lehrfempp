/**
 * @file
 * @brief Abstract base class for finite element meshes
 * @author Raffael Casagrande
 * @date   2018
 * @copyright MIT License
 */
#ifndef __b86b0c9cb0fd48da931a1f24421b8842
#define __b86b0c9cb0fd48da931a1f24421b8842

#include <lf/base/base.h>

#include "entity.h"

namespace lf::mesh {
/**
 * @brief Abstract interface for objects representing a single mesh
 *
 * This abstract base class describes the basic functionality of objects
 * that manage single-level conforming finite element meshes. These objects
 * essentially boil down to containers for mesh entities of different
 * co-dimensions. Thus they allow sequential traversal of these entities.
 *
 * Another important functionality concerns the management of entity indices,
 * which have to provide a consecutive numbering of entities of a specific
 * co-dimension starting from zero.
 *
 * #### Mesh object variables
 *
 * Mesh objects are usually accessed through [C++ *shared
 pointers*](https://en.cppreference.com/book/intro/smart_pointers); declare and
 initialize variable for a mesh as follows (second line optional)
 * ~~~
  std::shared_ptr<lf::mesh::Mesh> mesh_p =  .... ;
  lf::mesh::Mesh &mesh(*mesh_p);
 * ~~~
 * In most contexts mesh objects will be 'read-only' (immutable). Then the
 * declaration should look like
 * ~~~
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = .... ;
  const lf::mesh::Mesh &mesh(*mesh_p);
 * ~~~
 * Why do we use shared pointers? Usually in a finite element code a mesh object
 * is allocated dynamically in the beginning and disposed only before the code
 * terminates. Many other objects rely on information about the mesh. If the
 * mesh was destroyed before all the objects depending on it, those would be
 * plunged into an undefined state. Shared pointers ensure that a mesh object is
 * destroyed only after every object holding a pointer to it has been destroyed.
 *
 * #### Creating a mesh object
 *
 * The standard ways to create a @ref Mesh object are:
 * -# Direct generation via the @ref mesh::MeshFactory::Build() method.
 *    (This method is mainly used internally by LehrFEM++)
 * -# Calling LehrFEM++'s generator of test meshes @ref
 lf::mesh::test_utils::GenerateHybrid2DTestMesh().
 * -# Reading a mesh from file, see @ref lf::io::GmshReader.
 * -# Refining an existing mesh, see @ref lf::refinement::MeshHierarchy.
 *
 * #### Use case
 *
 * The following function demonstrates the use of all pertinent methods of
 * the class @ref Entity. In particular is shows how to loop over all entities
 * of a mesh.
 * @snippet meshuse.cc usage
 *
 */
class Mesh {
 protected:
  Mesh() = default;
  Mesh(const Mesh&) = default;
  Mesh(Mesh&&) = default;
  Mesh& operator=(const Mesh&) = default;
  Mesh& operator=(Mesh&&) = default;

 public:
  /** Auxiliary types */
  using size_type = lf::base::size_type;
  using dim_t = lf::base::dim_t;
  /**
   * @brief The dimension of the manifold described by the mesh, or
   *        equivalently the maximum dimension of the reference elements
   *        present in the mesh.
   */
  [[nodiscard]] virtual unsigned DimMesh() const = 0;

  /**
   * @brief The dimension of the Euclidean space in which the mesh is
   *        embedded.
   */
  [[nodiscard]] virtual unsigned DimWorld() const = 0;

  /**
   * @brief All entities of a given codimension.
   * @param codim The codimension of the entities that should be returned.
   * @return A [span](https://en.cppreference.com/w/cpp/container/span) of
   pointers to entities that enumerate all entities of the given codimension.
   *
   * @sa Entity
   *
   * Principal access method for entities distinguished only by their
   * co-dimension. Hence, all cells of a mesh are covered by the range returned
   * when giving co-dimension 0, regardless of their concrete shape.
   *
   * The typical loop for entity traversal looks like this, where `mesh` is a
   * variable containing a reference to a @ref Mesh object.
   * ~~~
     for (const lf::mesh::Entity* entity : mesh.Entities(codim)) {
       ....
     }
   * ~~~
   * Inside the loop body the variable `entity` contains a pointer to an
   * immutable object of type @ref Entity whose co-dimension is `codim`.
   *
   * @note The pointer remains valid for as long as the mesh data structure
   * remains valid.
   */
  [[nodiscard]] virtual nonstd::span<const Entity* const> Entities(
      unsigned codim) const = 0;

  /**
   * @brief The number of Entities that have the given codimension.
   * @param codim The codimension of the entities that should be counted.
   * @return That number of entities that have the given codimension.
   */
  [[nodiscard]] virtual size_type NumEntities(unsigned codim) const = 0;

  /**
   * @brief Tells number of entities of a particular topological/geometric type
   * @param ref_el_type topological/geometric type
   * @return number of entities of that type
   */
  [[nodiscard]] virtual size_type NumEntities(
      lf::base::RefEl ref_el_type) const = 0;

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
  [[nodiscard]] virtual size_type Index(const Entity& e) const = 0;

  /**
   * @brief Method for accessing an entity through its index
   * @param codim codimension of the entity. Remember that indices are supposed
   *        to be unique and contiguous for a given co-dimension
   * @param index an integer between 0 and number of entities of the given
   * co-dimension -1. It passes the index.
   * @return pointer to the entity object with the given index.
   *
   * Based on the bijectition between entities of a given co-dimension and an
   * integer range. The following expression should evaluate to `true`, if
   * `mesh` is a reference to a @ref Mesh object, and `idx` a valid index
   * ~~~
      ((idx < mesh.NumEntities(codim)) &&
   mesh.Index(*mesh.EntityByIndex(codim,idx)) == idx)
   * ~~~
   *
   * @note O(1) access complexity due to table lookup.
   */
  [[nodiscard]] virtual const mesh::Entity* EntityByIndex(
      dim_t codim, base::glb_idx_t index) const = 0;

  /**
   * @brief Check if the given entity is a part of this mesh.
   * @param e The entity that should be checked.
   * @return true if the entity belongs to the mesh.
   */
  [[nodiscard]] virtual bool Contains(const Entity& e) const = 0;

  /**
   * @brief virtual destructor
   */
  virtual ~Mesh() = default;
};

}  // namespace lf::mesh

#endif  // __b86b0c9cb0fd48da931a1f24421b8842
