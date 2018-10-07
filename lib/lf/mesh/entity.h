#ifndef __37e385afbd3b4b1dba8611fb71787822
#define __37e385afbd3b4b1dba8611fb71787822

#include <lf/base/base.h>
#include <lf/base/static_vars.h>
#include <lf/geometry/geometry.h>

namespace lf::mesh {
/** type for length of arrays/vectors */
using size_type = lf::base::size_type;
/** type for (co-)dimensions */
using dim_t = lf::base::dim_t;
/** type for local index numbers on the entity level */
using sub_idx_t = lf::base::sub_idx_t;
/** type for global index numbers referring to the mesh */
using glb_idx_t = lf::base::glb_idx_t;

/** @brief Relative orientation of a sub-entity
 *
 * Relative orientation is a fundamental topological concept describing the
 * relationship of of an entity and its subentities.
 *
 * For 2D meshes it is mainly relevant for edges: If the intrinsic orientation
 * of an edge defined by the ordering of its vertices agrees with the local
 * orientation of that edge in a cell, its relative orientation is positive,
 * otherwise negative.
 */
enum class Orientation : int { positive = 1, negative = -1 };

int to_sign(Orientation o);
char to_char(Orientation o);

/** @brief Representation of a topological entity in a cellular complex
 *
 * Example: a 2D hybrid mesh consists of cells (entities of co-dimension 0),
 * edges (entities of co-dimension 1), and nodes (entities of co-dimension 2).
 * The set of these entities endowed with _incidence relations_ defines the
 * topology of the mesh.
 *
 * The core purpose of this class is to provide methods for accessing incidence
 * relations.
 */
class Entity {
 protected:
  Entity() = default;
  Entity(const Entity&) = default;
  Entity(Entity&&) = default;
  Entity& operator=(const Entity&) = default;
  Entity& operator=(Entity&&) = default;

 public:
  /**
   * @brief The codimension of this entity w.r.t. the Mesh.dimMesh()
   * of the owning mesh manager.
   */
  virtual char Codim() const = 0;

  /**
   * @brief Return all sub entities of this entity that have the given
   *        codimension (w.r.t. this entity!)
   * @param rel_codim The _relative co-dimension_ w.r.t. this entity
   * @return a range containing all subentities of the specified
             _relative co-dimension_

   * @note For a mesh covering a manifold of dimension 2, we have the following
   cases
     - For a cell (co-dimension 0 entity), the cell itself is a subentity
       of relative co-dimension 0, the edges have relative co-dimension 1, and
   the vertices relative co-dimension 2: in this case the usual co-dimension
   agrees with the relative co-dimension.
     - For an edge (co-dimension 1 entity), the edge itself is the only
   sub-entity with relative co-dimension 0, and the endpoints are the
   sub-entitities of relative co-dimension 1.
   */
  virtual base::RandomAccessRange<const Entity> SubEntities(
      char rel_codim) const = 0;

  /**
   * @brief return array of relative orientations of sub-entities
   *        of the next hight co-dimension.
   */
  virtual base::RandomAccessRange<const Orientation> RelativeOrientations()
      const = 0;

  /**
   * @brief Describes the geometry of this entity.
   * @return A pointer to a Geometry object that will remain valid for as long
   *         as the Mesh remains valid.
   */
  virtual geometry::Geometry* Geometry() const = 0;

  /**
   * @brief Describes the reference element type of this entity.
   * @return An object of type base::RefEl.
   */
  virtual base::RefEl RefEl() const = 0;

  /**
   * @brief Check if two entities are the same
   * @param rhs Check if this entity is the same as the rhs entity.
   *
   * @note The behavior of this method is undefined if the rhs entity belongs
           to a different Mesh.
   */
  virtual bool operator==(const Entity& rhs) const = 0;

  /**
   * @brief Check if two entities are different.
   * @sa Entity::operator==
   */
  bool operator!=(const Entity& rhs) const { return !operator==(rhs); }

  /**
   * @brief Virtual Destructor.
   */
  virtual ~Entity() = default;

  // Add global output control
  /** @brief Diagnostics control variable */
  static int output_ctrl_;

};  // class entity

/**
 * @brief Operator overload to print a `Entity` to a stream, such as `std::cout`
 * @param stream The stream to which this function should output
 * @param entity The entity to write to `stream`.
 * @return The stream itself.
 *
 * - If Entity::output_ctrl_ == 0, type of reference element of entity is sent
 * as output to stream
 * - If Entity::output_ctrl_ > 0, then lf::mesh::utils::PrintInfo(const
 * lf::mesh::Entity& e, std::ostream& stream) is called.
 *
 */
std::ostream& operator<<(std::ostream& stream, const lf::mesh::Entity& entity);

}  // namespace lf::mesh

#endif  // __37e385afbd3b4b1dba8611fb71787822
