#ifndef INCG37e385afbd3b4b1dba8611fb71787822
#define INCG37e385afbd3b4b1dba8611fb71787822

#include <fmt/ostream.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>

namespace lf::mesh {

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

/** @brief **Interface class** representing a topological entity in a
 * cellular complex
 *
 * Example: a 2D hybrid mesh consists of cells (entities of co-dimension 0),
 * edges (entities of co-dimension 1), and nodes (entities of co-dimension 2).
 * The set of these entities endowed with _incidence relations_ defines the
 * topology of the mesh.
 *
 * The core purpose of this class is to provide interface definitions for
 * accessing incidence relations and geometry information. This interface
 * applies  to all (topological) types of entities.
 *
 * Further information can be found in @lref{sss:meshtopo} of the
 * [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
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
  [[nodiscard]] virtual unsigned Codim() const = 0;

  /**
   * @brief Return all sub entities of this entity that have the given
   *        codimension (w.r.t. this entity!)
   * @param rel_codim The _relative co-dimension_ w.r.t. this entity
   * @return A [span](https://en.cppreference.com/w/cpp/container/span) of
   * pointers to the sub-entities with the specified _relative co-dimension_
   *
   * Implicitly this function defines the numbering of sub-entities, see @ref
   * lf::base::RefEl and [Lecture
   * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
   * @lref{par:subentnum}.
   *
   * @note For a mesh covering a manifold of dimension 2, we have the following
   * cases
   * - For a cell (co-dimension 0 entity), the cell itself is a subentity
   *   of relative co-dimension 0, the edges have relative co-dimension 1, and
   *   the vertices relative co-dimension 2: in this case the usual co-dimension
   *   agrees with the relative co-dimension.
   * - For an edge (co-dimension 1 entity), the edge itself is the only
   *   sub-entity with relative co-dimension 0, and the endpoints are the
   *   sub-entitities of relative co-dimension 1.
   *
   * @note The lifetime of the returned span equals the lifetime of the Parent
   * Entity.
   *
   * #### Demonstration of usage
   * @snippet entityuse.cc usage
   *
   * Use of this method is also shown in [Lecture
   * Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
   * @lref{ex:subent}.
   */
  [[nodiscard]] virtual nonstd::span<const Entity* const> SubEntities(
      unsigned rel_codim) const = 0;

  /**
   * @brief return [span](https://en.cppreference.com/w/cpp/container/span) of
   * relative orientations of sub-entities of the next higher co-dimension.
   *
   * The corners of every entity are numbered and, thus, define a _local
   * orientation_ of the sub-entitites of co-dimension+1 as explained in
   * [Lecture Document](https://www.sam.math.ethz.ch/~grsam/NUMPDEFL/NUMPDE.pdf)
   * @lref{rem:ori}. This local orientation need not agree with the _intrinsic
   * orientation_ of that sub-entity as given by its own numbering of corners.
   * Such a possible mismatch is detected by this function.
   */
  [[nodiscard]] virtual nonstd::span<const Orientation> RelativeOrientations()
      const = 0;

  /**
   * @brief Describes the geometry of this entity.
   * @return A _pointer_ to a Geometry object that will remain valid for as long
   *         as the Mesh remains valid.
   *
   * Why does this member function return a pointer instead of a reference? One
   * reason is that entities without a geometric shape are conceivable as
   * building block of a "mesh graph", where only topological information
   * matters. Another reason is that during the construction of a mesh, it turns
   * out to be convenient to build entities "without geometry" first and then
   * endow them with geometric information. A `nullptr` is a good way to
   * indicate missing geometric information.
   *
   */
  [[nodiscard]] virtual const geometry::Geometry* Geometry() const = 0;

  /**
   * @brief Describes the reference element type of this entity.
   * @return An object of type @ref lf::base::RefEl.
   */
  [[nodiscard]] virtual base::RefEl RefEl() const = 0;

  /**
   * @brief Check if two entities are the same
   * @param rhs Check if this entity is the same as the rhs entity.
   *
   * @note The behavior of this method is undefined if the rhs entity belongs
           to a different Mesh.
   */
  [[nodiscard]] virtual bool operator==(const Entity& rhs) const = 0;

  /**
   * @brief Check if two entities are different.
   * @sa Entity::operator==
   */
  [[nodiscard]] bool operator!=(const Entity& rhs) const {
    return !operator==(rhs);
  }

  /**
   * @brief Virtual Destructor.
   */
  virtual ~Entity() = default;

};  // class entity

/**
 * @brief Operator overload to print the reference element of `Entity` to a
 * stream, such as `std::cout`.
 *
 * @param stream The stream to which this function should output
 * @param entity The entity to write to `stream`.
 * @return The stream itself.
 *
 */
std::ostream& operator<<(std::ostream& stream, const lf::mesh::Entity& entity);

}  // namespace lf::mesh

/**
 * @brief Make lf::mesh::Entity formattable with fmt, see
 * https://fmt.dev/latest/api.html#ostream-api
 */
template <> struct fmt::formatter<lf::mesh::Entity> : ostream_formatter {};

#endif  // INCG37e385afbd3b4b1dba8611fb71787822
