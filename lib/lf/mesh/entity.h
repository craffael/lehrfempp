#ifndef __37e385afbd3b4b1dba8611fb71787822
#define __37e385afbd3b4b1dba8611fb71787822
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>

namespace lf::mesh {

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
   * @param codim The codim w.r.t. this entity
   * @return
   */
  virtual base::RandomAccessRange<const Entity> SubEntities(
      char codim) const = 0;

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
   * @note The behavior of this method is undefined if the rhs entity belongs to
   * a different Mesh.
   */
  virtual bool operator==(const Entity& rhs) const = 0;

  /**
   * @brief Check if two entities are different.
   */
  bool operator!=(const Entity& rhs) const { return !operator==(rhs); }

  /**
   * @brief Virtual Destructor.
   */
  virtual ~Entity() = default;
};

}  // namespace lf::mesh

#endif  // __37e385afbd3b4b1dba8611fb71787822
