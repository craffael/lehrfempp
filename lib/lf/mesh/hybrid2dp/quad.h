/**
 * @file
 * @brief Entity implementation for Quadrilaterals
 * @author Raffael Casagrande
 * @date   2018-06-22 04:02:41
 * @copyright MIT License
 */

#ifndef __9dc22ac9eb6645d3a27f60a7abcd52a4
#define __9dc22ac9eb6645d3a27f60a7abcd52a4

#include <lf/mesh/mesh.h>

namespace lf::mesh::hybrid2dp {

class Point;
class Segment;

/**
 * @brief Describes a quadrilateral cell for a 2D hybrid mesh
 *
 * A quadrilateral cell stores and ordered list of four nodes and
 * of four edges, which have to be compatible.
 * @note Every `Segment` object owns a smart pointer to an associated geometry
 * object.
 *
 */
class Quadrilateral : public mesh::Entity {
  using size_type = mesh::Mesh::size_type;

 public:
  /** @brief default constructors, needed by std::vector */
  Quadrilateral() = default;

  /** @defgroup Default and disabled constructors
   * @{ */
  Quadrilateral(const Quadrilateral&) = delete;
  Quadrilateral(Quadrilateral&&) noexcept = default;
  Quadrilateral& operator=(const Quadrilateral&) = delete;
  Quadrilateral& operator=(Quadrilateral&&) noexcept = default;
  /** @} */

  /**
   * @brief constructor, is called from MeshFactory
   * @param index index of the entity to be created; will usually be
   * retrieved via the `Index()` method of `Mesh`
   * @param geometry pointer to a geometry object providing the shape of the
   * cell
   * @param corner0 pointer to first node
   * @param corner1 pointer to second node
   * @param corner2 pointer to third node
   * @param corner3 pointer to fourth node
   * @param edge0 pointer to first edge
   * @param edge1 pointer to second edge
   * @param edge2 pointer to third edge
   * @param edge3 pointer to fourth edge
   *
   * The sub-entities have to be consistent according to the conventions
   * fixed for a reference element of type `kTria`.
   * @note Note that you need to create a suitable geometry object for the
   * entity before you can initialize the entity object itseld.
   */
  explicit Quadrilateral(size_type index,
                         std::unique_ptr<geometry::Geometry>&& geometry,
                         const Point* corner0, const Point* corner1,
                         const Point* corner2, const Point* corner3,
                         const Segment* edge0, const Segment* edge1,
                         const Segment* edge2, const Segment* edge3)
      : index_(index),
        geometry_(std::move(geometry)),
        nodes_({corner0, corner1, corner2, corner3}),
        edges_({edge0, edge1, edge2, edge3}) {
    LF_VERIFY_MSG(corner0 != nullptr, "Invalid pointer to corner 0");
    LF_VERIFY_MSG(corner1 != nullptr, "Invalid pointer to corner 1");
    LF_VERIFY_MSG(corner2 != nullptr, "Invalid pointer to corner 2");
    LF_VERIFY_MSG(corner3 != nullptr, "Invalid pointer to corner 3");
    LF_VERIFY_MSG(edge0 != nullptr, "Invalid pointer to edge 0");
    LF_VERIFY_MSG(edge1 != nullptr, "Invalid pointer to edge 1");
    LF_VERIFY_MSG(edge2 != nullptr, "Invalid pointer to edge 2");
    LF_VERIFY_MSG(edge3 != nullptr, "Invalid pointer to edge 3");
    LF_VERIFY_MSG(geometry->DimLocal() == 2,
                  "Geometry must describe a 2D cell");
    LF_VERIFY_MSG(geometry->RefEl() == base::RefEl::kQuad(),
                  "Cell geometry must fit a quad");
    /*
       TODO: consistency check
    */
  }

  /** @brief an edge is an entity of co-dimension 1 */
  char Codim() const override { return 0; }

  /** @brief Access to all subentities selected by **relative** co-dimension
   * @param rel_codim if 1 select edges, if 2 select nodes, if 0 select cell
   itself
   * @return
     - for rel_codim == 1: return range with quadrilateral itself as only
   element
     - for rel_codim == 1: return 4-range covering corners
     - for rel_codim == 2: return 4-range containing the edges
   */
  base::RandomAccessRange<const mesh::Entity> SubEntities(
      char rel_codim) const override;

  /** @defgroup Standard methods of an Entity object
   * @sa mesh::Entity
   * @{
   */
  geometry::Geometry* Geometry() const override { return geometry_.get(); }
  base::RefEl RefEl() const override { return base::RefEl::kQuad(); }
  bool operator==(const mesh::Entity& rhs) const override {
    return this == &rhs;
  }
  /** @} */

  ~Quadrilateral() override = default;

 private:
  size_type index_ = -1;  // zero-based index of this entity.
  std::unique_ptr<geometry::Geometry> geometry_;  // shape information
  std::array<const Point*, 4> nodes_{};           // nodes = corners of quad
  std::array<const Segment*, 4> edges_{};         // edges of quad
};

}  // namespace lf::mesh::hybrid2dp

#endif  // __9dc22ac9eb6645d3a27f60a7abcd52a4
