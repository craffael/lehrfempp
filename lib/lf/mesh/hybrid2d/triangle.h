/**
 * @file
 * @brief Entity implementation for triangles
 * @author Raffael Casagrande
 * @date   2018-06-22 03:59:23
 * @copyright MIT License
 */

#ifndef __6f934bca210e4020914d12e0910fab42
#define __6f934bca210e4020914d12e0910fab42

#include <lf/mesh/mesh.h>

namespace lf::mesh::hybrid2d {

class Point;
class Segment;

/**
 * @brief Describes a trilateral cell for a 2D hybrid mesh
 *
 * A trilateral cell is defined by ordered lists of references to its nodes
 * and its edges; internal consistency is required
 * @note Every `Segment` object owns a smart pointer to an associated geometry
 * object.
 *
 */
class Triangle : public mesh::Entity {
  using size_type = mesh::Mesh::size_type;

 public:
  /** @brief default constructors, needed by std::vector */
  Triangle() = default;

  /** @name Default and disabled constructors
   * @{ */
  Triangle(const Triangle&) = delete;
  Triangle(Triangle&&) noexcept = default;
  Triangle& operator=(const Triangle&) = delete;
  Triangle& operator=(Triangle&&) noexcept = default;
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
   * @param edge0 pointer to first edge
   * @param edge1 pointer to second edge
   * @param edge2 pointer to third edge
   *
   * The sub-entities have to be consistent according to the conventions
   * fixed for a reference element of type `kTria`.
   * @note Note that you need to create a suitable geometry object for the
   * entity before you can initialize the entity object itseld.
   */
  explicit Triangle(size_type index,
                    std::unique_ptr<geometry::Geometry>&& geometry,
                    const Point* corner0, const Point* corner1,
                    const Point* corner2, const Segment* edge0,
                    const Segment* edge1, const Segment* edge2);

  /** @brief an edge is an entity of co-dimension 1 */
  [[nodiscard]] unsigned Codim() const override { return 0; }

  /** @brief access to index of an entity */
  [[nodiscard]] size_type index() const { return index_; }

  /** @brief Access to all subentities selected by **relative** co-dimension
   * @param rel_codim if 1 select edges, if 2 select nodes, if 0 select cell
   itself
   * @return
     - for rel_codim == 1: return 3-range covering corners
     - for rel_codim == 2: return 3-range containing the edges
   */
  [[nodiscard]] nonstd::span<const Entity* const> SubEntities(
      unsigned rel_codim) const override;

  /** @brief Access to relative orientations of edges
   * @sa mesh::Orientation
   */
  [[nodiscard]] base::RandomAccessRange<const lf::mesh::Orientation>
  RelativeOrientations() const override {
    return base::RandomAccessRange<const lf::mesh::Orientation>(
        edge_ori_.begin(), edge_ori_.end());
  }

  /** @name Standard methods of an Entity object
   * @sa mesh::Entity
   * @{
   */
  [[nodiscard]] geometry::Geometry* Geometry() const override {
    return geometry_.get();
  }
  [[nodiscard]] base::RefEl RefEl() const override {
    return base::RefEl::kTria();
  }
  [[nodiscard]] bool operator==(const mesh::Entity& rhs) const override {
    return this == &rhs;
  }
  /** @} */

  ~Triangle() override = default;

 private:
  size_type index_ = -1;  // zero-based index of this entity.
  std::unique_ptr<geometry::Geometry> geometry_;  // shape information
  std::array<const Point*, 3> nodes_{};           // nodes = corners of cell
  std::array<const Segment*, 3> edges_{};         // edges of the cells
  std::array<lf::mesh::Orientation, 3>
      edge_ori_{};  // orientation of edges (set in constructor)
  Entity* this_;    // needed for SubEntity()
};

}  // namespace lf::mesh::hybrid2d

#endif  // __6f934bca210e4020914d12e0910fab42
