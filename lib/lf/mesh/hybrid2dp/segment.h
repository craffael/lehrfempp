/**
 * @file
 * @brief Implementation for 1d Segments for hybrid2dp mesh manager
 * @author Raffael Casagrande
 * @date   2018-06-22 03:56:25
 * @copyright MIT License
 */

#ifndef __feff908010fa4d75a9c006c02b4fafe7
#define __feff908010fa4d75a9c006c02b4fafe7

#include <lf/mesh/mesh.h>

namespace lf::mesh::hybrid2dp {

// Forward declaration:
class Point;

/**
 * @brief An edge object for a 2D hybrid mesh
 *
 * An topological edge object is define through two distinct references
 * to node objects of the mesh. Their ordering reflects the intrinsic
 * orientation of the mesh
 * @note Every `Segment` object owns a smart pointer to an associated geometry
 * object.
 *
 */
class Segment : public mesh::Entity {
  using size_type = mesh::Mesh::size_type;

 public:
  /** @brief default constructors, needed by std::vector */
  Segment() = default;

  /** @defgroup Default and disabled constructors
   * @{ */
  Segment(const Segment&) = delete;
  Segment(Segment&&) noexcept = default;
  Segment& operator=(const Segment&) = delete;
  Segment& operator=(Segment&&) noexcept = default;
  /** @} */

  /**
   * @brief constructor, is called from MeshFactory
   * @param index index of the entity to be created; will usually be
   * retrieved via the `Index()` method of `Mesh`
   * @param geometry pointer to a geometry object providing the shape of the
   * edge
   * @param endpoint0 pointer to the first node
   * @param endpoint1 pointer to the second node
   *
   * @note Note that you need to create a suitable geometry object for the
   * entity before you can initialize the entity object itseld.
   */
  explicit Segment(size_type index,
                   std::unique_ptr<geometry::Geometry>&& geometry,
                   const Point* endpoint0, const Point* endpoint1)
      : index_(index),
        geometry_(std::move(geometry)),
        nodes_({endpoint0, endpoint1}) {
    LF_VERIFY_MSG((endpoint0 != nullptr) && (endpoint1 != nullptr),
                  "Invalid pointer to endnode of edge");
    if (geometry_) {
      LF_VERIFY_MSG(geometry_->DimLocal() == 1,
                    "Geometry must describe a curve");
      LF_VERIFY_MSG(geometry_->RefEl() == base::RefEl::kSegment(),
                    "Segment geometry must fit a segment");
    }
  }

  /** @brief an edge is an entity of co-dimension 1 */
  char Codim() const override { return 1; }

  /** @brief Access to all subentities selected by **relative** co-dimension
   * @param rel_codim if 1 select endnodes, if 0 select edge itself
   * @return
     - for rel_codim == 1: return 2-range covering endnodes
     - for rel_codim == 0: return the Segment entity itself
   */
  base::RandomAccessRange<const mesh::Entity> SubEntities(
      char rel_codim) const override;

  /** @brief Access to relative orientations of endpoints
   *
   * This method just returns {+,-}, because points always have the intrinsic
   * orientation + and the orientation of an edge is defined through the
   * ordering of its vertices. */
  base::RandomAccessRange<const lf::mesh::Orientation> RelativeOrientations()
      const override {
    return base::RandomAccessRange<const lf::mesh::Orientation>(
        endpoint_ori_.begin(), endpoint_ori_.end());
  }

  /** @brief access to index of an entity */
  size_type index() const { return index_; }

  /** @defgroup Standard methods of an Entity object
   * @sa mesh::Entity
   * @{
   */
  geometry::Geometry* Geometry() const override { return geometry_.get(); }
  base::RefEl RefEl() const override { return base::RefEl::kSegment(); }
  bool operator==(const mesh::Entity& rhs) const override {
    return this == &rhs;
  }
  /** @} */

  ~Segment() override = default;

 private:
  size_type index_ = -1;  // zero-based index of this entity.
  std::unique_ptr<geometry::Geometry> geometry_;  // shape information
  std::array<const Point*, 2> nodes_{};           // nodes connected by edge
  static constexpr std::array<lf::mesh::Orientation, 2> endpoint_ori_{
      lf::mesh::Orientation::negative,
      lf::mesh::Orientation::positive};  // orientation of endpoints
};

}  // namespace lf::mesh::hybrid2dp

#endif  // __feff908010fa4d75a9c006c02b4fafe7
