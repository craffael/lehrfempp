/**
 * @file
 * @brief Implementation of the 0-dimensional node entity for the hybrid2dp
 *        manager
 * @author Raffael Casagrande
 * @date   2018-06-22 03:54:24
 * @copyright MIT License
 */

#ifndef __818709b0104548a7b5e6f47bdba89f69
#define __818709b0104548a7b5e6f47bdba89f69

#include <lf/mesh/mesh.h>

namespace lf::mesh::hybrid2dp {

/**
 * @brief A node object for a 2D hybrid mesh
 *
 * @note Every `Entity` object owns a smart pointer to an associated geometry
 * object.
 *
 * Due to the unidirectional storage scheme for incidence information the node
 * object does not have much functionality, except for storing its index.
 */
class Point : public mesh::Entity {
  using size_type = mesh::Mesh::size_type;

 public:
  /** @brief default constructors, needed by std::vector */
  Point() = default;

  /** @ brief Default and disabled constructors
   * @{ */
  Point(const Point&) = delete;
  Point(Point&&) noexcept = default;
  Point& operator=(const Point&) = delete;
  Point& operator=(Point&&) noexcept = default;
  /** @} */

  /**
   * @brief constructor, is called from MeshFactory
   * @param index index of the entity to be created; will usually be
   * retrieved via the `Index()` method of `Mesh`
   * @param geometry pointer to a geometry object providing the shape of the
   * entity
   *
   * @note Note that you need to create a suitable geometry object for the
   * entity before you can initialize the entity object itseld.
   */
  explicit Point(size_type index,
                 std::unique_ptr<geometry::Geometry>&& geometry)
      : index_(index), geometry_(std::move(geometry)) {
    LF_VERIFY_MSG(geometry->DimLocal() == 0,
                  "Geometry must be that of a point");
    LF_VERIFY_MSG(geometry->RefEl() == base::RefEl::kPoint(),
                  "Geometry must fit point");
  }

  char Codim() const override { return 2; }

  base::RandomAccessRange<const mesh::Entity> SubEntities(
      char rel_codim) const override {
    LF_ASSERT_MSG(rel_codim == 0, "A point has only codim = 0 sub-entities");
    return base::RandomAccessRange<const mesh::Entity>(this, this + 1);
  }

  /** @brief return _pointer_ to associated geometry object */
  geometry::Geometry* Geometry() const override { return geometry_.get(); }

  base::RefEl RefEl() const override { return base::RefEl::kPoint(); }

  bool operator==(const mesh::Entity& rhs) const override {
    return this == &rhs;
  }

  ~Point() override = default;

 private:
  size_type index_ = -1;  // zero-based index of this entity.
  std::unique_ptr<geometry::Geometry> geometry_;  // shape information
};

}  // namespace lf::mesh::hybrid2dp

#endif  // __818709b0104548a7b5e6f47bdba89f69
