/**
 * @file
 * @brief Abstract base class for generators of lf::mesh::Mesh objects
 * @author Raffael Casagrande
 * @date   2018
 * @copyright MIT License
 */
#ifndef __5ac8f981f27e45d3b9d15fc9d52f7136
#define __5ac8f981f27e45d3b9d15fc9d52f7136

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>

#include "mesh_interface.h"

namespace lf::mesh {

class MeshFactory {
 protected:
  MeshFactory() = default;
  MeshFactory(const MeshFactory&) = default;
  MeshFactory(MeshFactory&&) = default;
  MeshFactory& operator=(const MeshFactory&) = default;
  MeshFactory& operator=(MeshFactory&&) = default;

 public:
  /** @copydoc Mesh::size_type */
  using size_type = unsigned int;

  /**
   * @brief Coordinate type of a point.
   */
  using coord_t = Eigen::VectorXd;

  using dim_t = base::dim_t;

  /**
   * @brief Return the Mesh::DimWorld() of the mesh that will be returned.
   */
  [[nodiscard]] virtual dim_t DimWorld() const = 0;

  /**
   * @brief Return the Mesh::DimMesh() of the mesh that will be returned.
   */
  [[nodiscard]] virtual dim_t DimMesh() const = 0;

  /**
   * @brief Add a point to the mesh.
   * @param coord The coordinate of the point (should have DimWorld() rows)
   * @return The 0-based index of the entity that will be created.
   *         The first call to this method will return 0, the second call 1,
   *         ...
   */
  // NOLINTNEXTLINE(modernize-use-nodiscard)
  virtual size_type AddPoint(coord_t coord) = 0;

  /**
   * @brief Add a point to the mesh.
   * @param geometry unique pointer to geometry object giving location of point
   * @return The 0-based index of the entity that will be created.
   *         The first call to this method will return 0, the second call 1,
   *         ...
   */
  // NOLINTNEXTLINE(modernize-use-nodiscard)
  virtual size_type AddPoint(
      std::unique_ptr<geometry::Geometry>&& geometry) = 0;

  /**
   * @brief Add an an entity (codim>0) to the mesh.
   * @param ref_el The reference element of the entity.
   * @param nodes The 0-based indices of the nodes that make up this entity
   *              (as returned from AddPoint()). This range should contain
   *              exactly `ref_el.NumNodes()`.
   * @param geometry The geometric description of the mesh entity.
   *                 Can be a nullptr (see below).
   * @return The index of the entity that will be created. It is guaranteed
   *         that the indices are consecutive (per dimension).
   *
   * Use this method to add mesh elements (codim=0) and optionally entities
   * with `codim>0` if this entity
   * - should use a particular geometry object that cannot be deduced from its
   *   super-entities, or
   * - if you need to know the index of this entity at a later stage, e.g.
   *   to assign flags.
   *
   * #### If `geometry == nullptr`
   * If the geometry object is not specified, the mesh will try to
   * deduce the geometry from super entities (i.e. entities that contain
   * the given entity as a sub-entity). E.g. if you add an entity with
   * `ref_el=RefEl::kSegment` to a mesh with `dimMesh=2` and you don't
   * specify a `geometry`, then the geometry of the entity will be deduced
   * from any triangle/quadrilateral that contains this entity
   * (using the method geometry::Geometry::SubGeometry()).
   * If the MeshFactory cannot deduce the geometry from a father entity,
   * an error will be raised when `Build()` is called.
   *
   * Since the geometry object is taken from a super entity, there is no
   * guarantee that the node order of the created entity and the one passed
   * trough the parameter `nodes` are the same!
   *
   * #### If `geometry != nullptr`
   * If you specify a geometry explicitly, the created entity will have exactly
   * this geometry object as a reference. Hence the node ordering of the
   * created entity object agrees with the ordering of `nodes`.
   *
   * @note The node indices passed with the parameter `nodes` must have been
   * obtained with AddPoint() before calling this method.
   */
  // NOLINTNEXTLINE(modernize-use-nodiscard)
  virtual size_type AddEntity(
      base::RefEl ref_el, const nonstd::span<const size_type>& nodes,
      std::unique_ptr<geometry::Geometry>&& geometry) = 0;

  /**
   * @brief Construct a mesh out of the specified nodes and elements.
   * @return The created mesh.
   *
   * @note All data supplied to the MeshFactory will be cleared after
   *       a mesh has been built successfully. I.e. the mesh factory can
   *       be used to construct another mesh after calling `Build()`.
   */
  [[nodiscard]] virtual std::shared_ptr<Mesh> Build() = 0;

  /// @brief Virtual destructor.
  virtual ~MeshFactory() = default;
};

}  // namespace lf::mesh

#endif  // __5ac8f981f27e45d3b9d15fc9d52f7136
