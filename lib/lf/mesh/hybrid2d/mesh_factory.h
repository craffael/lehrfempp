#ifndef __e98a803fac5b430a8ff634ceb2f809aX
#define __e98a803fac5b430a8ff634ceb2f809aX

#include <lf/mesh/mesh.h>
#include "mesh.h"

#include <iostream>

namespace lf::mesh::hybrid2d {

/**
 * @brief Implements mesh::MeshFactory interface and can be used to construct
 *        a hybrid mesh with `dimMesh=2`.
 *
 * A planar triangular mesh with affine triangles as cells can be completely
 * specified by giving the list of vertex coordinates, an optional list
 * of edges, and a list of cells (triangles, quadrilaterals). All entities
 * can be supplied with a geometry. If this is missing, the mesh builder
 * tries to infer it from sub-entities or super-entities. If this is not
 * possible, an affine entity is built.
 */
class MeshFactory : public mesh::MeshFactory {
 public:
  MeshFactory(const MeshFactory&) = delete;
  MeshFactory(MeshFactory&&) = delete;
  MeshFactory& operator=(const MeshFactory&) = delete;
  MeshFactory& operator=(MeshFactory&&) = delete;

  /**
   * @brief Construct a new builder that can be used to construct a new hybrid2d
   *        mesh.
   * @param dim_world The dimension of the euclidean space in which the
   *                  mesh is embedded.
   * @param check_completeness If set to true, calling Build() will check that
   * the mesh is topologically complete. That means that every entity with
   * codimension `codim>0` is a subentity of at least one entity with
   * codimension `codim-1`. If `check_completeness = true` and the mesh is not
   * complete, an assert will fail.
   */
  explicit MeshFactory(dim_t dim_world, bool check_completeness = true)
      : dim_world_(dim_world), check_completeness_(check_completeness) {}

  [[nodiscard]] dim_t DimWorld() const override { return dim_world_; }

  [[nodiscard]] dim_t DimMesh() const override { return 2; }

  // NOLINTNEXTLINE(modernize-use-nodiscard)
  size_type AddPoint(coord_t coord) override;

  // NOLINTNEXTLINE(modernize-use-nodiscard)
  size_type AddPoint(std::unique_ptr<geometry::Geometry>&& geometry) override;

  // NOLINTNEXTLINE(modernize-use-nodiscard)
  size_type AddEntity(base::RefEl ref_el,
                      const nonstd::span<const size_type>& nodes,
                      std::unique_ptr<geometry::Geometry>&& geometry) override;

  [[nodiscard]] std::shared_ptr<mesh::Mesh> Build() override;

  /** @brief output function printing assembled lists of entity information */
  void PrintLists(std::ostream& o = std::cout) const;

  ~MeshFactory() override = default;

 private:
  dim_t dim_world_;  // dimension of ambient space
  hybrid2d::Mesh::NodeCoordList nodes_;
  hybrid2d::Mesh::EdgeList edges_;
  hybrid2d::Mesh::CellList elements_;

  // If set to true, the Build() method will check whether all sub-entities
  // belong to at least one entity */
  bool check_completeness_;

 public:
  /**
   * @brief logger that is used by the build method to output additional
   * information to the command line.
   */
  static std::shared_ptr<spdlog::logger>& Logger();
};

inline std::ostream& operator<<(std::ostream& stream,
                                const MeshFactory& /*mesh_factory*/) {
  stream << "mesh factory object";
  return stream;
}

}  // namespace lf::mesh::hybrid2d

#endif  // __e98a803fac5b430a8ff634ceb2f809aX
