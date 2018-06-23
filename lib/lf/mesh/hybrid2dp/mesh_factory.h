#ifndef __e98a803fac5b430a8ff634ceb2f809aX
#define __e98a803fac5b430a8ff634ceb2f809aX

#include <lf/mesh/mesh.h>
#include "mesh.h"

namespace lf::mesh::hybrid2dp {

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
  /**
   * @brief Construct a new builder that can be used to construct a new hybrid2d
   *        mesh.
   * @param dim_world The dimension of the euclidean space in which the
   *                  mesh is embedded.
   */
  explicit MeshFactory(dim_t dim_world)
      : dim_world_(dim_world), built_(false) {}

  /** @copydoc MeshFactory::DimWorld() */
  dim_t DimWorld() const override { return dim_world_; }

  /** @copydoc Mesh::DimMesh() */
  dim_t DimMesh() const override { return 2; }

  /** @copydoc MeshFactory::AddPoint() */
  size_type AddPoint(coord_t coord) override;

  /** @copydoc MeshFactory::AddEntity() */
  size_type AddEntity(base::RefEl ref_el,
                      const base::ForwardRange<const size_type>& nodes,
                      std::unique_ptr<geometry::Geometry>&& geometry) override;

  /** @copydoc MeshFactory::Build() */
  std::shared_ptr<mesh::Mesh> Build() override;

 private:
  dim_t dim_world_;  // dimension of ambient space
  bool built_;
  hybrid2dp::Mesh::NodeCoordList nodes_;
  hybrid2dp::Mesh::EdgeList edges_;
  hybrid2dp::Mesh::CellList elements_;
};

}  // namespace lf::mesh::hybrid2dp

#endif  // __e98a803fac5b430a8ff634ceb2f809aX
