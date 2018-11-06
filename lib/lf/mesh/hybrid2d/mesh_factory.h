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
   */
  explicit MeshFactory(dim_t dim_world) : dim_world_(dim_world) {}

  dim_t DimWorld() const override { return dim_world_; }

  dim_t DimMesh() const override { return 2; }

  size_type AddPoint(coord_t coord) override;

  size_type AddPoint(std::unique_ptr<geometry::Geometry>&& geometry) override;

  size_type AddEntity(base::RefEl ref_el,
                      const base::ForwardRange<const size_type>& nodes,
                      std::unique_ptr<geometry::Geometry>&& geometry) override;

  std::shared_ptr<mesh::Mesh> Build() override;

  /** @brief output function printing asssembled lists of entity information */
  void PrintLists(std::ostream& o = std::cout) const;

  ~MeshFactory() override = default;

 private:
  dim_t dim_world_;  // dimension of ambient space
  hybrid2d::Mesh::NodeCoordList nodes_;
  hybrid2d::Mesh::EdgeList edges_;
  hybrid2d::Mesh::CellList elements_;

 public:
  // Switch for verbosity level of output
  /** @brief Diagnostics control variable */
  static int output_ctrl_;
};

inline std::ostream& operator<<(std::ostream& stream,
                                const MeshFactory& /*mesh_factory*/) {
  stream << "mesh factory object";
  return stream;
}

}  // namespace lf::mesh::hybrid2d

#endif  // __e98a803fac5b430a8ff634ceb2f809aX
