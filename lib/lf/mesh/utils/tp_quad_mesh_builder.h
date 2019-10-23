#ifndef TP_QUAD_MESH_BUILDER_H
#define TP_QUAD_MESH_BUILDER_H

#include "lf/mesh/mesh_factory.h"
#include "structured_mesh_builder.h"

namespace lf::mesh::hybrid2d {

/**
 * @brief Implements a Builder for a tensor product grid (with rectangular
 * cells)
 *
 * The generated grid will be uniform, which means that all cells are congruent.
 * The geometry of the rectangular domain and the number of cells can be
 * specified by setting Builder state parameters.
 *
 */
class TPQuadMeshBuilder : public StructuredMeshBuilder {
 public:
  /**
   * @brief Constructor: set factory object to be used by the builder
   *
   */
  explicit TPQuadMeshBuilder(std::unique_ptr<mesh::MeshFactory> mesh_factory)
      : StructuredMeshBuilder(std::move(mesh_factory)) {
    LF_ASSERT_MSG(
        mesh_factory_->DimWorld() == 2,
        "TPQuadMeshBuilder can only construct meshes with DimWorld==2");
  }

  /**
   * @brief actual construction of the mesh
   *
   */
  [[nodiscard]] std::shared_ptr<mesh::Mesh> Build() override;

 private:
  /**
   * @brief vertex index from grid position
   */
  [[nodiscard]] size_type VertexIndex(size_type i, size_type j) const {
    return (i + j * (no_of_x_cells_ + 1));
  }

 public:
  /** Diagnostics control variable */
  static unsigned int output_ctrl_;
};  // end class definition TPQuadMeshBuilder

}  // namespace lf::mesh::hybrid2d

#endif /* TP_QUAD_MESH_BUILDER_H */
