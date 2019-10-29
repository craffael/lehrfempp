#ifndef TP_TRIAG_MESH_BUILDER_H
#define TP_TRIAG_MESH_BUILDER_H

#include "lf/mesh/mesh_factory.h"
#include "structured_mesh_builder.h"

namespace lf::mesh::hybrid2d {

/**
 * @brief Implements a `MeshBuilder` that generates a triangular structured mesh
 *
 * Mesh generator for a triangular tensor product mesh covering a rectangle.
 *
 * A triangular tensor product mesh of a rectangular domain is built by
 * subdividing the domain into equal squares and splitting each of them
 * along the diagonal.
 *
 * The nodes and cells of the resulting mesh are numbered lexicographically
 * from bottom left to  top right. The super-diagonal cells have even indices
 * the sub-diagonal cells odd indices.
 *
 * The horizontal edges are numbered first, the vertical edges next; both groups
 * lexicographically.
 */
class TPTriagMeshBuilder : public StructuredMeshBuilder {
 public:
  /**
   * @brief Constructor: set factory object to be used by the builder
   *
   */
  explicit TPTriagMeshBuilder(std::unique_ptr<mesh::MeshFactory> mesh_factory)
      : StructuredMeshBuilder(std::move(mesh_factory)) {
    LF_ASSERT_MSG(
        mesh_factory_->DimWorld() == 2,
        "TPTriagMeshBuilder can only construct meshes with DimWorld==2");
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
};  // end class definition TPTriagMeshBuilder

}  // namespace lf::mesh::hybrid2d

#endif /* TP_TRIAG_MESH_BUILDER_H */
