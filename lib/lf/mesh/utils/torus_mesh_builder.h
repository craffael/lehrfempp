/**
 * @file
 * @brief Declares a mesh builder for a tensor product grid of a torus
 * @author Anian Ruoss
 * @date   2018-10-22 16:06:17
 * @copyright MIT License
 */

#ifndef TORUS_MESH_BUILDER_H
#define TORUS_MESH_BUILDER_H

#include "lf/mesh/mesh_factory.h"
#include "structured_mesh_builder.h"

namespace lf::mesh::hybrid2d {

/**
 * @brief Implements a MeshBuilder for a tensor product grid of a torus
 *
 * The torus is built from a tensor product mesh covering a rectangle with
 * opposite sides identified. The major radius of the torus is proportional
 * to the rectangle's width (x-axis), the minor radius is proportional to the
 * rectangle's height (y-axis). The parametrization is described here:
 * https://en.wikipedia.org/wiki/Torus#Geometry
 *
 */
class TorusMeshBuilder : public StructuredMeshBuilder {
 public:
  /**
   * @brief Constructor: set factory object to be used by the builder
   */
  explicit TorusMeshBuilder(std::shared_ptr<mesh::MeshFactory> mesh_factory)
      : StructuredMeshBuilder(std::move(mesh_factory)) {
    LF_ASSERT_MSG(
        mesh_factory_->DimWorld() == 3,
        "TorusMeshBuilder can only construct meshes with DimWorld==3");
  }

  /**
   * @brief actual construction of the mesh
   */
  std::shared_ptr<mesh::Mesh> Build() override;

 private:
  /**
   * @brief vertex index from grid position
   */
  size_type VertexIndex(size_type i, size_type j) const {
    return i + j * no_of_x_cells_;
  }

 public:
  /** Diagnostics control variable */
  static unsigned int output_ctrl_;
};  // end class definition TorusMeshBuilder

}  // namespace lf::mesh::hybrid2d

#endif /* TORUS_MESH_BUILDER_H */
