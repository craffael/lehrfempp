#ifndef THESIS_MESH_ANNULUS_TRIAG_MESH_BUILDER_H
#define THESIS_MESH_ANNULUS_TRIAG_MESH_BUILDER_H

/**
 * @file annulus_triag_mesh_builder.h
 * @brief Constructs a mesh on an annulus
 */

#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/structured_mesh_builder.h>

namespace projects::ipdg_stokes {

namespace mesh {

/**
 * @brief A mesh builder for disks with a hole in the middle
 *
 * The diameter of the hole is controlled by the y-value passed to
 * `setBottomLeftCorner(x, y)` while the diameter of the disk is controlled by
 * the y-value passed to `setTopRightCorner(x, y)`.
 */
class AnnulusTriagMeshBuilder
    : public lf::mesh::hybrid2d::StructuredMeshBuilder {
 public:
  /**
   * @brief Constructor
   * @param mesh_factory A shared pointer to a mesh factory object used for the
   * assembly of the mesh
   */
  explicit AnnulusTriagMeshBuilder(
      std::shared_ptr<lf::mesh::MeshFactory> mesh_factory)
      : lf::mesh::hybrid2d::StructuredMeshBuilder(std::move(mesh_factory)) {
    // Nothing to do here
  }

  /**
   * @brief Build the mesh
   * @returns A shared pointer to a mesh in the form of a disk with a hole in
   * the middle
   */
  std::shared_ptr<lf::mesh::Mesh> Build() override;
};

}  // end namespace mesh

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_MESH_ANNULUS_TRIAG_MESH_BUILDER_H
