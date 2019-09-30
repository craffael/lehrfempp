#ifndef THESIS_MESH_MINIMAL_MESH_BUILDER_H
#define THESIS_MESH_MINIMAL_MESH_BUILDER_H

#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/structured_mesh_builder.h>

namespace projects::ipdg_stokes {

namespace mesh {

class MinimalMeshBuilder : public lf::mesh::hybrid2d::StructuredMeshBuilder {
 public:
  explicit MinimalMeshBuilder(
      std::shared_ptr<lf::mesh::MeshFactory> mesh_factory)
      : lf::mesh::hybrid2d::StructuredMeshBuilder(std::move(mesh_factory)) {
    // Nothing to do here
  }

  std::shared_ptr<lf::mesh::Mesh> Build() override;
};

}  // end namespace mesh

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_MESH_MINIMAL_MESH_BUILDER_H
