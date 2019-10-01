#include "minimal_mesh_builder.h"

#include <array>

namespace projects::ipdg_stokes::mesh {

std::shared_ptr<lf::mesh::Mesh> MinimalMeshBuilder::Build() {
  Eigen::Matrix<double, 2, 5> vertices;
  vertices << 0, 1, 0, -1, 0, 0, 0, 1, 0, -1;
  std::array<lf::base::size_type, 5> idx;

  // Add the vertices
  for (Eigen::Index i = 0; i < vertices.cols(); ++i) {
    idx[i] = mesh_factory_->AddPoint(vertices.col(i));
  }

  // Add the triangles
  auto add_triangle = [&](unsigned v1, unsigned v2, unsigned v3) {
    lf::base::ForwardRange<const lf::base::size_type> indices({v1, v2, v3});
    Eigen::Matrix<double, 2, 3> verts;
    verts << vertices.col(v1), vertices.col(v2), vertices.col(v3);
    std::unique_ptr<lf::geometry::Geometry> geom =
        std::make_unique<lf::geometry::TriaO1>(verts);
    mesh_factory_->AddEntity(lf::base::RefEl::kTria(), indices,
                             std::move(geom));
  };
  add_triangle(0, 1, 2);
  add_triangle(0, 2, 3);
  add_triangle(0, 3, 4);
  add_triangle(0, 4, 1);

  // Build the mesh
  return mesh_factory_->Build();
}

}  // end namespace projects::ipdg_stokes::mesh
