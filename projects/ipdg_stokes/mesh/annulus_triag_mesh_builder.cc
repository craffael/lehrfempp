#include "annulus_triag_mesh_builder.h"

#include <lf/base/base.h>
#include <lf/base/forward_range.h>
#include <lf/geometry/geometry.h>

#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

namespace projects::ipdg_stokes {

namespace mesh {

std::shared_ptr<lf::mesh::Mesh> AnnulusTriagMeshBuilder::Build() {
  const lf::base::size_type nx = no_of_x_cells_;
  const lf::base::size_type ny = no_of_y_cells_;
  const lf::base::size_type vertex_count = nx * (ny + 1);
  const double r = bottom_left_corner_[1];
  const double R = top_right_corner_[1];
  const double dr = (R - r) / ny;
  const double dphi = 2 * M_PI / nx;
  auto node_position = [&](lf::base::size_type r_idx,
                           lf::base::size_type phi_idx) {
    Eigen::Vector2d pos;
    const double pos_r = r + r_idx * dr;
    const double pos_phi = phi_idx * dphi;
    pos << pos_r * std::cos(pos_phi), pos_r * std::sin(pos_phi);
    return pos;
  };

  // Add all vertices to the mesh
  std::vector<lf::base::size_type> v_idx(vertex_count);
  lf::base::size_type node_count = 0;
  for (lf::base::size_type r_idx = 0; r_idx < ny + 1; ++r_idx)
    for (lf::base::size_type phi_idx = 0; phi_idx < nx; ++phi_idx)
      v_idx[node_count++] =
          mesh_factory_->AddPoint(node_position(r_idx, phi_idx));

  // Construct the triangles
  const lf::base::RefEl ref_el = lf::base::RefEl::kTria();
  for (lf::base::size_type r_idx = 0; r_idx < ny; ++r_idx) {
    for (lf::base::size_type phi_idx = 0; phi_idx < nx; ++phi_idx) {
      const lf::base::size_type v1 = v_idx[nx * r_idx + phi_idx];
      const lf::base::size_type v2 = v_idx[nx * (r_idx + 1) + phi_idx];
      const lf::base::size_type v3 =
          v_idx[nx * (r_idx + 1) + ((phi_idx + 1) % nx)];
      const lf::base::size_type v4 = v_idx[nx * r_idx + ((phi_idx + 1) % nx)];
      const lf::base::ForwardRange<const lf::base::size_type> trig1(
          {v1, v2, v4});
      const lf::base::ForwardRange<const lf::base::size_type> trig2(
          {v2, v3, v4});
      Eigen::Matrix<double, 2, 3> verts1;
      Eigen::Matrix<double, 2, 3> verts2;
      verts1 << node_position(r_idx, phi_idx),
          node_position(r_idx + 1, phi_idx),
          node_position(r_idx, (phi_idx + 1) % nx);
      verts2 << node_position(r_idx + 1, phi_idx),
          node_position(r_idx + 1, (phi_idx + 1) % nx),
          node_position(r_idx, (phi_idx + 1) % nx);
      std::unique_ptr<lf::geometry::Geometry> geom1 =
          std::make_unique<lf::geometry::TriaO1>(verts1);
      std::unique_ptr<lf::geometry::Geometry> geom2 =
          std::make_unique<lf::geometry::TriaO1>(verts2);
      mesh_factory_->AddEntity(ref_el, trig1, std::move(geom1));
      mesh_factory_->AddEntity(ref_el, trig2, std::move(geom2));
    }
  }

  return mesh_factory_->Build();
}

}  // end namespace mesh

}  // end namespace projects::ipdg_stokes
