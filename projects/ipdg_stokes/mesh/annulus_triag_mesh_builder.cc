#include "annulus_triag_mesh_builder.h"

#include <lf/base/base.h>
#include <lf/base/span.h>
#include <lf/geometry/geometry.h>

#include <array>
#include <cmath>
#include <vector>

namespace projects::ipdg_stokes::mesh {

std::shared_ptr<lf::mesh::Mesh> AnnulusTriagMeshBuilder::Build() {
  const lf::base::size_type vertex_count = num_angular_cells_ * (num_radial_cells_ + 1);
  const double r = inner_radius_;
  const double R = outer_radius_;
  const double dr = (R - r) / num_radial_cells_;
  const double dphi = 2 * M_PI / num_angular_cells_;
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
  for (lf::base::size_type r_idx = 0; r_idx < num_radial_cells_ + 1; ++r_idx) {
    for (lf::base::size_type phi_idx = 0; phi_idx < num_angular_cells_; ++phi_idx) {
      v_idx[node_count++] =
          mesh_factory_->AddPoint(node_position(r_idx, phi_idx));
    }
  }

  // Construct the triangles
  const lf::base::RefEl ref_el = lf::base::RefEl::kTria();
  for (lf::base::size_type r_idx = 0; r_idx < num_radial_cells_; ++r_idx) {
    for (lf::base::size_type phi_idx = 0; phi_idx < num_angular_cells_; ++phi_idx) {
      const lf::base::size_type v1 = v_idx[num_angular_cells_ * r_idx + phi_idx];
      const lf::base::size_type v2 = v_idx[num_angular_cells_ * (r_idx + 1) + phi_idx];
      const lf::base::size_type v3 =
          v_idx[num_angular_cells_ * (r_idx + 1) + ((phi_idx + 1) % num_angular_cells_)];
      const lf::base::size_type v4 = v_idx[num_angular_cells_ * r_idx + ((phi_idx + 1) % num_angular_cells_)];
      const std::array<lf::base::size_type, 3> trig1 = {v1, v2, v4};
      const std::array<lf::base::size_type, 3> trig2 = {v2, v3, v4};
      Eigen::Matrix<double, 2, 3> verts1;
      Eigen::Matrix<double, 2, 3> verts2;
      verts1 << node_position(r_idx, phi_idx),
          node_position(r_idx + 1, phi_idx),
          node_position(r_idx, (phi_idx + 1) % num_angular_cells_);
      verts2 << node_position(r_idx + 1, phi_idx),
          node_position(r_idx + 1, (phi_idx + 1) % num_angular_cells_),
          node_position(r_idx, (phi_idx + 1) % num_angular_cells_);
      std::unique_ptr<lf::geometry::Geometry> geom1 =
          std::make_unique<lf::geometry::TriaO1>(verts1);
      std::unique_ptr<lf::geometry::Geometry> geom2 =
          std::make_unique<lf::geometry::TriaO1>(verts2);
      mesh_factory_->AddEntity(ref_el, nonstd::span(trig1.data(), 3),
                               std::move(geom1));
      mesh_factory_->AddEntity(ref_el, nonstd::span(trig2.data(), 3),
                               std::move(geom2));
    }
  }

  return mesh_factory_->Build();
}

}  // end namespace projects::ipdg_stokes::mesh
