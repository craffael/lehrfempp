#include "sphere_triag_mesh_builder.h"

#include <lf/base/base.h>
#include <lf/base/span.h>
#include <lf/geometry/geometry.h>

#include <array>
#include <cmath>
#include <vector>

namespace projects::hldo_sphere::mesh {

std::shared_ptr<lf::mesh::Mesh> SphereTriagMeshBuilder::Build() {
  const double r = radius_;
  const lf::base::size_type l = refinement_level_;

  // the order in which we create mesh elements will be top down
  // in rings, of which the number is computed based on the
  // refinement_level_ = l
  const lf::base::size_type ring_count = std::pow(2, l + 1) + 1;

  // number of overall vertices
  const lf::base::size_type vertex_count =
      2 + 4 * (ring_count / 2) +
      4 * ((ring_count - 2) / 2 + 1) * ((ring_count - 2) / 2);

  // difference in angle between the rings
  const double d_theta = M_PI / (ring_count - 1);

  // stores the mappings form the global vertex indices to the local ones on the
  // rings including the corresponding ring
  // The vector will be filled lazy, hence
  // USE `local_vert_idx` intstead of querying the vector
  std::vector<std::pair<lf::base::size_type, lf::base::size_type>>
      vertex_loc_glob_map_vector(vertex_count, std::make_pair(0, 0));

  // returns the ring index with respect to the middle.
  // This means the first and the last ring will have index 0 and the second and
  // second last index 1 ...
  // requires 0 <= `r_idx` < `ring_count`
  auto r_tilde = [&](lf::base::size_type r_idx) -> lf::base::size_type {
    LF_ASSERT_MSG(0 <= r_idx && r_idx < ring_count, "ring index out of range");
    if (r_idx <= (ring_count / 2)) {
      return r_idx;
    } else {
      return ring_count - r_idx - 1;
    }
  };

  // returns the number of vertices on a ring
  // requires 0 <= `r_idx` < `ring_count`
  auto vert_count_on_ring =
      [&](lf::base::size_type r_idx) -> lf::base::size_type {
    LF_ASSERT_MSG(0 <= r_idx && r_idx < ring_count, "ring index out of range");
    if (r_tilde(r_idx) == 0) {
      return 1;
    } else {
      return r_tilde(r_idx) * 4;
    }
  };

  // returns the ring index on which the vertex is located as well as
  // the local vertex index on that ring
  //
  // Internally computes the values only on the first access and then stores
  // them in a vector for later calls of the function
  //
  // requires 0 <= `r_idx` < `ring_count`
  auto local_vert_idx = [&](lf::base::size_type v_idx)
      -> std::pair<lf::base::size_type, lf::base::size_type> {
    LF_ASSERT_MSG(0 <= v_idx && v_idx < vertex_count,
                  "vertex index out of range");

    std::pair<lf::base::size_type, lf::base::size_type> zero_pair =
        std::make_pair(0, 0);

    if (v_idx == 0) {
      return zero_pair;
    }

    if (vertex_loc_glob_map_vector[v_idx] == zero_pair) {
      // subtract the one because it is not the first index
      lf::base::size_type v_idx_rest = v_idx - 1;
      lf::base::size_type r_idx = 1;
      while (vert_count_on_ring(r_idx) <= v_idx_rest) {
        v_idx_rest -= vert_count_on_ring(r_idx);
        r_idx++;
      }

      // fill the vector
      vertex_loc_glob_map_vector[v_idx] = std::make_pair(r_idx, v_idx_rest);
    }

    return vertex_loc_glob_map_vector[v_idx];
  };

  // returns the global index of a vertex given the ring_index on which the
  // vertex is located and the local vertex index on the ring
  //
  // requires 0 <= `r_idx` < `ring_count`
  // requires 0 <= `v_local_idx` < `vert_count_on_ring(r_idx)`
  auto global_vert_idx =
      [&](lf::base::size_type r_idx,
          lf::base::size_type v_local_idx) -> lf::base::size_type {
    LF_ASSERT_MSG(0 <= r_idx && r_idx < ring_count, "ring index out of range");
    LF_ASSERT_MSG(0 <= v_local_idx && v_local_idx < vert_count_on_ring(r_idx),
                  "vertex index out of range");

    // TODO implement with closed form formula
    lf::base::size_type v_idx = 0;
    for (lf::base::size_type r = 0; r < r_idx; ++r) {
      v_idx += vert_count_on_ring(r);
    }
    v_idx += v_local_idx;
    return v_idx;
  };

  // returns the position of a vertex in the space $\mathbb{R}^3$, given the
  // ring index and the vertex index on that ring
  //
  // requires 0 <= `r_idx` < `ring_count`
  // requires 0 <= `v_idx` < `vert_count_on_ring(r_idx)`
  auto node_position = [&](lf::base::size_type r_idx,
                           lf::base::size_type v_idx) -> Eigen::Vector3d {
    LF_ASSERT_MSG(0 <= r_idx && r_idx < ring_count, "ring index out of range");
    LF_ASSERT_MSG(0 <= v_idx && v_idx < vert_count_on_ring(r_idx),
                  "vertex index out of range");
    const double theta = r_idx * d_theta;
    const double d_phi = 2.0 * M_PI / vert_count_on_ring(r_idx);
    const double phi = v_idx * d_phi;
    Eigen::Vector3d pos;
    pos << r * std::cos(phi) * std::sin(theta),
        r * std::sin(phi) * std::sin(theta), r * std::cos(theta);
    return pos;
  };

  // returns the position of a vertex in the space $\mathbb{R}^3$, given the
  // global index of the vertex
  //
  // requires 0 <= `v_idx` < `vertex_count`
  auto global_node_position =
      [&](lf::base::size_type v_idx) -> Eigen::Vector3d {
    LF_ASSERT_MSG(0 <= v_idx && v_idx < vertex_count,
                  "vertex index out of range");
    lf::base::size_type r_idx;
    lf::base::size_type v_local_idx;
    std::tie(r_idx, v_local_idx) = local_vert_idx(v_idx);
    return node_position(r_idx, v_local_idx);
  };

  // returns the index \in \{0,1,2,3\} which quarter a given local index
  // is to be found.
  //
  // requires be 0 < `r_idx` < `ring_count`
  // requires 0 <= `v_idx` < `vert_count_on_ring(r_idx)`
  auto get_quarter = [&](lf::base::size_type r_idx,
                         lf::base::size_type v_idx) -> lf::base::size_type {
    LF_ASSERT_MSG(0 < r_idx && r_idx < ring_count, "ring index out of range");
    LF_ASSERT_MSG(0 <= v_idx && v_idx < vert_count_on_ring(r_idx),
                  "vertex index out of range");
    lf::base::size_type num_v_per_quarter = vert_count_on_ring(r_idx) / 4;
    return v_idx / num_v_per_quarter;
  };

  // Add all vertices to the mesh
  for (lf::base::size_type i = 0; i < vertex_count; ++i) {
    lf::base::size_type r_idx;
    lf::base::size_type v_local_idx;
    std::tie(r_idx, v_local_idx) = local_vert_idx(i);
    mesh_factory_->AddPoint(node_position(r_idx, v_local_idx));
  }

  // Construct the triangles
  // two for every vertex except the ones the first and last ring
  const lf::base::RefEl ref_el = lf::base::RefEl::kTria();
  for (lf::base::size_type i = 1; i < vertex_count - 1; ++i) {
    lf::base::size_type r_idx;
    lf::base::size_type v_local_idx;
    std::tie(r_idx, v_local_idx) = local_vert_idx(i);
    const int quarter = get_quarter(r_idx, v_local_idx);

    lf::base::size_type v1;
    lf::base::size_type v2;
    lf::base::size_type v3;

    // construct upper triangle
    if (vert_count_on_ring(r_idx) > vert_count_on_ring(r_idx - 1)) {
      v1 = i;
      v2 =
          global_vert_idx(r_idx, (v_local_idx + 1) % vert_count_on_ring(r_idx));
      v3 = global_vert_idx(
          r_idx - 1, (v_local_idx - quarter) % vert_count_on_ring(r_idx - 1));
    } else {
      v1 = i;
      v2 = global_vert_idx(r_idx - 1, v_local_idx + 1 + quarter);
      v3 =
          global_vert_idx(r_idx, (v_local_idx + 1) % vert_count_on_ring(r_idx));
    }

    const std::array<lf::base::size_type, 3> upper_trig = {v1, v2, v3};

    Eigen::Matrix<double, 3, 3> upper_verts;
    upper_verts << node_position(r_idx, v_local_idx), global_node_position(v2),
        global_node_position(v3);

    std::unique_ptr<lf::geometry::Geometry> upper_geom =
        std::make_unique<lf::geometry::TriaO1>(upper_verts);
    mesh_factory_->AddEntity(ref_el, nonstd::span(upper_trig.data(), 3),
                             std::move(upper_geom));

    // construct lower triangle
    if (vert_count_on_ring(r_idx) > vert_count_on_ring(r_idx + 1)) {
      v1 = i;
      v2 = global_vert_idx(
          r_idx + 1, (v_local_idx - quarter) % vert_count_on_ring(r_idx + 1));
      v3 =
          global_vert_idx(r_idx, (v_local_idx + 1) % vert_count_on_ring(r_idx));
    } else {
      v1 = i;
      v2 = global_vert_idx(r_idx + 1, v_local_idx + 1 + quarter);
      v3 =
          global_vert_idx(r_idx, (v_local_idx + 1) % vert_count_on_ring(r_idx));
    }

    const std::array<lf::base::size_type, 3> lower_trig = {v1, v2, v3};

    Eigen::Matrix<double, 3, 3> lower_verts;
    lower_verts << node_position(r_idx, v_local_idx), global_node_position(v2),
        global_node_position(v3);

    std::unique_ptr<lf::geometry::Geometry> lower_geom =
        std::make_unique<lf::geometry::TriaO1>(lower_verts);
    mesh_factory_->AddEntity(ref_el, nonstd::span(lower_trig.data(), 3),
                             std::move(lower_geom));
  }

  return mesh_factory_->Build();
}

}  // namespace projects::hldo_sphere::mesh
