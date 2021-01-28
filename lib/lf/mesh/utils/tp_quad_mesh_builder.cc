#include "tp_quad_mesh_builder.h"

#include <lf/geometry/geometry.h>
#include <spdlog/spdlog.h>

#include <iostream>

#include "lf/mesh/mesh_interface.h"

namespace lf::mesh::utils {

std::shared_ptr<spdlog::logger>& TPQuadMeshBuilder::Logger() {
  static auto logger =
      base::InitLogger("lf::mesh::utils::TPQuadMeshBuilder::Logger");
  return logger;
}

std::shared_ptr<mesh::Mesh> TPQuadMeshBuilder::Build() {
  using coord_t = Eigen::Vector2d;
  const size_type nx = num_of_x_cells_;
  const size_type ny = num_of_y_cells_;
  // Total number of entities in the mesh
  // Each triangle is split into two squares
  const unsigned no_of_cells = nx * ny;  // no of rectangles
  const unsigned no_of_edges = (nx + 1) * ny + nx * (ny + 1);
  const unsigned no_of_vertices = (nx + 1) * (ny + 1);
  // No mesh to build
  if (no_of_cells == 0) {
    return nullptr;
  }
  // define rectangle; return if none
  const double x_size = top_right_corner_[0] - bottom_left_corner_[0];
  const double y_size = top_right_corner_[1] - bottom_left_corner_[1];
  if ((x_size <= 0.0) || (y_size <= 0.0)) {
    return nullptr;
  }
  // meshwidths
  const double hx = x_size / nx;
  const double hy = y_size / ny;

  SPDLOG_LOGGER_DEBUG(
      Logger(),
      "TPQuadmesh: {} cells, {} edges, {} vertices, meshwidths hx/hy = {}/{}",
      no_of_cells, no_of_edges, no_of_vertices, hx, hy);

  // Initialize vertices
  std::vector<size_type> v_idx(no_of_vertices);
  int node_cnt = 0;  // index of current vertex: lexicographic numbering
  for (size_type j = 0; j <= ny; ++j) {
    for (size_type i = 0; i <= nx; ++i, ++node_cnt) {
      // Tensor-product node locations
      coord_t node_coord(2);
      node_coord << bottom_left_corner_[0] + i * hx,
          bottom_left_corner_[1] + j * hy;
      // Diagnostics
      SPDLOG_LOGGER_TRACE(Logger(), "Adding vertex {}: {}", node_cnt,
                          node_coord.transpose());

      // Register vertex
      v_idx[node_cnt] = mesh_factory_->AddPoint(node_coord);
    }
  }
  // Set quadrilateral cells
  // Index remapping for cells
  std::vector<size_type> t_idx(no_of_cells);

  size_type quad_cnt = 0;  // index of current triangle
  for (size_type i = 0; i < nx; ++i) {
    for (size_type j = 0; j < ny; ++j, quad_cnt++) {
      // Indices of vertices of quadrilateral (i,j)
      std::vector<size_type> vertex_index_list{
          v_idx[VertexIndex(i, j)], v_idx[VertexIndex(i + 1, j)],
          v_idx[VertexIndex(i + 1, j + 1)], v_idx[VertexIndex(i, j + 1)]};
      // Construct geometry of rectangle
      Eigen::Matrix<double, 2, 4> quad_geo(2, 4);
      quad_geo << bottom_left_corner_[0] + i * hx,
          bottom_left_corner_[0] + (i + 1) * hx,
          bottom_left_corner_[0] + (i + 1) * hx,
          bottom_left_corner_[0] + i * hx, bottom_left_corner_[1] + j * hy,
          bottom_left_corner_[1] + j * hy,
          bottom_left_corner_[1] + (j + 1) * hy,
          bottom_left_corner_[1] + (j + 1) * hy;
      // Diagnostics
      SPDLOG_LOGGER_TRACE(Logger(), "Adding quad {}:\n{}", quad_cnt, quad_geo);

      // Request production of the cell from the MeshFactory
      // Straight edges will be created as needed
      t_idx[quad_cnt] = mesh_factory_->AddEntity(
          lf::base::RefEl::kQuad(), vertex_index_list,
          std::make_unique<geometry::QuadO1>(quad_geo));
    }
  }
  return mesh_factory_->Build();
}

}  // namespace lf::mesh::utils
