/**
 * @file
 * @brief Defines the Build method of the TorusMeshBuilder
 * @author Anian Ruoss
 * @date   2018-10-22 16:06:17
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES

#include <lf/geometry/geometry.h>
#include <cmath>

#include <iostream>
#include "lf/mesh/mesh_interface.h"
#include "torus_mesh_builder.h"

#include <spdlog/spdlog.h>

namespace lf::mesh::utils {

std::shared_ptr<spdlog::logger>& TorusMeshBuilder::Logger() {
  static auto logger =
      base::InitLogger("lf::mesh::utils::TorusMeshBuilder::Logger");
  return logger;
}

std::shared_ptr<mesh::Mesh> TorusMeshBuilder::Build() {
  using coord_t = Eigen::Vector3d;

  const size_type nx = num_of_x_cells_;
  const size_type ny = num_of_y_cells_;

  // opposite edges and vertices are identical and only counted once
  const unsigned no_of_cells = nx * ny;
  const unsigned no_of_edges = 2 * nx * ny;
  const unsigned no_of_vertices = nx * ny;

  if (no_of_cells == 0) {
    return nullptr;
  }

  const double x_size = top_right_corner_[0] - bottom_left_corner_[0];
  const double y_size = top_right_corner_[1] - bottom_left_corner_[1];

  if ((x_size <= 0.0) || (y_size <= 0.0)) {
    return nullptr;
  }

  // calculate mesh width of uniform grid on rectangle
  const double hx = x_size / nx;
  const double hy = y_size / ny;

  // parametrize torus: https://en.wikipedia.org/wiki/Torus#Geometry
  const double r = x_size / (2. * M_PI);
  const double R = y_size / (2. * M_PI);
  auto theta = [r, hx](double i) -> double { return (i * hx) / r; };
  auto phi = [R, hy](double j) -> double { return (j * hy) / R; };

  SPDLOG_LOGGER_DEBUG(
      Logger(),
      "TorusMesh: {} cells, {} edges, {} vertices, mesh widths hx/hy = {}/{}",
      no_of_cells, no_of_edges, no_of_vertices, hx, hy);

  // compute vertices of mesh on torus with lexicographic numbering
  std::vector<size_type> v_idx(no_of_vertices);
  int node_cnt = 0;

  for (size_type j = 0; j < ny; ++j) {
    for (size_type i = 0; i < nx; ++i, ++node_cnt) {
      // compute vertex coordinates
      coord_t node_coord;
      node_coord << (R + r * std::cos(theta(i))) * std::cos(phi(j)),
          (R + r * std::cos(theta(i))) * std::sin(phi(j)),
          r * std::sin(theta(i));

      SPDLOG_LOGGER_TRACE(Logger(), "Adding vertex {}: {}", node_cnt,
                          node_coord.transpose());
      // register vertex
      v_idx[node_cnt] = mesh_factory_->AddPoint(node_coord);
    }
  }

  // compute cells of mesh on torus
  std::vector<size_type> t_idx(no_of_cells);
  size_type quad_cnt = 0;

  for (size_type i = 0; i < nx; ++i) {
    for (size_type j = 0; j < ny; ++j, quad_cnt++) {
      // gather vertex indices of given cell wrapping around edges of rectangle
      std::vector<size_type> vertex_index_list{
          v_idx[VertexIndex(i, j)], v_idx[VertexIndex((i + 1) % nx, j)],
          v_idx[VertexIndex((i + 1) % nx, (j + 1) % ny)],
          v_idx[VertexIndex(i, (j + 1) % ny)]};

      Eigen::Matrix<double, 3, 4> quad_geo;
      quad_geo << (R + r * std::cos(theta(i))) * std::cos(phi(j)),
          (R + r * std::cos(theta(i + 1))) * std::cos(phi(j)),
          (R + r * std::cos(theta(i + 1))) * std::cos(phi(j + 1)),
          (R + r * std::cos(theta(i))) * std::cos(phi(j + 1)),
          (R + r * std::cos(theta(i))) * std::sin(phi(j)),
          (R + r * std::cos(theta(i + 1))) * std::sin(phi(j)),
          (R + r * std::cos(theta(i + 1))) * std::sin(phi(j + 1)),
          (R + r * std::cos(theta(i))) * std::sin(phi(j + 1)),
          r * std::sin(theta(i)), r * std::sin(theta(i + 1)),
          r * std::sin(theta(i + 1)), r * std::sin(theta(i));

      SPDLOG_LOGGER_TRACE(Logger(), "Adding quad {}:\n{}", quad_cnt, quad_geo);

      // request cell production from MeshFactory (straight edges will be
      // created as needed)
      t_idx[quad_cnt] = mesh_factory_->AddEntity(
          lf::base::RefEl::kQuad(), vertex_index_list,
          std::make_unique<geometry::QuadO1>(quad_geo));
    }
  }

  return mesh_factory_->Build();
}  // TorusMeshBuilder::Build()

}  // namespace lf::mesh::utils
