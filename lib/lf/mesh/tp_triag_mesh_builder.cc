#include <lf/geometry/segment_o1.h>
#include <lf/geometry/tria_o1.h>

#include <iostream>
#include "mesh_interface.h"
#include "tp_triag_mesh_builder.h"

namespace lf::mesh::hybrid2d {

CONTROLDECLARECOMMENT(TPTriagMeshBuilder, output_ctrl_, "tpquad_output_ctrl",
                      "Diagnostics control for TPTriagMeshBuilder");

std::shared_ptr<mesh::Mesh> TPTriagMeshBuilder::Build() {
  using coord_t = Eigen::Vector2d;
  const size_type nx = no_of_x_cells_;
  const size_type ny = no_of_y_cells_;
  // Total number of entities in the mesh
  // Each triangle is split into two squares
  const int no_of_cells = 2 * nx * ny;
  const int no_of_edges = no_of_cells + (nx + 1) * ny + nx * (ny + 1);
  const int no_of_vertices = (nx + 1) * (ny + 1);
  // Diagnostics
  if (output_ctrl_) {
    std::cout << "TPmesh: " << no_of_cells << " cells, " << no_of_edges
              << " edges " << no_of_vertices << " vertices" << std::endl;
  }
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

  // Initialize vertices
  std::vector<size_type> v_idx(no_of_vertices);
  int node_cnt = 0;  // index of current vertex: lexikographic numbering
  for (size_type j = 0; j <= ny; ++j) {
    for (size_type i = 0; i <= nx; ++i, ++node_cnt) {
      // Tensor-product node locations
      coord_t node_coord(2);
      node_coord << i * hx, j * hy;
      // Diagnostics
      if (output_ctrl_) {
        std::cout << "Adding vertex " << node_cnt << ": " << node_coord
                  << std::endl;
      }
      // Enlist vertex
      v_idx[node_cnt] = mesh_factory_->AddPoint(node_coord);
    }
  }

  // Initialize edges
  // Just in case of index permutations
  std::vector<size_type> e_idx(no_of_edges);
  int edge_cnt = 0;  // index of currrent edge
  // First horizontal edges
  for (size_type i = 0; i < nx; ++i) {
    for (size_type j = 0; j <= ny; ++j, ++edge_cnt) {
      // Indices of the two endpoints of the edge
      auto first_endpoint_idx = v_idx[VertexIndex(i, j)];
      auto second_endpoint_idx = v_idx[VertexIndex(i + 1, j)];
      // Diagnostics
      if (output_ctrl_) {
        std::cout << "horizontal edge " << edge_cnt << ": "
                  << first_endpoint_idx << " <-> " << second_endpoint_idx
                  << std::endl;
      }
      lf::base::ForwardRange<const size_type> nodes_index_list{
          first_endpoint_idx, second_endpoint_idx};
      // Coordinates of endpoints a columns of a 2x2 matrix
      Eigen::Matrix<double, Eigen::Dynamic, 2> edge_geo(2, 2);
      edge_geo << i * hx, (i + 1) * hx, j * hy, j * hy;
      e_idx[edge_cnt] = mesh_factory_->AddEntity(
          lf::base::RefEl::kSegment(), nodes_index_list,
          std::make_unique<geometry::SegmentO1>(edge_geo));
    }
  }
  // Next vertical edges
  for (size_type i = 0; i <= nx; ++i) {
    for (size_type j = 0; j < ny; ++j, ++edge_cnt) {
      // Indices of the two endpoints of the edge
      const size_type first_endpoint_idx = v_idx[VertexIndex(i, j)];
      const size_type second_endpoint_idx = v_idx[VertexIndex(i, j + 1)];
      // Diagnostics
      if (output_ctrl_) {
        std::cout << "vertical edge " << edge_cnt << ": " << first_endpoint_idx
                  << " <-> " << second_endpoint_idx << std::endl;
      }
      lf::base::ForwardRange<const size_type> nodes_index_list{
          first_endpoint_idx, second_endpoint_idx};
      // Coordinates of endpoints a columns of a 2x2 matrix
      Eigen::Matrix<double, Eigen::Dynamic, 2> edge_geo(2, 2);
      edge_geo << i * hx, i * hx, j * hy, (j + 1) * hy;
      e_idx[edge_cnt] = mesh_factory_->AddEntity(
          lf::base::RefEl::kSegment(), nodes_index_list,
          std::make_unique<geometry::SegmentO1>(edge_geo));
    }
  }
  // Then the skew edges (diagonals of squares)
  for (size_type i = 0; i < nx; ++i) {
    for (size_type j = 0; j < ny; ++j, ++edge_cnt) {
      // Indices of the two endpoints of the edge
      const size_type first_endpoint_idx = v_idx[VertexIndex(i, j)];
      const size_type second_endpoint_idx = v_idx[VertexIndex(i + 1, j + 1)];
      // Diagnostics
      if (output_ctrl_) {
        std::cout << "diagonal edge " << edge_cnt << ": " << first_endpoint_idx
                  << " <-> " << second_endpoint_idx << std::endl;
      }
      lf::base::ForwardRange<const size_type> nodes_index_list{
          first_endpoint_idx, second_endpoint_idx};
      // Coordinates of endpoints a columns of a 2x2 matrix
      Eigen::Matrix<double, Eigen::Dynamic, 2> edge_geo(2, 2);
      edge_geo << i * hx, (i + 1) * hx, j * hy, (j + 1) * hy;
      e_idx[edge_cnt] = mesh_factory_->AddEntity(
          lf::base::RefEl::kSegment(), nodes_index_list,
          std::make_unique<geometry::SegmentO1>(edge_geo));
    }
  }

  // Finally initialize the triangles
  // Index remapping for triangles
  std::vector<size_type> t_idx(no_of_cells);

  size_type tria_cnt = 0;  // index of currrent triangle
  for (size_type i = 0; i < nx; ++i) {
    for (size_type j = 0; j < ny; ++j, tria_cnt += 2) {
      // Triangle above the diagonal
      // Indices of the vertices
      lf::base::ForwardRange<const size_type> vertex_index_list_up{
          v_idx[VertexIndex(i, j)], v_idx[VertexIndex(i + 1, j + 1)],
          v_idx[VertexIndex(i, j + 1)]};
      // Construct geometry
      Eigen::Matrix<double, Eigen::Dynamic, 3> tria_geo_up(2, 3);
      tria_geo_up << i * hx, (i + 1) * hx, i * hx, j * hy, (j + 1) * hy,
          (j + 1) * hy;
      // Enroll the triangle entity
      t_idx[tria_cnt] = mesh_factory_->AddEntity(
          lf::base::RefEl::kTria(), vertex_index_list_up,
          std::make_unique<geometry::TriaO1>(tria_geo_up));
      // Triangle below the diagonal
      // Indices of the vertices
      lf::base::ForwardRange<const size_type> vertex_index_list_low{
          v_idx[VertexIndex(i, j)], v_idx[VertexIndex(i + 1, j)],
          v_idx[VertexIndex(i + 1, j + 1)]};
      // Construct geometry
      Eigen::Matrix<double, Eigen::Dynamic, 3> tria_geo_low(2, 3);
      tria_geo_low << i * hx, (i + 1) * hx, (i + 1) * hx, j * hy, j * hy,
          (j + 1) * hy;
      auto tria_geo_low_ptr = std::make_unique<geometry::TriaO1>(tria_geo_low);
      // Generate the triangle entity
      t_idx[tria_cnt + 1] = mesh_factory_->AddEntity(
          lf::base::RefEl::kTria(), vertex_index_list_low,
          std::make_unique<geometry::TriaO1>(tria_geo_low));
    }
  }
  return mesh_factory_->Build();
}  // end Build()

CONTROLDECLARECOMMENT(TPQuadMeshBuilder, output_ctrl_, "tpquad_output_ctrl",
                      "Diagnostics control for TPQuadMeshBuilder");

std::shared_ptr<mesh::Mesh> TPQuadMeshBuilder::Build() {
  using coord_t = Eigen::Vector2d;
  const size_type nx = no_of_x_cells_;
  const size_type ny = no_of_y_cells_;
  // Total number of entities in the mesh
  // Each triangle is split into two squares
  const int no_of_cells = nx * ny;  // no of rectangles
  const int no_of_edges = (nx + 1) * ny + nx * (ny + 1);
  const int no_of_vertices = (nx + 1) * (ny + 1);
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

  if (output_ctrl_ > 0) {
    std::cout << "TPQuadmesh: " << no_of_cells << " cells, " << no_of_edges
              << " edges " << no_of_vertices << " vertices, "
              << "meshwidths hx/hy = " << hx << "/" << hy << std::endl;
  }

  // Initialize vertices
  std::vector<size_type> v_idx(no_of_vertices);
  int node_cnt = 0;  // index of current vertex: lexikographic numbering
  for (size_type j = 0; j <= ny; ++j) {
    for (size_type i = 0; i <= nx; ++i, ++node_cnt) {
      // Tensor-product node locations
      coord_t node_coord(2);
      node_coord << i * hx, j * hy;
      // Diagnostics
      if (output_ctrl_ > 0) {
        std::cout << "Adding vertex " << node_cnt << ": "
                  << node_coord.transpose() << std::endl;
      }
      // Register vertex
      v_idx[node_cnt] = mesh_factory_->AddPoint(node_coord);
    }
  }
  // Set quadrilateral cells
  // Index remapping for cells
  std::vector<size_type> t_idx(no_of_cells);

  size_type quad_cnt = 0;  // index of currrent triangle
  for (size_type i = 0; i < nx; ++i) {
    for (size_type j = 0; j < ny; ++j, quad_cnt++) {
      // Indices of vertices of quadrilaterial (i,j)
      lf::base::ForwardRange<const size_type> vertex_index_list{
          v_idx[VertexIndex(i, j)], v_idx[VertexIndex(i + 1, j)],
          v_idx[VertexIndex(i + 1, j + 1)], v_idx[VertexIndex(i, j + 1)]};
      // Construct geometry of rectangle
      Eigen::Matrix<double, 2, 4> quad_geo(2, 4);
      quad_geo << i * hx, (i + 1) * hx, (i + 1) * hx, i * hx, j * hy, j * hy,
          (j + 1) * hy, (j + 1) * hy;
      // Diagnostics
      if (output_ctrl_ > 0) {
        std::cout << "Adding quad " << quad_cnt << ": " << quad_geo
                  << std::endl;
      }

      // Request production of the cell from the MeshFactory
      // Straight edges will be created as needed
      t_idx[quad_cnt] = mesh_factory_->AddEntity(
          lf::base::RefEl::kQuad(), vertex_index_list,
          std::make_unique<geometry::QuadO1>(quad_geo));
    }
  }
  return mesh_factory_->Build();
}

}  // namespace lf::mesh::hybrid2d
