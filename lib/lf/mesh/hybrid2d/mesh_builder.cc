#include "mesh_builder.h"
#include "entity.h"
#include <lf/geometry/geometry.h>
#include <lf/geeometry/point.h>
#include <lf/geeometry/segment_o1.h>
#include <lf/geeometry/tria_o1.h>
#include <lf/geeometry/quad_o1.h>

namespace lf::mesh::hybrid2d {

  /** @brief Implementation of complex methods of hybrid2D::Meshbuilder
   * @{
   */
  MeshBuilder::size_type MeshBuilder::AddPoint(coord_t coord) {
    LF_ASSERT_MSG(!built_, "Build() already called.");
    LF_ASSERT_MSG(coord.rows() == dim_world_,
		  "coord has incompatible number of rows.");
    nodes_.emplace_back(std::move(coord));
    return nodes_.size() - 1;
  }

  MeshBuilder::size_type
  MeshBuilder::AddElement(const base::ForwardRange<const size_type>& nodes,
			  std::unique_ptr<geometry::Geometry>&& geometry) {
    LF_ASSERT_MSG(!built_, "Build() already called.");
    LF_ASSERT_MSG(geometry->DimGlobal() == dim_world_,
		  "geometry->DimGlobal() != dim_world_");
    LF_ASSERT_MSG(geometry->DimLocal() == 2, "geometry->DimLocal() != 2");

    std::vector<size_type> ns;
    ns.reserve(4);
    for (auto& n : nodes) {
      LF_ASSERT_MSG(n < nodes_.size(), "node " << n <<
      " specified in call to AddElement must be inserted with AddNode() first."
		    );
      ns.push_back(n);
    }
    LF_ASSERT_MSG(geometry->RefEl().NumNodes() == ns.size(),
		  "mismatch between number of nodes and RefEl of geometry.");

    elements_.emplace_back(std::move(ns), std::move(geometry));
    return elements_.size() - 1;
  }

  std::unique_ptr<mesh::Mesh> MeshBuilder::Build() {
    built_ = true;
    return std::make_unique<Mesh>(dim_world_, std::move(nodes_),
				  std::move(elements_));
  }
  /** @} */

  /** 
   * @brief Generate 2D triangular tensor product mesh 
   *
   */
  std::unique_ptr<mesh::Mesh>  TPTriagMeshBuilder::Build() {
    // Total number of entities in the mesh
    // Each triangle is split into two squares
    const int no_of_cells = 2*no_of_x_cells*no_of_y_cells;
    const int no_of_edges =  (no_of_cells + (no_of_x_cells_+1)*no_of_y_cells +
			      no_of_x_cells_*(no_of_y_cells_+1));
    const int no_of_vertices = (no_of_x_cells_+1)*(no_of_y_cells_+1);
    // No mesh to build 
    if (no_of_cells == 0) return nullptr;
    // define rectangle; return if none
    const double x_size = top_right_corner_[0] - bottom_left_corner[0];
    const double y_size = top_right_corner_[1] - bottom_left_corner[1];
    if ((xsize <= 0.0) || (ysize <= 0.0)) return nullptr;
    // meshwidths
    const double hx = x_size/no_of_x_cells_;
    const double hx = y_size/no_of_y_cells_;

    // Instantiate empty mesh object
    auto mesh = std::make_unique<Mesh>(2);
    
    // Sizes of arrays for mesh entities
    mesh->entities0_.reserve(no_of_vertices);
    mesh->entities1_.reserve(no_of_edges);
    mesh->entities2_.reserve(no_of_cells);
    

    // Initialize vertices
    // A vertex does not have any sub-entities; dummy argument
    std::array<std::vector<size_type>, 0> sub_vertex;
    int node_cnt = 0; // index of current vertex: lexikographic numbering
    for(int i=0;i<=no_of_x_cells_;++i)
      for(int j=0;j<=no_of_y_cells_;++j,++node_cnt) {
	// Tensor-product node locations
	const Eigen::Vector2d node_coord{i*hx,j*hy};
	// Create suitable geometry object
	auto point_ptr = std::make_unique<geometry::Point>(node_coord);
	mesh->entities0_.push_back(mesh,node_cnt,point_ptr,sub_vertex);
      }
    // Initialize edges
    // An edge has its two endpoints as sub-entities
    std::array<std::vector<size_type>, 1> sub_edge;
    int edge_cnt = 0; // index of currrent edge
    // First horizontal edges
    for(int i=0;i<no_of_x_cells_;++i)
      for(int j=0;j<=no_of_y_cells_;++j,++edge_cnt) {
	// Indices of the two endpoints of the edge
	const int first_endpoint_idx = VertexIndex(i,j);
	const int second_endpoint_idx = VertexIndex(i+1,j);
	sub_edge[0] = std::vector<size_type>{first_endpoint_idx,second_endpoint_idx};
	// Coordinates of endpoints a columns of a 2x2 matrix
	Eigen::Matrix<double, Eigen::Dynamic, 2>
	  edge_geo << i*hx,(i+1)*hx,j*hy,j*hy;
	auto seg_ptr = std::make_unique<geometry::SegmentO1>(edge_geo);
	mesh->entities1_.push_back(mesh,edge_cnt,sub_edge);
      }
    // Next vertical edges
    for(int i=0;i<=no_of_x_cells_;++i)
      for(int j=0;j<no_of_y_cells_;++j,++edge_cnt) {
	// Indices of the two endpoints of the edge
	const int first_endpoint_idx = VertexIndex(i,j);
	const int second_endpoint_idx = VertexIndex(i,j+1);
	sub_edge[0] = std::vector<size_type>{first_endpoint_idx,second_endpoint_idx};
	// Coordinates of endpoints a columns of a 2x2 matrix
	Eigen::Matrix<double, Eigen::Dynamic, 2>
	  edge_geo << i*hx,i*hx,j*hy,(j+1)*hy;
	auto seg_ptr = std::make_unique<geometry::SegmentO1>(edge_geo);
	mesh->entities1_.push_back(mesh,edge_cnt,sub_edge);
      }
    // Then the skew edges (diagonals of squares)
    for(int i=0;i<no_of_x_cells_;++i)
      for(int j=0;j<no_of_y_cells_;++j,++edge_cnt) {
	// Indices of the two endpoints of the edge
	const int first_endpoint_idx = VertexIndex(i,j);
	const int second_endpoint_idx = VertexIndex(i+1,j+1);
	sub_edge[0] = std::vector<size_type>{first_endpoint_idx,second_endpoint_idx};
	// Coordinates of endpoints a columns of a 2x2 matrix
	Eigen::Matrix<double, Eigen::Dynamic, 2>
	  edge_geo << i*hx,(i+1)*hx,j*hy,(j+1)*hy;
	auto seg_ptr = std::make_unique<geometry::SegmentO1>(edge_geo);
	mesh->entities1_.push_back(mesh,edge_cnt,sub_edge);
      }
    
    
    // Finally initialize the triangles
    // A triangle expects two arrays of sub-entity indices
    std::array<std::vector<size_type>, 2> sub_tria;
    int tria_cnt = 0; // index of currrent triangle
    for(int i=0;i<no_of_x_cells_;++i)
      for(int j=0;j<no_of_y_cells_;++j,++tria_cnt) {
	// Indices of the vertices
	sub_tria[1] = std::vector<size_type>
	  {VertexIndex(i,j),VertexIndex(i+1,j+1),VertexIndex(i,j+1)};
	// Indices of the adjacent edges
	sub_tria[0] = std:vector<size_type>{/* TODO */};
      }
  }
  
}
