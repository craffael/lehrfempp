#include "mesh_builder.h"
#include "entity.h"
#include <lf/geometry/geometry.h>
#include <lf/geometry/point.h>
#include <lf/geometry/segment_o1.h>
#include <lf/geometry/tria_o1.h>
#include <lf/geometry/quad_o1.h>

namespace lf::mesh::hybrid2d {

  std::unique_ptr<mesh::Mesh>  TPTriagMeshBuilder::Build() {
    const size_type nx = no_of_x_cells_;
    const size_type ny = no_of_y_cells_;
    // Total number of entities in the mesh
    // Each triangle is split into two squares
    const int no_of_cells = 2*nx*ny;
    const int no_of_edges =  no_of_cells + (nx+1)*ny + nx*(ny+1);
    const int no_of_vertices = (nx+1)*(ny+1);
    // No mesh to build 
    if (no_of_cells == 0) return nullptr;
    // define rectangle; return if none
    const double x_size = top_right_corner_[0] - bottom_left_corner_[0];
    const double y_size = top_right_corner_[1] - bottom_left_corner_[1];
    if ((x_size <= 0.0) || (y_size <= 0.0)) return nullptr;
    // meshwidths
    const double hx = x_size/nx;
    const double hy = y_size/ny;

    // Initialize vertices
    std::vector<size_type> v_idx(no_of_vertices);
    int node_cnt = 0; // index of current vertex: lexikographic numbering
    for(int i=0;i<=nx;++i)
      for(int j=0;j<=ny;++j,++node_cnt) {
	// Tensor-product node locations
	coord_t node_coord(2); node_coord << i*hx,j*hy;
	// Create suitable geometry object
	v_idx[node_cnt] = AddPoint(node_coord);
      }

    // Initialize edges
    // Just in case of index permutations
    std::vector<size_type> e_idx(no_of_edges);
    int edge_cnt = 0; // index of currrent edge
    // First horizontal edges
    for(int i=0;i<nx;++i)
      for(int j=0;j<=ny;++j,++edge_cnt) {
	// Indices of the two endpoints of the edge
	const size_type first_endpoint_idx = v_idx[VertexIndex(i,j)];
	const size_type second_endpoint_idx = v_idx[VertexIndex(i+1,j)];
	lf::base::ForwardRange<const size_type>
	  nodes_index_list{first_endpoint_idx,second_endpoint_idx};
	// Coordinates of endpoints a columns of a 2x2 matrix
	Eigen::Matrix<double, Eigen::Dynamic, 2> edge_geo(2,2);
	edge_geo << i*hx,(i+1)*hx,j*hy,j*hy;
	e_idx[edge_cnt] = AddEntity(lf::base::RefEl::kSegment(),
				    nodes_index_list,
				    std::make_unique<geometry::SegmentO1>(edge_geo));
      }
    // Next vertical edges
    for(int i=0;i<=nx;++i)
      for(int j=0;j<ny;++j,++edge_cnt) {
	// Indices of the two endpoints of the edge
	const size_type first_endpoint_idx = v_idx[VertexIndex(i,j)];
	const size_type second_endpoint_idx = v_idx[VertexIndex(i,j+1)];
	lf::base::ForwardRange<const size_type>
	  nodes_index_list{first_endpoint_idx,second_endpoint_idx};
	// Coordinates of endpoints a columns of a 2x2 matrix
	Eigen::Matrix<double, Eigen::Dynamic, 2> edge_geo(2,2);
	edge_geo << i*hx,i*hx,j*hy,(j+1)*hy;
	e_idx[edge_cnt] = AddEntity(lf::base::RefEl::kSegment(),
				    nodes_index_list,
				    std::make_unique<geometry::SegmentO1>(edge_geo));
      }
    // Then the skew edges (diagonals of squares)
    for(int i=0;i<nx;++i)
      for(int j=0;j<ny;++j,++edge_cnt) {
	// Indices of the two endpoints of the edge
	const size_type first_endpoint_idx = v_idx[VertexIndex(i,j)];
	const size_type second_endpoint_idx = v_idx[VertexIndex(i+1,j+1)];
	lf::base::ForwardRange<const size_type>
	  nodes_index_list{first_endpoint_idx,second_endpoint_idx};
	// Coordinates of endpoints a columns of a 2x2 matrix
	Eigen::Matrix<double, Eigen::Dynamic, 2> edge_geo(2,2);
	edge_geo << i*hx,(i+1)*hx,j*hy,(j+1)*hy;
	e_idx[edge_cnt] = AddEntity(lf::base::RefEl::kSegment(),
				    nodes_index_list,
				    std::make_unique<geometry::SegmentO1>(edge_geo));
      }

    // Finally initialize the triangles
    // Index remapping for triangles
    std::vector<size_type> t_idx(no_of_cells);

    int tria_cnt = 0; // index of currrent triangle
    for(int i=0;i<nx;++i)
      for(int j=0;j<ny;++j,tria_cnt+=2) {
	// Triangle above the diagonal
	// Indices of the vertices
	lf::base::ForwardRange<const size_type>
	vertex_index_list_up { v_idx[VertexIndex(i,j)],
                            v_idx[VertexIndex(i+1,j+1)],
	                    v_idx[VertexIndex(i,j+1)] };
	// Construct geometry
	Eigen::Matrix<double, Eigen::Dynamic, 3> tria_geo_up(2,3);
	tria_geo_up << i*hx,(i+1)*hx,i*hx,j*hy,(j+1)*hy,(j+1)*hy;
	// Enroll the triangle entity
	t_idx[tria_cnt] =  AddEntity(lf::base::RefEl::kTria(),
				     vertex_index_list_up,
				     std::make_unique<geometry::TriaO1>(tria_geo_up));
	// Triangle below the diagonal
	// Indices of the vertices
	lf::base::ForwardRange<const size_type>
	  vertex_index_list_low { v_idx[VertexIndex(i,j)],
	                          v_idx[VertexIndex(i+1,j)],
	                          v_idx[VertexIndex(i+1,j+1)] };
	// Construct geometry
	Eigen::Matrix<double, Eigen::Dynamic, 3> tria_geo_low(2,3);
	tria_geo_low << i*hx,(i+1)*hx,(i+1)*hx,j*hy,j*hy,(j+1)*hy;
	auto tria_geo_low_ptr = std::make_unique<geometry::TriaO1>(tria_geo_low);
	// Generate the triangle entity
	t_idx[tria_cnt+1] =  AddEntity(lf::base::RefEl::kTria(),
				       vertex_index_list_low,
				       std::make_unique<geometry::TriaO1>(tria_geo_low));
      }
  } // end Build()
  
}
