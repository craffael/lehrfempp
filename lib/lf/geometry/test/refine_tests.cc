#include <gtest/gtest.h>
#include <lf/geometry/geometry.h>
#include <iostream>

namespace lf::geometry::test {

  std::ostream  &
  operator << (std::ostream &o,const std::vector<std::unique_ptr<Geometry>> &geo_uptrs) {
    const int n_uptrs =  geo_uptrs.size();
    o << "Array of size " << n_uptrs << " of Geometry pointers: " << std::endl;
    for (int l = 0 ; l < n_uptrs; l++) {
      if (geo_uptrs[l] != nullptr) {
	const lf::base::RefEl ref_el = geo_uptrs[l]->RefEl();
	o << "\t ===> Entity " << l << " = " << ref_el.ToString()
	  << " with node positions " << std::endl;
	const Eigen::MatrixXd ref_el_vertices(ref_el.NodeCoords());
	o << geo_uptrs[l]->Global(ref_el_vertices) << std::endl;
      }
      else o << "!!! Pointer " << l << " not valid!" << std::endl;
    }
    return o;
  }
  
  TEST(RefineTest,SegRef) {
    // Create edge geometry
    Eigen::Matrix<double, Eigen::Dynamic, 2> ed_coords(2,2);
    ed_coords << 1,0,2,1;
    auto edge_ptr = std::make_unique<SegmentO1>(ed_coords);

    // Refine
    std::vector<std::unique_ptr<Geometry>> 
      ed_copy(edge_ptr->ChildGeometry((int)RefinementPattern::rp_copy,-1,-1));
    std::vector<std::unique_ptr<Geometry>> 
      ed_split(edge_ptr->ChildGeometry((int)RefinementPattern::rp_split,-1,-1));
    EXPECT_EQ(ed_copy.size(),1);
    EXPECT_NE(ed_copy[0],nullptr);
    EXPECT_EQ(ed_split.size(),2);
    EXPECT_NE(ed_split[0],nullptr);
    EXPECT_NE(ed_split[1],nullptr);

    // Retrieve coordinates of endpoints of children
    const Eigen::MatrixXd ref_seg_endpoints(lf::base::RefEl::kSegment().NodeCoords());
    std::cout << "Geometry of parent = " << std::endl
	      << edge_ptr->Global(ref_seg_endpoints) << std::endl;
    std::cout << "Geometry of copy = " << std::endl
	      << ed_copy[0]->Global(ref_seg_endpoints) << std::endl;
    std::cout << "Geometry of child 0 = "  << std::endl
	      << ed_split[0]->Global(ref_seg_endpoints) << std::endl;
    std::cout << "Geometry of child 1 = "  << std::endl
	      << ed_split[1]->Global(ref_seg_endpoints) << std::endl;

    std::cout << ed_split;
  }

  TEST(RefineTest,TriaRef) {
    // Reference coordinate of triangle corners
    const Eigen::Matrix<double,2,3> ref_tria_corners(lf::base::RefEl::kTria().NodeCoords());
    // Create a triangle
    Eigen::Matrix<double,2,3> tria_coords;
    tria_coords.col(0) = Eigen::Vector2d({1,1});
    tria_coords.col(1) = Eigen::Vector2d({3,-1});
    tria_coords.col(2) = Eigen::Vector2d({3,3});
    std::unique_ptr<TriaO1> tria_geo_ptr = std::make_unique<TriaO1>(tria_coords);
    EXPECT_NE(tria_geo_ptr,nullptr);
    std::cout << "Parent triangle " << std::endl
     	      << tria_geo_ptr->Global(ref_tria_corners)
     	      << std::endl;
		 
    // Check the various refinements
    std::cout << "***** rp_copy(0) : child " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_copy,-1,0))
	      << std::endl;
    // Bisection refinement
    for (int anchor = 0; anchor < 3; anchor++) {
      std::cout << "*** rp_bisect(" << anchor << ") child "
		<< (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_bisect,anchor,-1))
		<< std::endl;
    }
    // Trisection refinement
    for (int anchor = 0; anchor < 3; anchor++) {
      std::cout << "*** rp_trisect(anchor = " << anchor << ") : children " << std::endl 
		<< (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect,anchor,-1))
		<< std::endl;
      std::cout << "*** rp_trisect_left(anchor = " << anchor << ") : children " << std::endl 
		<< (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_left,anchor,-1))
		<< std::endl;
    }
    // Splitting into four triangles
    for (int anchor = 0; anchor < 3; anchor++) {
      std::cout << "*** rp_quadsect(anchor = " << anchor << ") : children " << std::endl 
		<< (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect,anchor,-1))
		<< std::endl;
    }
    // Regular refinement
    std::cout << "*** rp_regular: children " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_regular,-1,-1))
	      << std::endl;
    // Bayrcentric refinement
    std::cout << "*** rp_baryccentric, children " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_barycentric,-1,-1))
	      << std::endl;
  }

  TEST(RefineTest,QuadRef) {
    // Reference coordinate of corners of quads and triangles
    const Eigen::MatrixXd ref_quad_corners(lf::base::RefEl::kQuad().NodeCoords());
    const Eigen::MatrixXd ref_tria_corners(lf::base::RefEl::kTria().NodeCoords());
    // Create a triangle
    Eigen::Matrix<double,2,4> quad_coords;
    quad_coords.col(0) = Eigen::Vector2d({1,1});
    quad_coords.col(1) = Eigen::Vector2d({5,-1});
    quad_coords.col(2) = Eigen::Vector2d({5,3});
    quad_coords.col(3) = Eigen::Vector2d({3,5});
    std::unique_ptr<QuadO1> quad_geo_ptr = std::make_unique<QuadO1>(quad_coords);
    EXPECT_NE(quad_geo_ptr,nullptr);
    std::cout << "Parent quadngle " << std::endl
     	      << quad_geo_ptr->Global(ref_quad_corners)
     	      << std::endl;
    
    // Check the various refinements
    std::cout << "rp_copy(0) : child " << std::endl 
	      << (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_copy,-1,0))
	      << std::endl;

    for (int anchor = 0; anchor < 4; anchor++) {
      std::cout << "** rp_trisect(anchor = " << anchor << ") : children " << std::endl 
		<< (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect,anchor,-1))
		<< std::endl;
    }
    for (int anchor = 0; anchor < 4; anchor++) {
      std::cout << "** rp_quadsect(anchor = " << anchor << ") : children " << std::endl 
		<< (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect,anchor,-1))
		<< std::endl;
    }
    for (int anchor = 0; anchor < 4; anchor++) {
      std::cout << "** rp_split(anchor = " << anchor << ") : children "
		<< -1 << " : " << std::endl 
		<< (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_split,anchor,-1))
		<< std::endl;
    }
   for (int anchor = 0; anchor < 4; anchor++) {
     std::cout << "*** rp_threeedge(anchor = " << anchor << ") : children " << std::endl
	       << (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_threeedge,anchor,-1))
	       << std::endl;
   }
   std::cout << "*** rp_regular, children " << std::endl 
	     << (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_regular,-1,-1))
	     << std::endl;
  }
}  // namespace lf::geometry::test
