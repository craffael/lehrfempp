#include <gtest/gtest.h>
#include <lf/geometry/geometry.h>
#include <iostream>

namespace lf::geometry::test {

  TEST(RefineTest,SegRef) {
    // Create edge geometry
    Eigen::Matrix<double, Eigen::Dynamic, 2> ed_coords(2,2);
    ed_coords << 1,0,2,1;
    auto edge_ptr = std::make_unique<SegmentO1>(ed_coords);

    // Refine
    
    auto ed_copy = edge_ptr->ChildGeometry((int)RefinementPattern::rp_copy,-1,0);
    auto child0 = edge_ptr->ChildGeometry((int)RefinementPattern::rp_split,-1,0);
    auto child1 = edge_ptr->ChildGeometry((int)RefinementPattern::rp_split,-1,1);
    EXPECT_NE(child0,nullptr);
    EXPECT_NE(child1,nullptr);
    EXPECT_NE(ed_copy,nullptr);

    // Retrieve coordinates of endpoints of children
    const Eigen::MatrixXd ref_seg_endpoints(lf::base::RefEl::kSegment().NodeCoords());
    std::cout << "Geometry of parent = " << std::endl
	      << edge_ptr->Global(ref_seg_endpoints) << std::endl;
    std::cout << "Geometry of copy = " << std::endl
	      << ed_copy->Global(ref_seg_endpoints) << std::endl;
    std::cout << "Geometry of child 0 = "  << std::endl
	      << child0->Global(ref_seg_endpoints) << std::endl;
    std::cout << "Geometry of child 1 = "  << std::endl
	      << child1->Global(ref_seg_endpoints) << std::endl;
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
    std::cout << "rp_copy(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_copy,-1,0))
      ->Global(ref_tria_corners) << std::endl;
    // Bisection refinement
    for (int anchor = 0; anchor < 3; anchor++) {
      for (int selector = 0; selector < 2; selector++) {
	std::cout << "rp_bisect(" << anchor << ") child "
		  << selector << " : " << std::endl 
		  << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_bisect,anchor,selector))->Global(ref_tria_corners) << std::endl;
      }
    }
    // Trisection refinement
    for (int anchor = 0; anchor < 3; anchor++) {
      for (int selector = 0; selector < 3; selector++) {
	std::cout << "rp_trisect(" << anchor << ") child "
		  << selector << " : " << std::endl 
		  << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect,anchor,selector))->Global(ref_tria_corners) << std::endl;
	std::cout << "rp_trisect_left(" << anchor << ") child "
		  << selector << " : " << std::endl 
		  << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_left,anchor,selector))->Global(ref_tria_corners) << std::endl;
      }
    }
    // Splitting into four triangles
    for (int anchor = 0; anchor < 3; anchor++) {
      for (int selector = 0; selector < 4; selector++) {
	std::cout << "rp_quadsect(" << anchor << ") child "
		  << selector << " : " << std::endl 
		  << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect,anchor,selector))->Global(ref_tria_corners) << std::endl;
      }
    }
    // Regular refinement
    for (int selector = 0; selector < 4; selector++) {
      std::cout << "rp_regular, child "
		<< selector << " : " << std::endl 
		<< (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_regular,-1,selector))->Global(ref_tria_corners) << std::endl;
    }
    // Bayrcentric refinement
    for (int selector = 0; selector < 6; selector++) {
      std::cout << "rp_baryccentric, child "
		<< selector << " : " << std::endl 
		<< (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_barycentric,-1,selector))->Global(ref_tria_corners) << std::endl;
    }
    
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
    std::cout << "rp_copy(0) : " << std::endl 
	      << (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_copy,-1,0))
      ->Global(ref_quad_corners) << std::endl;

    for (int anchor = 0; anchor < 4; anchor++) {
      for (int selector = 0; selector < 3; selector++) {
	std::cout << "rp_trisect(" << anchor << ") child "
		  << selector << " : " << std::endl 
		  << (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect,anchor,selector))->Global(ref_tria_corners) << std::endl;
      }
    }
    for (int anchor = 0; anchor < 4; anchor++) {
      for (int selector = 0; selector < 4; selector++) {
	std::cout << "rp_quadsect(" << anchor << ") child "
		  << selector << " : " << std::endl 
		  << (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect,anchor,selector))->Global(ref_tria_corners) << std::endl;
      }
    }
    for (int anchor = 0; anchor < 4; anchor++) {
      for (int selector = 0; selector < 2; selector++) {
	std::cout << "rp_split(" << anchor << ") child "
		  << selector << " : " << std::endl 
		  << (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_split,anchor,selector))->Global(ref_quad_corners) << std::endl;
      }
    }
    for (int anchor = 0; anchor < 4; anchor++) {
      for (int selector = 0; selector < 4; selector++) {
	std::cout << "rp_threeedge(" << anchor << ") child "
		  << selector << " : " << std::endl 
		  << (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_threeedge,anchor,selector))->Global((selector==0)?ref_quad_corners:ref_tria_corners) << std::endl;
      }
    }
    for (int selector = 0; selector < 4; selector++) {
      std::cout << "rp_regular, child "
		<< selector << " : " << std::endl 
		<< (quad_geo_ptr->ChildGeometry((int)RefinementPattern::rp_regular,-1,selector))->Global(ref_quad_corners) << std::endl;
    }
  }
}  // namespace lf::geometry::test
