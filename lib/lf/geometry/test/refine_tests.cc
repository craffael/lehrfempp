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
    
    auto ed_copy = edge_ptr->ChildGeometry((int)RefinementPattern::rp_copy,0);
    auto child0 = edge_ptr->ChildGeometry((int)RefinementPattern::rp_split,0);
    auto child1 = edge_ptr->ChildGeometry((int)RefinementPattern::rp_split,1);
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
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_copy,0))
      ->Global(ref_tria_corners) << std::endl;
    std::cout << "rp_bisect_0(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_bisect_0,0))
      ->Global(ref_tria_corners) << std::endl;
    std::cout << "rp_bisect_0(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_bisect_0,1))
      ->Global(ref_tria_corners) << std::endl;
    std::cout << "rp_bisect_1(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_bisect_1,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_bisect_1(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_bisect_1,1))
      ->Global(ref_tria_corners) << std::endl;
    std::cout << "rp_bisect_2(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_bisect_2,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_bisect_2(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_bisect_2,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_01(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_01,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_01(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_01,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_01(2) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_01,2))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_02(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_02,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_02(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_02,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_02(2) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_02,2))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_10(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_10,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_10(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_10,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_10(2) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_10,2))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_12(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_12,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_12(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_12,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_12(2) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_12,2))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_20(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_20,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_20(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_20,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_20(2) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_20,2))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_21(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_21,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_21(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_21,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_trisect_21(2) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_trisect_21,2))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_0(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_0,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_0(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_0,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_0(2) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_0,2))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_0(3) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_0,3))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_1(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_1,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_1(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_1,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_1(2) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_1,2))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_1(3) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_1,3))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_2(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_2,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_2(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_2,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_2(2) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_2,2))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_quadsect_2(3) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_quadsect_2,3))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_regular(0) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_regular,0))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_regular(1) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_regular,1))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_regular(2) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_regular,2))
      ->Global(ref_tria_corners) << std::endl;
   std::cout << "rp_regular(3) : " << std::endl 
	      << (tria_geo_ptr->ChildGeometry((int)RefinementPattern::rp_regular,3))
      ->Global(ref_tria_corners) << std::endl;
  }
  
}  // namespace lf::geometry::test
