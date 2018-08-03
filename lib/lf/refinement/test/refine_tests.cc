/** @file geo_ref_test.cc
 * tests for the topological refinement of reference entities
 */
#include <gtest/gtest.h>
#include <lf/refinement/refinement.h>
#include <stdlib.h>
#include <iostream>

namespace lf::refinement::test {

using lt_polys_t =
    std::vector<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>;

// Output of lattice polygons
void PrintLatticePolygons(const lt_polys_t &polygons);

void PrintLatticePolygons(std::ostream &o, const lt_polys_t &polygons) {
  int n_poly = polygons.size();
  if (n_poly > 0) {
    int dim = polygons[0].rows();
    if (dim > 0) {
      o << n_poly << " polygons in dimension " << dim << std::endl;
      for (int l = 0; l < n_poly; l++) {
        LF_VERIFY_MSG(polygons[l].rows() == dim,
                      "Polygon " << l << " dimension mismatch: " << dim
                                 << " <-> " << polygons[l].rows());
        o << "\t polygon " << l << ": " << std::endl
          << polygons[l] << std::endl;
      }
    } else
      o << "Point polygon";
  } else
    o << "Empty polygon";
}

std::ostream &operator<<(std::ostream &o, const lt_polys_t &polygons) {
  PrintLatticePolygons(o, polygons);
  return o;
}

// Computes the area of a _convex_ integer lattice polygon
//  template<typename INTEGERMATRIX>
using INTEGERMATRIX = Eigen::MatrixXi;
unsigned int TwiceArea2DLatticePolygon(const INTEGERMATRIX &lt_poly) {
  LF_VERIFY_MSG(lt_poly.rows() == 2, "Only available in dimension 2");
  int n_vertices = lt_poly.cols();
  if (n_vertices < 3)
    return 0;
  else if (n_vertices == 3) {
    // Determinant area formula for triangles
    Eigen::Vector2i d1((lt_poly.col(1) - lt_poly.col(0)));
    Eigen::Vector2i d2((lt_poly.col(2) - lt_poly.col(0)));
    return std::abs((int)(d1(0) * d2(1) - d1(1) * d2(0)));
  } else {
    // recursion: split off a triangle
    unsigned int area_triangle =
        TwiceArea2DLatticePolygon(lt_poly.template block<2, 3>(0, 0));
    unsigned int area_remainder =
        TwiceArea2DLatticePolygon(lt_poly.block(0, 1, 2, n_vertices - 1));
    return area_triangle + area_remainder;
  }
  return 0;
}

// Summing 2x area of several polygons
unsigned int SumTwiceAreaPolygons(const lt_polys_t &polygons) {
  int n_poly = polygons.size();
  if (n_poly > 0) {
    unsigned int area_sum = 0;
    for (int l = 0; l < n_poly; l++) {
      area_sum += TwiceArea2DLatticePolygon(polygons[l]);
    }
    return area_sum;
  }
  return 0;
}

// Test of topological refinement of a triangle
TEST(TopRefTest, TopSegRef) {
  // Set up refinement pattern
  Hybrid2DRefinementPattern rp_copy(lf::base::RefEl::kSegment(),
                                    RefPat::rp_copy);
  Hybrid2DRefinementPattern rp_split(lf::base::RefEl::kSegment(),
                                     RefPat::rp_split);

  std::cout << "COPY(N=" << rp_copy.LatticeConst()
            << "): " << rp_copy.ChildPolygons(0) << std::endl;
  std::cout << "SPLIT(N=" << rp_copy.LatticeConst()
            << "): " << rp_split.ChildPolygons(0) << std::endl;
}


TEST(TopRefTest, TopTriaRef) {
  // Copy "refinement"
  {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kTria(), RefPat::rp_copy);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "COPY: " << child_polygons << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
  // Bisection refinement
  for (int anchor = 0; anchor < 3; anchor++) {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kTria(), RefPat::rp_bisect,
                                 anchor);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "BISECT(anchor = " << anchor << ") : " << child_polygons
              << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
  // Trisection refinement
  for (int anchor = 0; anchor < 3; anchor++) {
    Hybrid2DRefinementPattern rp1(lf::base::RefEl::kTria(), RefPat::rp_trisect,
                                  anchor);
    Hybrid2DRefinementPattern rp2(lf::base::RefEl::kTria(),
                                  RefPat::rp_trisect_left, anchor);
    lt_polys_t child_polygons1(rp1.ChildPolygons(0));
    lt_polys_t child_polygons2(rp2.ChildPolygons(0));
    std::cout << "TRISECT(anchor=" << anchor << ") : " << child_polygons1
              << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons1)
              << std::endl;

    std::cout << "TRISECT_LEFT(anchor=" << anchor << ") : " << child_polygons2
              << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons2)
              << std::endl;
  }
  // Splitting into four triangles
  for (int anchor = 0; anchor < 3; anchor++) {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kTria(), RefPat::rp_quadsect,
                                 anchor);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "QUADSECT(anchor=" << anchor << ") : " << child_polygons
              << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
  // Regular refinement
  {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kTria(), RefPat::rp_regular);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "REGULAR: " << child_polygons << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
  // Bayrcentric refinement
  {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kTria(),
                                 RefPat::rp_barycentric);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "BARYCENTRIC: " << child_polygons << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
}

TEST(TopRefTest, TopQuadRef) {
  // Check the various refinements
  {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kQuad(), RefPat::rp_copy);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "COPY: " << child_polygons << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
  for (int anchor = 0; anchor < 4; anchor++) {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kQuad(), RefPat::rp_trisect,
                                 anchor);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "TRISECT(anchor=" << anchor << ") : " << child_polygons
              << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
  for (int anchor = 0; anchor < 4; anchor++) {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kQuad(), RefPat::rp_quadsect,
                                 anchor);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "QUADSECT(anchor=" << anchor << ") : " << child_polygons
              << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
  for (int anchor = 0; anchor < 4; anchor++) {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kQuad(), RefPat::rp_split,
                                 anchor);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "Split(anchor=" << anchor << ") : " << child_polygons
              << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
  for (int anchor = 0; anchor < 4; anchor++) {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kQuad(), RefPat::rp_threeedge,
                                 anchor);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "THREEEDGE(anchor=" << anchor << ") : " << child_polygons
              << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
  {
    Hybrid2DRefinementPattern rp(lf::base::RefEl::kQuad(), RefPat::rp_regular);
    lt_polys_t child_polygons(rp.ChildPolygons(0));
    std::cout << "BARYCENTRIC: " << child_polygons << std::endl;
    std::cout << "Area = " << SumTwiceAreaPolygons(child_polygons) << std::endl;
  }
  } 

}  // namespace lf::refinement::test
