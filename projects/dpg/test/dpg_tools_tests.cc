/**
 * @file
 * @brief Test of dpg tools
 * @author Philippe Peter
 * @date June 2019
 * @copyright MIT License
 */

#include <gtest/gtest.h>
#include <lf/io/io.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "../dpg_tools.h"

namespace projects::dpg::test {

TEST(OuterNormalVectors, TRIA) {
  std::cout << ">>> Testing computation of outer normal vectors on triangular "
               "elements: \n";

  // right angled triangle, counter clockwise
  // numbering of nodes.
  {
    Eigen::MatrixXd global_nodes(2, 3);
    global_nodes << 1, 0, 1, 1, 1, 0;

    Eigen::MatrixXd exact_normals(2, 3);
    double sqrt2h = std::sqrt(2) / 2.0;
    exact_normals << 0, -sqrt2h, 1, 1, -sqrt2h, 0;

    lf::geometry::TriaO1 geom(global_nodes);
    Eigen::MatrixXd computed_normals = OuterNormals(geom);
    EXPECT_TRUE(exact_normals.isApprox(computed_normals));

    std::cout << "exact normals: \n" << exact_normals << std::endl;
    std::cout << "computed normals: \n" << computed_normals << std::endl;
  }

  // right angled triangle, clockwise
  // numbering of nodes.
  {
    Eigen::MatrixXd global_nodes(2, 3);
    global_nodes << 1, 1, 0, 1, 0, 1;

    Eigen::MatrixXd exact_normals(2, 3);
    double sqrt2h = std::sqrt(2) / 2.0;
    exact_normals << 1, -sqrt2h, 0, 0, -sqrt2h, 1;

    lf::geometry::TriaO1 geom(global_nodes);
    Eigen::MatrixXd computed_normals = OuterNormals(geom);
    EXPECT_TRUE(exact_normals.isApprox(computed_normals));

    std::cout << "exact normals: \n" << exact_normals << std::endl;
    std::cout << "computed normals: \n" << computed_normals << std::endl;
  }

  // right angled triangle, counter clockwise numbering with offset:
  {
    Eigen::MatrixXd global_nodes(2, 3);
    global_nodes << 11, 10, 11, 11, 11, 10;

    Eigen::MatrixXd exact_normals(2, 3);
    double sqrt2h = std::sqrt(2) / 2.0;
    exact_normals << 0, -sqrt2h, 1, 1, -sqrt2h, 0;

    lf::geometry::TriaO1 geom(global_nodes);
    Eigen::MatrixXd computed_normals = OuterNormals(geom);
    EXPECT_TRUE(exact_normals.isApprox(computed_normals));

    std::cout << "exact normals: \n" << exact_normals << std::endl;
    std::cout << "computed normals: \n" << computed_normals << std::endl;
  }
}

TEST(OuterNormalVectors, QUAD) {
  std::cout << ">>> Testing computation of outer normal vectors on "
               "quadrilateral elements: \n";

  // convex quad, counter clockwise
  // numbering of nodes.
  {
    Eigen::MatrixXd global_nodes(2, 4);
    global_nodes << 1, 2, 3, 2, 1, 1, 2, 2;

    Eigen::MatrixXd exact_normals(2, 4);
    double sqrt2h = std::sqrt(2) / 2.0;
    exact_normals << 0, sqrt2h, 0, -sqrt2h, -1, -sqrt2h, 1, sqrt2h;

    lf::geometry::QuadO1 geom(global_nodes);
    Eigen::MatrixXd computed_normals = OuterNormals(geom);
    EXPECT_TRUE(exact_normals.isApprox(computed_normals));

    std::cout << "exact normals: \n" << exact_normals << std::endl;
    std::cout << "computed normals: \n" << computed_normals << std::endl;
  }

  // convex quad, clockwise
  // numbering of nodes.
  {
    Eigen::MatrixXd global_nodes(2, 4);
    global_nodes << 1, 2, 3, 2, 1, 2, 2, 1;

    Eigen::MatrixXd exact_normals(2, 4);
    double sqrt2h = std::sqrt(2) / 2.0;
    exact_normals << -sqrt2h, 0, sqrt2h, 0, sqrt2h, 1, -sqrt2h, -1;

    lf::geometry::QuadO1 geom(global_nodes);
    Eigen::MatrixXd computed_normals = OuterNormals(geom);
    EXPECT_TRUE(exact_normals.isApprox(computed_normals));

    std::cout << "exact normals: \n" << exact_normals << std::endl;
    std::cout << "computed normals: \n" << computed_normals << std::endl;
  }

  // non-convex quad, counter clockwise
  // numbering of nodes.
  {
    Eigen::MatrixXd global_nodes(2, 4);
    global_nodes << 1, 0, 0, 0.5, 0, 1, -1, 0;

    Eigen::MatrixXd exact_normals(2, 4);
    double sqrt2h = std::sqrt(2) / 2.0;
    double sqrt5h = std::sqrt(5) / 2.0;
    exact_normals << sqrt2h, -1, 1 / sqrt5h, 0, sqrt2h, 0, -0.5 / sqrt5h, -1;

    lf::geometry::QuadO1 geom(global_nodes);
    Eigen::MatrixXd computed_normals = OuterNormals(geom);
    EXPECT_TRUE(exact_normals.isApprox(computed_normals));

    std::cout << "exact normals: \n" << exact_normals << std::endl;
    std::cout << "computed normals: \n" << computed_normals << std::endl;
  }

  // non-convex quad, clockwise
  // numbering of nodes.
  {
    Eigen::MatrixXd global_nodes(2, 4);
    global_nodes << 1, 0.5, 0, 0, 0, 0, -1, 1;

    Eigen::MatrixXd exact_normals(2, 4);
    double sqrt2h = std::sqrt(2) / 2.0;
    double sqrt5h = std::sqrt(5) / 2.0;
    exact_normals << 0, 1 / sqrt5h, -1, sqrt2h, -1, -0.5 / sqrt5h, 0, sqrt2h;

    lf::geometry::QuadO1 geom(global_nodes);
    Eigen::MatrixXd computed_normals = OuterNormals(geom);
    EXPECT_TRUE(exact_normals.isApprox(computed_normals));

    std::cout << "exact normals: \n" << exact_normals << std::endl;
    std::cout << "computed normals: \n" << computed_normals << std::endl;
  }
}

TEST(InOutFlowBoundary, assignement_test_1) {
  std::cout << ">>>InOutFlowBoundary: Testing unique assignement of edges to "
               "in/outflow boundary \n";
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  auto beta = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> Eigen::Vector2d {
        return (Eigen::VectorXd(2) << 1.0, 0.0).finished();
      });

  auto inflow_boundary = flagEntitiesOnInflowBoundary(mesh_p, beta);
  auto outflow_boundary = flagEntitiesOnOutflowBoundary(mesh_p, beta);

  for (const lf::mesh::Entity* const edge : mesh_p->Entities(1)) {
    LF_ASSERT_MSG(
        !(inflow_boundary(*edge) && outflow_boundary(*edge)),
        "edge " << edge << " assigned to inflow- and outflow boundary");
  }
}

/*
TEST(BoundaryQuadRule, consistency_test){
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  int degree = 3;
  lf::quad::QuadRule segment_qr =
lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), degree); std::cout <<
"segment quadrule, " << "Points: \n" << segment_qr.Points() << std::endl <<
               "Weights: \n" << segment_qr.Weights() << std::endl;

  std::array<std::vector<lf::quad::QuadRule>,5> qrs;

  for(auto ref_el: {lf::base::RefEl::kTria(), lf::base::RefEl::kQuad()}){
    std::cout << "Boundary quadrule on ref_el " << ref_el << std::endl;
    auto boundary_qr = internal::BoundaryQuadRule(ref_el, segment_qr);
    int num_segments = ref_el.NumSubEntities(1);
    int idx = ref_el.Id();
    qrs[ref_el.Id()].resize(num_segments);
    for(int segment = 0; segment < num_segments; segment++){
      qrs[ref_el.Id()][segment] = boundary_qr[segment];
      std::cout <<  "segment " << segment  <<
                 std::endl << "Points: \n" << qrs[ref_el.Id()][segment].Points()
<< std::endl
                << "Weights: \n" << qrs[ref_el.Id()][segment].Weights() <<
std::endl;
    }
  }
  for(const lf::mesh::Entity& cell: mesh_p->Entities(0)){
    lf::geometry::Geometry* geo_ptr = cell.Geometry();

    for(int segment = 0; segment < cell.RefEl().NumSubEntities(1);segment++){
      auto& qr = qrs[cell.RefEl().Id()][segment];
      auto segment_geo_ptr = geo_ptr->SubGeometry(1,segment);
      const Eigen::VectorXd cell_dets =
geo_ptr->IntegrationElement(qr.Points()); const Eigen::VectorXd segment_dets =
segment_geo_ptr->IntegrationElement(segment_qr.Points()); std::cout << "expected
eaulity: \n " << cell_dets.cwiseProduct(qr.Weights()) << "\n <-> \n"
                << segment_dets.cwiseProduct(segment_qr.Weights()) << "\n";

    }
  }
}
*/

}  // namespace projects::dpg::test
