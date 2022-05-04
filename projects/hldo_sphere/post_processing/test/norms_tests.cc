#include <gtest/gtest.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <norms.h>

#include <Eigen/Dense>

TEST(projects_hldo_sphere_post_processing, norm_l2_scalar_const) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 4);
  // clang-format off
  vertices << 1, 0, 0, 0,
              0, 1, 0,-1,
              0, 0, 1, 0;
  // clang-format on
  lf::mesh::hybrid2d::MeshFactory factory(3);
  factory.AddPoint(vertices.col(0));
  factory.AddPoint(vertices.col(1));
  factory.AddPoint(vertices.col(2));
  factory.AddPoint(vertices.col(3));

  // Build Segments
  // First define directions
  std::vector<std::array<lf::mesh::MeshFactory::size_type, 2>> edge_nodes(5);
  edge_nodes[0] = {0, 1};
  edge_nodes[1] = {1, 2};
  edge_nodes[2] = {2, 0};
  edge_nodes[3] = {2, 3};
  edge_nodes[4] = {3, 0};

  // Build segments based on directions
  std::vector<Eigen::MatrixXd> edge_endpoints(5);
  std::vector<std::unique_ptr<lf::geometry::Geometry>> edge_geom(5);
  for (int i = 0; i < 5; i++) {
    edge_endpoints[i] = Eigen::MatrixXd::Zero(3, 2);
    edge_endpoints[i].col(0) = vertices.col(edge_nodes[i][0]);
    edge_endpoints[i].col(1) = vertices.col(edge_nodes[i][1]);
    edge_geom[i] = std::make_unique<lf::geometry::SegmentO1>(edge_endpoints[i]);
    factory.AddEntity(seg, nonstd::span(edge_nodes[i].data(), 2),
                      std::move(edge_geom[i]));
  }

  // Build triangles
  Eigen::MatrixXd triag1 = vertices.leftCols(3);
  std::unique_ptr<lf::geometry::Geometry> geom1 =
      std::make_unique<lf::geometry::TriaO1>(triag1);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes1 = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes1.data(), 3), std::move(geom1));
  std::make_unique<lf::geometry::TriaO1>(triag1);

  Eigen::MatrixXd triag2(3, 3);
  triag2.col(0) = vertices.col(0);
  triag2.col(1) = vertices.col(2);
  triag2.col(2) = vertices.col(3);
  std::unique_ptr<lf::geometry::Geometry> geom2 =
      std::make_unique<lf::geometry::TriaO1>(triag2);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes2 = {0, 2, 3};
  factory.AddEntity(trig, nonstd::span(nodes2.data(), 3), std::move(geom2));
  std::make_unique<lf::geometry::TriaO1>(triag2);

  // build mesh
  const auto mesh = factory.Build();

  // create function f
  auto f_zero = [](Eigen::Vector3d x) -> double { return 0; };
  lf::mesh::utils::MeshFunctionGlobal<decltype(f_zero)> mf_f_zero(f_zero);
  // create the square funciton of scalar values
  auto sq_f = [](double x) -> double { return x * x; };

  // create the quadrature rule which is exact for cosnt functions
  lf::quad::QuadRule qr = lf::quad::make_TriaQR_EdgeMidpointRule();

  // Compute the L2 norm of the function constant zero on the mesh
  const double L2_zero =
      projects::hldo_sphere::post_processing::L2norm<decltype(mf_f_zero),
                                                     decltype(sq_f)>(
          mesh, mf_f_zero, sq_f, qr);
  ASSERT_NEAR(L2_zero, 0, 1e-10);

  // create function f
  auto f_const = [](Eigen::Vector3d x) -> double { return 5.0; };
  lf::mesh::utils::MeshFunctionGlobal<decltype(f_const)> mf_f_const(f_const);
  // Compute the L2 norm of the function constant zero on the mesh
  const double L2_const =
      projects::hldo_sphere::post_processing::L2norm<decltype(mf_f_const),
                                                     decltype(sq_f)>(
          mesh, mf_f_const, sq_f, qr);
  ASSERT_NEAR(L2_const, sqrt(sqrt(0.75) * 2 * 25.0), 1e-10);
}
