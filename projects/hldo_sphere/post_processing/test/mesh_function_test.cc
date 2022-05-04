#include <gtest/gtest.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <mesh_function_whitney_one.h>
#include <mesh_function_whitney_two.h>
#include <mesh_function_whitney_zero.h>

#include <Eigen/Dense>

TEST(projects_hldo_sphere_post_processing, mesh_function_whitney_one_basic) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices << 1, 0, 0, 
              0, 1, 0,
              0, 0, 1;
  // clang-format on
  lf::mesh::hybrid2d::MeshFactory factory(3);
  factory.AddPoint(vertices.col(0));
  factory.AddPoint(vertices.col(1));
  factory.AddPoint(vertices.col(2));

  // Build Segments
  // First define directions
  std::vector<std::array<lf::mesh::MeshFactory::size_type, 2>> edge_nodes(3);
  edge_nodes[0] = {0, 1};
  edge_nodes[1] = {1, 2};
  edge_nodes[2] = {2, 0};

  // Build segments based on directions
  std::vector<Eigen::MatrixXd> edge_endpoints(3);
  std::vector<std::unique_ptr<lf::geometry::Geometry>> edge_geom(3);
  for (int i = 0; i < 3; i++) {
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

  // build mesh
  const auto mesh = factory.Build();

  // Create mesh function
  Eigen::VectorXd mu(3);
  mu << 1, 2, 3;
  projects::hldo_sphere::post_processing::MeshFunctionWhitneyOne mf_w_one(mu,
                                                                          mesh);

  Eigen::MatrixXd eval_points(2, 2);
  // clang-format off
  eval_points << 0, 0.5, 
                 0, 0.5;
  // clang-format on

  std::vector<Eigen::VectorXd> evaluations_tria0 =
      mf_w_one(*(mesh->Entities(0)[0]), eval_points);

  // We use the precomputed gradients off the tests in the assembly
  // for the point 0,0 we get on the first triangle the function
  // grad(labmda_1) = (-1/3, 2/3, -1/3),
  // 0,
  // -grad(lambda_2) = (1/3, 1/3, -2/3)
  // for  the point 0.5,0.5 we get something slightly more complicated
  // -0.5 grad(lambda_0) = (-2/6, 1/6, 1/6)
  // 0.5 grad(lambda_2) - 0.5 grad(lambda_1) = (0,-3/6, 3/6)
  // 0.5 grad(lambda_0) = (2/6, -1/6, -1/6)
  Eigen::MatrixXd eval_tria0_ana0(3, 1);
  // clang-format off
  eval_tria0_ana0 <<  2,
                      5,
                     -7;
  // clang-format on
  eval_tria0_ana0 *= 1 / 3.;

  Eigen::MatrixXd eval_tria0_ana1(3, 1);
  // clang-format off
  eval_tria0_ana1 <<  4,
                     -8,
                      4;
  eval_tria0_ana1 *= 1 / 6.;
  // clang-format on
  for (int c = 0; c < evaluations_tria0[0].cols(); c++) {
    ASSERT_NEAR(evaluations_tria0[0](c), eval_tria0_ana0(c), 1e-10)
        << "Entries (" << c << ") do not match\n";
  }

  for (int c = 0; c < evaluations_tria0[1].cols(); c++) {
    ASSERT_NEAR(evaluations_tria0[1](c), eval_tria0_ana1(c), 1e-10)
        << "Entries (" << c << ") do not match\n";
  }
}

TEST(projects_hldo_sphere_post_processing, mesh_function_whitney_zero_basic) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices << 1, 0, 0, 
              0, 1, 0,
              0, 0, 1;
  // clang-format on
  lf::mesh::hybrid2d::MeshFactory factory(3);
  factory.AddPoint(vertices.col(0));
  factory.AddPoint(vertices.col(1));
  factory.AddPoint(vertices.col(2));

  // Build Segments
  // First define directions
  std::vector<std::array<lf::mesh::MeshFactory::size_type, 2>> edge_nodes(3);
  edge_nodes[0] = {0, 1};
  edge_nodes[1] = {1, 2};
  edge_nodes[2] = {2, 0};

  // Build segments based on directions
  std::vector<Eigen::MatrixXd> edge_endpoints(3);
  std::vector<std::unique_ptr<lf::geometry::Geometry>> edge_geom(3);
  for (int i = 0; i < 3; i++) {
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

  // build mesh
  const auto mesh = factory.Build();

  // Create mesh function
  Eigen::VectorXd mu(3);
  mu << 1, 2, 3;
  projects::hldo_sphere::post_processing::MeshFunctionWhitneyZero mf_w_zero(
      mu, mesh);

  Eigen::MatrixXd eval_points(2, 3);
  // clang-format off
  eval_points << 0, 0.5, 0.25, 
                 0, 0.5, 0.25;
  // clang-format on

  std::vector<double> evaluations_tria0 =
      mf_w_zero(*(mesh->Entities(0)[0]), eval_points);

  Eigen::MatrixXd eval_tria0_ana(3, 1);
  // clang-format off
  eval_tria0_ana <<    1,
                     2.5, 
                    1.75;
  // clang-format on

  for (int c = 0; c < eval_tria0_ana.cols(); c++) {
    ASSERT_NEAR(evaluations_tria0[c], eval_tria0_ana(c), 1e-10)
        << "Entries (" << c << ") do not match\n";
  }
}

TEST(projects_hldo_sphere_post_processing, mesh_function_whitney_two_basic) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices << 1, 0, 0, 
              0, 1, 0,
              0, 0, 1;
  // clang-format on
  lf::mesh::hybrid2d::MeshFactory factory(3);
  factory.AddPoint(vertices.col(0));
  factory.AddPoint(vertices.col(1));
  factory.AddPoint(vertices.col(2));

  // Build Segments
  // First define directions
  std::vector<std::array<lf::mesh::MeshFactory::size_type, 2>> edge_nodes(3);
  edge_nodes[0] = {0, 1};
  edge_nodes[1] = {1, 2};
  edge_nodes[2] = {2, 0};

  // Build segments based on directions
  std::vector<Eigen::MatrixXd> edge_endpoints(3);
  std::vector<std::unique_ptr<lf::geometry::Geometry>> edge_geom(3);
  for (int i = 0; i < 3; i++) {
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

  // build mesh
  const auto mesh = factory.Build();

  // Create mesh function
  Eigen::VectorXd mu(1);
  mu << 3;
  projects::hldo_sphere::post_processing::MeshFunctionWhitneyTwo mf_w_two(mu,
                                                                          mesh);

  Eigen::MatrixXd eval_points(2, 3);
  // clang-format off
  eval_points << 0, 0.5, 0.25, 
                 0, 0.5, 0.25;
  // clang-format on

  std::vector<double> evaluations_tria0 =
      mf_w_two(*(mesh->Entities(0)[0]), eval_points);

  Eigen::MatrixXd eval_tria0_ana(3, 1);
  // clang-format off
  eval_tria0_ana <<    3,
                       3, 
                       3;
  // clang-format on

  for (int c = 0; c < eval_tria0_ana.cols(); c++) {
    ASSERT_NEAR(evaluations_tria0[c], eval_tria0_ana(c), 1e-10)
        << "Entries (" << c << ") do not match\n";
  }
}
