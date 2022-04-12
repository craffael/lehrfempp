#include <curl_element_vector_provider.h>
#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/base/base.h>
#include <lf/base/span.h>
#include <lf/geometry/geometry.h>
#include <lf/geometry/tria_o1.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/all_codim_mesh_data_set.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/quad/quad_rule.h>
#include <load_element_vector_provider.h>
#include <whitney_two_element_vector_provider.h>

#include <array>
#include <cmath>

TEST(projects_hldo_sphere_assembly, element_vector_provider_one_form) {
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
  edge_nodes[0] = {1, 0};
  edge_nodes[1] = {1, 2};
  edge_nodes[2] = {0, 2};

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

  // Build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));
  const auto mesh = factory.Build();

  // Define function f
  auto f = [&](const Eigen::Vector3d x) -> Eigen::Vector3d {
    Eigen::VectorXd res(3);
    res(0) = 1;
    res(1) = 1;
    res(2) = x(0);
    return res;
  };

  // define quad rule
  lf::quad::QuadRule quad{lf::quad::make_TriaQR_EdgeMidpointRule()};

  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element vec for the triangle
  const auto elem_vec_provider =
      projects::hldo_sphere::assemble::CurlElementVectorProvider(f, quad);
  const Eigen::VectorXd Ae = elem_vec_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::VectorXd Ae_anal(3);
  // clang-format off
  Ae_anal <<  std::sqrt(0.5), 
              3 * std::sqrt(0.5) - 6,
              5 * std::sqrt(0.5) - 6;
  // clang-format on
  Ae_anal *= std::sqrt(3) / 36.;
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    EXPECT_DOUBLE_EQ(Ae(i), Ae_anal(i)) << "mismatch in entry (" << i << ")";
  }
}

TEST(projects_hldo_sphere_assembly, element_vector_provider_two_form) {
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
  edge_nodes[0] = {1, 0};
  edge_nodes[1] = {1, 2};
  edge_nodes[2] = {0, 2};

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

  // Build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));
  const auto mesh = factory.Build();

  // Define function f
  auto f = [&](const Eigen::Vector3d x) -> double { return x(0); };

  // define quad rule
  lf::quad::QuadRule quad{lf::quad::make_TriaQR_EdgeMidpointRule()};

  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element vec for the triangle
  const auto elem_vec_provider =
      projects::hldo_sphere::assemble::WhitneyTwoElementVectorProvider(f, quad);
  const Eigen::VectorXd Ae = elem_vec_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::VectorXd Ae_anal(1);
  // clang-format off
  Ae_anal <<  2 * std::sqrt(0.5);
  // clang-format on
  Ae_anal *= std::sqrt(3) / 6;
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  EXPECT_DOUBLE_EQ(Ae(0), Ae_anal(0)) << "mismatch in entry (0)";
}

TEST(projects_hldo_sphere_assembly, load_vector_provider_zero_form) {
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
  edge_nodes[0] = {1, 0};
  edge_nodes[1] = {1, 2};
  edge_nodes[2] = {0, 2};

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

  // Build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));
  const auto mesh = factory.Build();

  // Define function f
  auto f = [&](const Eigen::Vector3d x) -> double { return x(0); };

  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element vec for the triangle
  const auto elem_vec_provider =
      projects::hldo_sphere::assemble::LoadElementVectorProvider(f);
  const Eigen::VectorXd Ae = elem_vec_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::VectorXd Ae_anal(3);
  // clang-format off
  Ae_anal <<  std::sqrt(0.5), 
              0.5 * std::sqrt(0.5),
              0.5 * std::sqrt(0.5);
  // clang-format on
  Ae_anal *= std::sqrt(3) / 6.;
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    EXPECT_DOUBLE_EQ(Ae(i), Ae_anal(i)) << "mismatch in entry (" << i << ")";
  }
}
