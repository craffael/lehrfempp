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
#include <whitney_one_form_curl_element_matrix_provider.h>
#include <whitney_one_form_grad_element_matrix_provider.h>

#include <array>
#include <cmath>

TEST(projects_hldo_sphere_assembly,
     whitney_one_form_curl_matrix_assembler_test_two) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices << 0, 5, 3,
              0, -1, 2,
              1,  1, 1;
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

  // Build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));
  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element matrix for the triangle
  const auto elem_mat_provider = projects::hldo_sphere::assemble::
      WhitneyOneFormCurlElementMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal <<  1, -1, -1,
             -1, 1,  1,
             -1, 1,  1;
  // clang-format on
  Ae_anal *= 2. / 13.;
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_DOUBLE_EQ(Ae(i, j), Ae_anal(i, j))
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

TEST(projects_hldo_sphere_assembly,
     whitney_one_form_curl_matrix_assembler_test_one) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices << 0, 5, 2,
              0, 0, 1,
              1, 1, 1;
  // clang-format on
  lf::mesh::hybrid2d::MeshFactory factory(3);
  factory.AddPoint(vertices.col(0));
  factory.AddPoint(vertices.col(1));
  factory.AddPoint(vertices.col(2));

  // Build Segments
  // First define directions
  std::vector<std::array<lf::mesh::MeshFactory::size_type, 2>> edge_nodes(3);
  edge_nodes[0] = {0, 1};
  edge_nodes[1] = {2, 1};
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

  // build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));

  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element matrix for the triangle
  const auto elem_mat_provider = projects::hldo_sphere::assemble::
      WhitneyOneFormCurlElementMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal <<  1, -1,  1,
             -1,  1, -1,
              1, -1,  1;
  // clang-format on
  Ae_anal *= 2. / 5.;
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_DOUBLE_EQ(Ae(i, j), Ae_anal(i, j))
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

TEST(projects_hldo_sphere_assembly,
     whitney_one_form_grad_matrix_assembler_test_ref_2d) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(2, 3);
  // clang-format off
  vertices << 0, 1, 0,
              0, 0, 1;
  // clang-format on
  lf::mesh::hybrid2d::MeshFactory factory(2);
  factory.AddPoint(vertices.col(0));
  factory.AddPoint(vertices.col(1));
  factory.AddPoint(vertices.col(2));

  // Build Segments
  // First define directions
  std::vector<std::array<lf::mesh::MeshFactory::size_type, 2>> edge_nodes(3);
  edge_nodes[0] = {0, 1};
  edge_nodes[1] = {1, 2};
  edge_nodes[2] = {0, 2};

  // Build segments based on directions
  std::vector<Eigen::MatrixXd> edge_endpoints(3);
  std::vector<std::unique_ptr<lf::geometry::Geometry>> edge_geom(3);
  for (int i = 0; i < 3; i++) {
    edge_endpoints[i] = Eigen::MatrixXd::Zero(2, 2);
    edge_endpoints[i].col(0) = vertices.col(edge_nodes[i][0]);
    edge_endpoints[i].col(1) = vertices.col(edge_nodes[i][1]);
    edge_geom[i] = std::make_unique<lf::geometry::SegmentO1>(edge_endpoints[i]);
    factory.AddEntity(seg, nonstd::span(edge_nodes[i].data(), 2),
                      std::move(edge_geom[i]));
  }

  // build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));
  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element matrix for the triangle
  const auto elem_mat_provider = projects::hldo_sphere::assemble::
      WhitneyOneFormGradElementMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal = Eigen::MatrixXd::Ones(3, 3);
  // clang-format off
  Ae_anal <<   -3,  2,  1,
                0, -1,  1,
               -3,  1,  2;
  // clang-format on
  Ae_anal *= 1. / 6.;
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_DOUBLE_EQ(Ae(i, j), Ae_anal(i, j))
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

TEST(projects_hldo_sphere_assembly,
     whitney_one_form_grad_matrix_assembler_test_3d) {
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

  // build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));

  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element matrix for the triangle
  const auto elem_mat_provider = projects::hldo_sphere::assemble::
      WhitneyOneFormGradElementMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal << -1,  1, 0,
              0, -1, 1,
             -1,  0, 1;
  // clang-format on
  Ae_anal *= 1. / 2. / std::sqrt(3);
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_DOUBLE_EQ(Ae(i, j), Ae_anal(i, j))
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}