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
#include <rot_w_one_form_div_element_matrix_provider.h>
#include <rot_w_one_form_dot_element_matrix_provider.h>
#include <whitney_two_mass_matrix_provider.h>

#include <array>
#include <cmath>

TEST(projects_hldo_sphere_assembly, rot_w_div_form_dot_matrix_assembler_test) {
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
  const auto elem_mat_provider =
      projects::hldo_sphere::assemble::RotWOneFormDotElementMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal <<  5, -1, 1,
             -1,  5, 1,
              1,  1, 5;
  // clang-format on
  Ae_anal *= std::sqrt(3) / 36;
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

TEST(projects_hldo_sphere_assembly, rot_w_one_form_div_matrix_assembler_test) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices << 5, 4, 10,
              1, 0, 30,
              0, 1,  2;
  // clang-format on
  lf::mesh::hybrid2d::MeshFactory factory(3);
  factory.AddPoint(vertices.col(0));
  factory.AddPoint(vertices.col(1));
  factory.AddPoint(vertices.col(2));

  // Build Segments
  // First define directions
  std::vector<std::array<lf::mesh::MeshFactory::size_type, 2>> edge_nodes(3);
  edge_nodes[0] = {1, 0};
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
  const auto elem_mat_provider =
      projects::hldo_sphere::assemble::RotWOneFormDivElementMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 1);
  // clang-format off
  Ae_anal <<  1, 1, -1;
  // clang-format on
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    EXPECT_DOUBLE_EQ(Ae(i), Ae_anal(i)) << "mismatch in entry (" << i << ")";
  }
}

TEST(projects_hldo_sphere_assembly, whitney_two_mass_matrix_assembler_test) {
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
  const auto elem_mat_provider =
      projects::hldo_sphere::assemble::WhitneyTwoMassMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal = Eigen::MatrixXd::Ones(1, 1);
  Ae_anal *= sqrt(3) / 2.;
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  EXPECT_DOUBLE_EQ(Ae(0), Ae_anal(0));
}
