#include <gtest/gtest.h>
#include <laplace_matrix_provider.h>
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
#include <load_vector_provider.h>
#include <mass_matrix_provider.h>

#include <array>
#include <cmath>

// helper functions to directly use mathematica output
auto Power = [](double a, double b) -> double { return pow(a, b); };
auto Sin = [](double a) -> double { return sin(a); };
auto Cos = [](double a) -> double { return cos(a); };
auto Sqrt = [](double a) -> double { return sqrt(a); };

/**
 *
 * @brief Test the element matrix provider for the mass matrix of
 * the whitney zero form
 *
 * Test the element vector provider on the triangle
 * ((1,0,0), (0,1,0), (0,0,1))
 *
 * with relative edge orientations
 *
 * s_0 = 1, s_1 = -1, s_2 = 1
 *
 * And the linear form
 *
 * @f[
 *  \int_K\ u \cdot v \,\mathrm{d}x, \quad u, v \in H^1
 * @f]
 *
 * Using the whitney zero forms (barycentric basis functions)
 *
 *
 * For analytical solution see assemby.pdf
 *
 */
TEST(projects_hldo_sphere_assembly, laplace_matrix_assembler_test_regular) {
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

  // Build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));
  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element matrix for the triangle
  const auto elem_mat_provider =
      projects::hldo_sphere::assemble::LaplaceMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal <<  2, -1, -1,
             -1, 2,  -1,
             -1, -1,  2;
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

/**
 *
 * @brief Test the element matrix provider for the mass matrix of
 * the whitney zero form
 *
 * Test the element vector provider on the triangle
 * ((1,0,0), (0,1,0), (0,0,1))
 *
 * with relative edge orientations
 *
 * s_0 = 1, s_1 = -1, s_2 = 1
 *
 * And the linear form
 *
 * @f[
 *  \int_K\ u \cdot v \,\mathrm{d}x, \quad u, v \in H^1
 * @f]
 *
 * Using the whitney zero forms (barycentric basis functions)
 *
 * For analytical solution see assemby.pdf
 *
 */
TEST(projects_hldo_sphere_assembly, mass_element_matrix_provider_test_regular) {
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
      projects::hldo_sphere::assemble::MassMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal <<  2, 1, 1,
              1, 2, 1,
              1, 1, 2;
  // clang-format on
  Ae_anal *= 1. / 8. / std::sqrt(3);
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

/**
 *
 * @brief Test the element vector provider for the zero-form
 *
 * Test the element vector provider on the triangle
 * ((1,0,0), (0,1,0), (0,0,1))
 *
 * with relative edge orientations
 *
 * s_0 = -1, s_1 = 1, s_2 = -1
 *
 * And the linear form
 *
 * @f[
 *  \int_K\ f \cdot v \,\mathrm{d}x, \quad f, v \in H^1
 * @f]
 *
 * using the load function
 *
 * @f[
 *  f((x_0,x_1, x_2)) = x_0 / \|x\|
 * @f]
 *
 * And the barycentric basis functions
 *
 * For analytical solution see assemby.pdf
 *
 */
TEST(projects_hldo_sphere_assembly, load_vector_provider_test_regular) {
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
  auto f = [&](const Eigen::Vector3d x) -> double { return x(0) / x.norm(); };

  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element vec for the triangle
  const auto elem_vec_provider =
      projects::hldo_sphere::assemble::LoadVectorProvider(f);
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

/**
 *
 * @brief Test the element matrix provider for the mass matrix of
 * the whitney zero form
 *
 * Test the element vector provider on the triangle
 * ((1,1/2,0), (0,1/2,1/2), (0,0,1/2))
 *
 * with relative edge orientations
 *
 * s_0 = 1, s_1 = -1, s_2 = 1
 *
 * And the linear form
 *
 * @f[
 *  \int_K\ u \cdot v \,\mathrm{d}x, \quad u, v \in H^1
 * @f]
 *
 * Using the whitney zero forms (barycentric basis functions)
 *
 *
 * For analytical solution see assemby.pdf
 *
 */
TEST(projects_hldo_sphere_assembly, laplace_matrix_assembler_test_flat) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices << 1,   0,   0,
            0.5, 0.5,   0,
              0, 0.5, 0.5;
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

  // Build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));
  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element matrix for the triangle
  const auto elem_mat_provider =
      projects::hldo_sphere::assemble::LaplaceMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);

  // auto generated with mathematica
  Ae_anal << 1 / (2. * Sqrt(5)), -0.5 * 1 / Sqrt(5), 0, -0.5 * 1 / Sqrt(5),
      3 / Sqrt(5), -0.5 * Sqrt(5), 0, -0.5 * Sqrt(5), Sqrt(5) / 2.;
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae(i, j), Ae_anal(i, j), 1e-13)
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

/**
 *
 * @brief Test the element matrix provider for the mass matrix of
 * the whitney zero form
 *
 * Test the element vector provider on the triangle
 * ((1,0,0), (1/2,1/2,0), (0,1/2,1/2))
 *
 * with relative edge orientations
 *
 * s_0 = 1, s_1 = -1, s_2 = 1
 *
 * And the linear form
 *
 * @f[
 *  \int_K\ u \cdot v \,\mathrm{d}x, \quad u, v \in H^1
 * @f]
 *
 * Using the whitney zero forms (barycentric basis functions)
 *
 * For analytical solution see assemby.pdf
 *
 */
TEST(projects_hldo_sphere_assembly, mass_element_matrix_provider_test_flat) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices << 1,   0,   0,
            0.5, 0.5,   0,
              0, 0.5, 0.5;
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
  const auto elem_mat_provider =
      projects::hldo_sphere::assemble::MassMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);

  // autogenerated with mathematica
  Ae_anal << Sqrt(5) / 48., Sqrt(5) / 96., Sqrt(5) / 96., Sqrt(5) / 96.,
      Sqrt(5) / 48., Sqrt(5) / 96., Sqrt(5) / 96., Sqrt(5) / 96., Sqrt(5) / 48.;

  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae(i, j), Ae_anal(i, j), 1e-13)
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

/**
 *
 * @brief Test the element vector provider for the zero-form
 *
 * Test the element vector provider on the triangle
 * ((1,0,0), (1/2,1/2,0), (0,1/2,1/2))
 *
 * with relative edge orientations
 *
 * s_0 = -1, s_1 = 1, s_2 = -1
 *
 * And the linear form
 *
 * @f[
 *  \int_K\ f \cdot v \,\mathrm{d}x, \quad f, v \in H^1
 * @f]
 *
 * using the load function
 *
 * @f[
 *  f((x_0,x_1, x_2)) = Sin[x + y]
 * @f]
 *
 * And the barycentric basis functions
 *
 * For analytical solution see assemby.pdf
 *
 */
TEST(projects_hldo_sphere_assembly, load_vector_provider_test_flat) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices << 1,   0,   0,
            0.5, 0.5,   0,
              0, 0.5, 0.5;
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
  auto f = [&](const Eigen::Vector3d x) -> double { return sin(x(0) + x(1)); };

  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element vec for the triangle
  const auto elem_vec_provider =
      projects::hldo_sphere::assemble::LoadVectorProvider(f);
  const Eigen::VectorXd Ae = elem_vec_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::VectorXd Ae_anal(3);

  // autogenerated by mathematica
  Ae_anal << (Sqrt(5) * (Sin(0.75) / 2. + Sin(1) / 2.)) / 24.,
      (Sqrt(5) * (Sin(0.25) / 2. + Sin(1) / 2.)) / 24.,
      (Sqrt(5) * (Sin(0.25) / 2. + Sin(0.75) / 2.)) / 24.;

  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    EXPECT_DOUBLE_EQ(Ae(i), Ae_anal(i)) << "mismatch in entry (" << i << ")";
  }
}
