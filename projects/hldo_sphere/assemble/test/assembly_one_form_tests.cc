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
#include <sphere_triag_mesh_builder.h>
#include <whitney_one_curl_curl_matrix_provider.h>
#include <whitney_one_grad_matrix_provider.h>
#include <whitney_one_mass_matrix_provider.h>
#include <whitney_one_vector_provider.h>

#include <array>
#include <cmath>

/*********************************************************************
 * TESTS ON REGULAR TRIANGLE ((1,0,0), (0, 1, 0), (0, 0, 1))
 *********************************************************************/

/**
 *
 * @brief Test the element matrix provider for the whitney one form
 *
 * Given by the expression
 *
 * @f[
 *    \int\limits_K \text{curl}_{\Gamma} u \ \text{curl}_{\Gamma} v \ dx
 * @f]
 *
 * Test the element vector provider on the triangle
 * ((1,0,0), (0,1,0), (0,0,1))
 *
 * with relative edge orientations
 *
 * s_0 = 1, s_1 = 1, s_2 = -1
 *
 * And the linear form
 *
 * Using the whitney one forms for u and u
 *
 */
TEST(projects_hldo_sphere_assembly,
     whitney_one_curl_curl_matrix_assembler_test_regular) {
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
      projects::hldo_sphere::assemble::WhitneyOneCurlCurlMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal << 1, 1,-1,
             1, 1,-1,
            -1,-1, 1;
  // clang-format on
  Ae_anal *= 2. / std::sqrt(3);
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
 * @brief Test the element matrix provider for the whitney one form
 *
 * Test the element vector provider on the triangle
 * ((1,0,0), (0,1,0), (0,0,1))
 *
 * with relative edge orientations
 *
 * s_0 = 1, s_1 = 1, s_2 = -1
 *
 * And the linear form
 *
 * @f[
 *  \int\limits_K\ \mathbf{u} grad_{\Gamma}v \,\mathrm{d}x, \quad u \in
 * \bm{H}(curl_{\Gamma}, \partial\mathbb{S})
 * , v \in H^1
 * @f]
 *
 * Using the whitney one forms form u and the barycentric basis function for v
 *
 */
TEST(projects_hldo_sphere_assembly,
     whitney_one_form_grad_matrix_assembler_test_regular) {
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
      projects::hldo_sphere::assemble::WhitneyOneGradMatrixProvider();
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

/**
 *
 * @brief Test the element matrix provider for the whitney one form
 *
 * Given by the expression
 *
 * @f[
 *    \int\limits_K u \ v \ dx
 * @f]
 *
 * Test the element vector provider on the triangle
 * ((1,0,0), (0,1,0), (0,0,1))
 *
 * with relative edge orientations
 *
 * s_0 = 1, s_1 = 1, s_2 = -1
 *
 * And the linear form
 *
 * Using the whitney one forms for u and v
 *
 */
TEST(projects_hldo_sphere_assembly,
     whitney_one_mass_matrix_assembler_test_regular_one) {
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
      projects::hldo_sphere::assemble::WhitneyOneMassMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal << 5,-1, 1,
            -1, 5, 1,
             1, 1, 5;
  // clang-format on
  Ae_anal *= 1. / 12. / std::sqrt(3);
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
 * @brief Test the element vector provider with the one-form
 *
 * Test the element vector provider on the triangle
 * ((1,0,0), (0,1,0), (0,0,1))
 *
 * with relative edge orientations
 *
 * s_0 = 1, s_1 = 1, s_2 = -1
 *
 * And the linear form
 *
 * @f[
 *  \int_K\ f \cdot v \,\mathrm{d}x, \quad f, v \in \bm{H}(div_{\Gamma},
 * \partial\mathbb{S})
 * @f]
 *
 * using the load function
 *
 * @f[
 *  f((x_0,x_1, x_2)) = (1, 1, x_0)
 * @f]
 *
 * and the whithey one form basis functions
 *
 */
TEST(projects_hldo_sphere_assembly, vector_provider_one_form_test_regular) {
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

  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element vec for the triangle
  const auto elem_vec_provider =
      projects::hldo_sphere::assemble::WhitneyOneVectorProvider(f);
  const Eigen::VectorXd Ae = elem_vec_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::VectorXd Ae_anal(3);
  // clang-format off
  Ae_anal <<
    -0.041666666666666664*1/sqrt(3),
    -0.125*sqrt(3),-7/(24.*sqrt(3));
  // clang-format on
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    EXPECT_DOUBLE_EQ(Ae(i), Ae_anal(i)) << "mismatch in entry (" << i << ")";
  }
}

// Additional  tests for the funciton
//
// @f[
// f(\vec{x}) = \begin{pmatrix} -k^2 x_2
// - \frac{2x_2}{\|\vec{x}\|} \\-k^2 x_1 - \frac{2x_1}{\|\vec{x}\|} \\ 0
// \begin{pmatrix}
// @f]
//
// for all triagles on the octacon with radius 1
TEST(projects_hldo_sphere_assembly, vector_provider_one_form_test_octacon) {
  // Build Mesh
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

  projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
      projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

  sphere.setRefinementLevel(0);
  sphere.setRadius(1);

  std::shared_ptr<const lf::mesh::Mesh> mesh = sphere.Build();

  // Define function f
  double k = 0.1;
  auto f = [&](const Eigen::Vector3d x) -> Eigen::Vector3d {
    Eigen::VectorXd res(3);
    res(0) = -k * k * x(1) - 2 * x(1) / x.squaredNorm();
    res(1) = k * k * x(0) + 2 * x(0) / x.squaredNorm();
    res(2) = 0;
    return res / x.norm();
  };

  // Compute the element vec for the triangle
  const auto elem_vec_provider =
      projects::hldo_sphere::assemble::WhitneyOneVectorProvider(f);

  // Get solutions
  // the mesh (octacon) conatains 8 cells
  std::vector<Eigen::VectorXd> Aes(8);
  for (int i = 0; i < 8; i++) {
    const auto element = mesh->EntityByIndex(0, i);
    Aes[i] = elem_vec_provider.Eval(*element);
  }

  // Construct the analytically computed solutions with mathematica
  std::vector<Eigen::Vector3d> Ae_anal(8);
  // clang-format off
  Ae_anal[0] << 1.3642297039667421,-0.2728459407933484,-0.2728459407933484;
  Ae_anal[1] << 0.2728459407933484,0.2728459407933484,1.3642297039667421;
  Ae_anal[2] << 1.3642297039667421,-0.2728459407933484,0.2728459407933484;
  Ae_anal[3] << -0.2728459407933484,0.2728459407933484,1.3642297039667421;
  Ae_anal[4] << 1.3642297039667421,-0.2728459407933484,0.2728459407933484;
  Ae_anal[5] << -0.2728459407933484,0.2728459407933484,1.3642297039667421;
  Ae_anal[6] << 1.3642297039667421,0.2728459407933484,0.2728459407933484;
  Ae_anal[7] << -0.2728459407933484,-0.2728459407933484,1.3642297039667421;
  // clang-format on
  // Assert that the two matrices are approximately equal
  for (int k = 0; k < 8; k++) {
    ASSERT_EQ(Aes[k].rows(), Ae_anal[k].rows());
    ASSERT_EQ(Aes[k].cols(), Ae_anal[k].cols());
    for (int i = 0; i < Ae_anal[k].rows(); ++i) {
      EXPECT_DOUBLE_EQ(Aes[k](i), Ae_anal[k](i))
          << "mismatch in entry (" << i << ") of triangle " << k;
    }
  }
}

/*********************************************************************
 * TESTS ON FLAT TRIANGLE ((3,0.5,0), (0, 1, 0.2), (0, 0, 1))
 *********************************************************************/

/**
 *
 * @brief Test the element matrix provider for the whitney one form
 *
 * Given by the expression
 *
 * @f[
 *    \int\limits_K \text{curl}_{\Gamma} u \ \text{curl}_{\Gamma} v \ dx
 * @f]
 *
 * Test the element vector provider on the triangle
 * ((3,0.5,0), (0, 1, 0.2), (0, 0, 1))
 *
 * with relative edge orientations
 *
 * s_0 = -1, s_1 = 1, s_2 = -1
 *
 * And the linear form
 *
 * Using the whitney one forms for u and u
 *
 */
TEST(projects_hldo_sphere_assembly,
     whitney_one_curl_curl_matrix_assembler_test_flat) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices <<   3,   0,   0,
              0.5,   1,   0,
                0, 0.2,   1;
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

  // build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));

  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element matrix for the triangle
  const auto elem_mat_provider =
      projects::hldo_sphere::assemble::WhitneyOneCurlCurlMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal << 
      0.5143444998736394,-0.5143444998736396,0.5143444998736395,-0.5143444998736396,
   0.5143444998736396,-0.5143444998736396,0.5143444998736395,-0.5143444998736396,0.5143444998736396;
  // clang-format on
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae(i, j), Ae_anal(i, j), 1e-14)
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

/**
 *
 * @brief Test the element matrix provider for the whitney one form
 *
 * Test the element vector provider on the triangle
 * ((3,0.5,0), (0, 1, 0.2), (0, 0, 1))
 *
 * with relative edge orientations
 *
 * s_0 = -1, s_1 = 1, s_2 = -1
 *
 * And the linear form
 *
 * @f[
 *  \int\limits_K\ \mathbf{u} grad_{\Gamma}v \,\mathrm{d}x, \quad u \in
 * \bm{H}(curl_{\Gamma}, \partial\mathbb{S})
 * , v \in H^1
 * @f]
 *
 * Using the whitney one forms form u and the barycentric basis function for v
 *
 */
TEST(projects_hldo_sphere_assembly,
     whitney_one_form_grad_matrix_assembler_test_flat) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices <<   3,   0,   0,
              0.5,   1,   0,
                0, 0.2,   1;
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

  // build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));

  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element matrix for the triangle
  const auto elem_mat_provider =
      projects::hldo_sphere::assemble::WhitneyOneGradMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal <<
    0.1260144024690417,-0.49505658112837825,0.3690421786593366,0.04114755998989116,
   -0.8229511997978233,0.7818036398079323,-0.08486684247915052,-0.32789461866944536,
   0.4127614611485959;
  // clang-format on
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae(i, j), Ae_anal(i, j), 1e-14)
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

/**
 *
 * @brief Test the element matrix provider for the whitney one form
 *
 * Given by the expression
 *
 * @f[
 *    \int\limits_K u \ v \ dx
 * @f]
 *
 * Test the element vector provider on the triangle
 * ((3,0.5,0), (0, 1, 0.2), (0, 0, 1))
 *
 * with relative edge orientations
 *
 * s_0 = -1, s_1 = 1, s_2 = -1
 *
 * And the linear form
 *
 * Using the whitney one forms for u and v
 *
 */
TEST(projects_hldo_sphere_assembly,
     whitney_one_mass_matrix_assembler_test_flat_one) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices <<   3,   0,   0,
              0.5,   1,   0,
                0, 0.2,   1;
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

  // build triangle
  std::unique_ptr<lf::geometry::Geometry> geom =
      std::make_unique<lf::geometry::TriaO1>(vertices);
  const std::array<lf::mesh::MeshFactory::size_type, 3> nodes = {0, 1, 2};
  factory.AddEntity(trig, nonstd::span(nodes.data(), 3), std::move(geom));

  const auto mesh = factory.Build();
  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element matrix for the triangle
  const auto elem_mat_provider =
      projects::hldo_sphere::assemble::WhitneyOneMassMatrixProvider();
  const Eigen::MatrixXd Ae = elem_mat_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::MatrixXd Ae_anal(3, 3);
  // clang-format off
  Ae_anal << 
    0.28267516472222104,0.2123814164061571,0.15666076225317957,0.2123814164061571,
   0.6105697833916662,0.17123385641626587,0.15666076225317957,0.17123385641626587,
   0.24152760473232993;
  // clang-format on
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    for (int j = 0; j < Ae_anal.cols(); ++j) {
      EXPECT_NEAR(Ae(i, j), Ae_anal(i, j), 1e-14)
          << "mismatch in entry (" << i << ", " << j << ")";
    }
  }
}

/**
 *
 * @brief Test the element vector provider with the one-form
 *
 * Test the element vector provider on the triangle
 * ((3,0.5,0), (0.5,1,0.2), (0,0,1))
 *
 * with relative edge orientations
 *
 * s_0 = -1, s_1 = 1, s_2 = -1
 *
 * And the linear form
 *
 * @f[
 *  \int_K\ f \cdot v \,\mathrm{d}x, \quad f, v \in \bm{H}(div_{\Gamma},
 * \partial\mathbb{S})
 * @f]
 *
 * using the load function
 *
 * @f[
 *  f((x_0,x_1, x_2)) = (1, 1, x_0)
 * @f]
 *
 * and the whithey one form basis functions
 *
 */
TEST(projects_hldo_sphere_assembly, vector_provider_one_form_test_flat) {
  // Build Mesh
  const auto trig = lf::base::RefEl::kTria();
  const auto seg = lf::base::RefEl::kSegment();

  // Define Vertices
  Eigen::MatrixXd vertices(3, 3);
  // clang-format off
  vertices <<   3,   0, 0,
              0.5,   1, 0,
                0, 0.2, 1;
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
    res(0) = cos(x(0) - x(2));
    res(1) = sin(x(0) + x(1));
    res(2) = x(2);
    return res;
  };

  const auto element = mesh->EntityByIndex(0, 0);
  // Compute the element vec for the triangle
  const auto elem_vec_provider =
      projects::hldo_sphere::assemble::WhitneyOneVectorProvider(f);
  const Eigen::VectorXd Ae = elem_vec_provider.Eval(*element);
  // Construct the analytically computed element matrix
  Eigen::VectorXd Ae_anal(3);
  // clang-format off
  Ae_anal <<
        -0.13142416550290498,-0.17929962791878098,
        -0.3782779903706242;
  // clang-format on
  // Assert that the two matrices are approximately equal
  ASSERT_EQ(Ae.rows(), Ae_anal.rows());
  ASSERT_EQ(Ae.cols(), Ae_anal.cols());
  for (int i = 0; i < Ae_anal.rows(); ++i) {
    EXPECT_DOUBLE_EQ(Ae(i), Ae_anal(i)) << "mismatch in entry (" << i << ")";
  }
}
