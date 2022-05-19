#include <gtest/gtest.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <norms.h>
#include <sphere_triag_mesh_builder.h>

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
  const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> out_zero =
      projects::hldo_sphere::post_processing::L2norm<decltype(mf_f_zero),
                                                     decltype(sq_f)>(
          mesh, mf_f_zero, sq_f, qr);
  double L2_zero = std::get<0>(out_zero);
  ASSERT_NEAR(L2_zero, 0, 1e-10);

  // create function f
  auto f_const = [](Eigen::Vector3d x) -> double { return 5.0; };
  lf::mesh::utils::MeshFunctionGlobal<decltype(f_const)> mf_f_const(f_const);
  // Compute the L2 norm of the function constant zero on the mesh
  const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> out_const =
      projects::hldo_sphere::post_processing::L2norm<decltype(mf_f_const),
                                                     decltype(sq_f)>(
          mesh, mf_f_const, sq_f, qr);
  double L2_const = std::sqrt(std::get<0>(out_const));
  ASSERT_NEAR(L2_const, sqrt(sqrt(0.75) * 2 * 25.0), 1e-12);
}

/**
 *
 * Tests the L2_norm on the Sphere with the function
 *
 * f(x) = (I - x x^T / x^Tx) (x(1)^2, x(2)^2, x(0)^2)
 *
 * For which the analytical squared integral over the sphere is
 * 72 \pi / 25 = 6.4627...
 *
 */
TEST(projects_hldo_sphere_post_processing, norm_l2_vector_field) {
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

  projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
      projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

  sphere.setRadius(1);

  sphere.setRefinementLevel(6);

  const std::shared_ptr<lf::mesh::Mesh> mesh = sphere.Build();

  // create function f
  auto Power = [](double a, double b) -> double { return pow(a, b); };
  auto Abs = [](double a) -> double { return abs(a); };
  auto f = [&](Eigen::Vector3d x_) -> Eigen::VectorXd {
    // scale x_ up to the sphere
    x_ = x_ / x_.norm();

    double x = x_(0);
    double y = x_(1);
    double z = x_(2);

    Eigen::VectorXd res(3);
    res << -((Power(x, 3) * z) /
             (Power(Abs(x), 2) + Power(Abs(y), 2) + Power(Abs(z), 2))) -
               (x * y * Power(z, 2)) /
                   (Power(Abs(x), 2) + Power(Abs(y), 2) + Power(Abs(z), 2)) +
               Power(y, 2) *
                   (1 - Power(x, 2) / (Power(Abs(x), 2) + Power(Abs(y), 2) +
                                       Power(Abs(z), 2))),
        -((x * Power(y, 3)) /
          (Power(Abs(x), 2) + Power(Abs(y), 2) + Power(Abs(z), 2))) -
            (Power(x, 2) * y * z) /
                (Power(Abs(x), 2) + Power(Abs(y), 2) + Power(Abs(z), 2)) +
            Power(z, 2) *
                (1 - Power(y, 2) / (Power(Abs(x), 2) + Power(Abs(y), 2) +
                                    Power(Abs(z), 2))),
        -((x * Power(y, 2) * z) /
          (Power(Abs(x), 2) + Power(Abs(y), 2) + Power(Abs(z), 2))) -
            (y * Power(z, 3)) /
                (Power(Abs(x), 2) + Power(Abs(y), 2) + Power(Abs(z), 2)) +
            Power(x, 2) *
                (1 - Power(z, 2) / (Power(Abs(x), 2) + Power(Abs(y), 2) +
                                    Power(Abs(z), 2)));
    return res;
  };
  lf::mesh::utils::MeshFunctionGlobal<decltype(f)> mf_f(f);

  // create the square funciton of scalar values
  auto sq_f = [](Eigen::Vector3d x) -> double { return x.squaredNorm(); };

  // create the quadrature rule which is exact for cosnt functions
  lf::quad::QuadRule qr = lf::quad::make_TriaQR_EdgeMidpointRule();

  // Compute the L2 norm of the function constant zero on the mesh
  const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> out_norm =
      projects::hldo_sphere::post_processing::L2norm<decltype(mf_f),
                                                     decltype(sq_f)>(mesh, mf_f,
                                                                     sq_f, qr);
  const double L2 = std::sqrt(std::get<0>(out_norm));
  ASSERT_NEAR(L2, 2.54218506, 1e-3);
}
