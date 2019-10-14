/**
 * @file
 * @brief We integrate a function over a square using quadrature rules of a
 * fixed degree and refine the mesh to achieve convergence.
 * @author Raffael Casagrande
 * @date   2018-09-02 05:17:45
 * @copyright MIT License
 */

#include <boost/math/constants/constants.hpp>
#include <boost/program_options.hpp>

#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/quad/quad.h>

#include "lf/refinement/test/refinement_test_utils.h"

template <class F>
double integrate(const lf::mesh::Mesh& mesh, lf::quad::quadDegree_t degree,
                 F f) {
  double result = 0.;
  auto qr_tria = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), degree);
  auto qr_quad = lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), degree);

  Eigen::MatrixXd points(1, 1);
  Eigen::VectorXd weights(1);

  for (auto e : mesh.Entities(0)) {
    if (e->RefEl() == lf::base::RefEl::kTria()) {
      points = qr_tria.Points();
      weights = qr_tria.Weights();
    } else {
      points = qr_quad.Points();
      weights = qr_quad.Weights();
    }
    auto mapped_points = e->Geometry()->Global(points);
    auto integration_elements = e->Geometry()->IntegrationElement(points);
    for (Eigen::Index j = 0; j < points.cols(); ++j) {
      result += f(mapped_points.col(j)) * weights(j) * integration_elements(j);
    }
  }
  return result;
}

int main(int argc, char** argv) {
  // define allowed command line arguments:
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()
  ("help", "produce this help message")
  ("quad_degree", po::value<int>()->default_value(3), "The degree of the local quadrature rule.")
  ("max_level", po::value<int>()->default_value(5), "The number of refinement levels")
  ;
  // clang-format on
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help") != 0U) {
    std::cout << desc << std::endl;
    return 1;
  }

  int max_level = vm["max_level"].as<int>();
  int quad_degree = vm["quad_degree"].as<int>();

  // Create a three element mesh that contains two triangles and one
  // quadrilateral
  lf::mesh::hybrid2d::MeshFactory mesh_factory(2);
  mesh_factory.AddPoint(Eigen::Vector2d{0, 0});
  mesh_factory.AddPoint(Eigen::Vector2d{0.5, 0});
  mesh_factory.AddPoint(Eigen::Vector2d{1, 0});
  mesh_factory.AddPoint(Eigen::Vector2d{1, 1});
  mesh_factory.AddPoint(Eigen::Vector2d{0.5, 1});
  mesh_factory.AddPoint(Eigen::Vector2d{0, 1});
  Eigen::MatrixXd node_coords(2, 3);
  node_coords << 0, 0.5, 0.5, 0, 0, 1;
  mesh_factory.AddEntity(lf::base::RefEl::kTria(),
                         std::vector<lf::base::size_type>{0, 1, 4},
                         std::make_unique<lf::geometry::TriaO1>(node_coords));
  node_coords << 0, 0.5, 0, 0, 1, 1;
  mesh_factory.AddEntity(lf::base::RefEl::kTria(),
                         std::vector<lf::base::size_type>{0, 4, 5},
                         std::make_unique<lf::geometry::TriaO1>(node_coords));
  node_coords = Eigen::MatrixXd(2, 4);
  node_coords << 0.5, 1, 1, 0.5, 0, 0, 1, 1;
  mesh_factory.AddEntity(lf::base::RefEl::kQuad(),
                         std::vector<lf::base::size_type>{1, 2, 3, 4},
                         std::make_unique<lf::geometry::QuadO1>(node_coords));

  auto base_mesh = mesh_factory.Build();

  // parameters:
  auto pi = boost::math::constants::pi<double>();
  auto f = [&](const Eigen::Vector2d& x) {
    return std::sin(pi * x(0)) * std::pow(std::cos(pi * x(1)), 2);
  };
  auto exact_integral = 1. / pi;

  auto mesh = base_mesh;
  auto errors = Eigen::VectorXd(max_level + 1);
  for (int level = 0; level <= max_level; ++level) {
    lf::io::VtkWriter vtk_writer(mesh,
                                 "level" + std::to_string(level) + ".vtk");

    auto approx = integrate(*mesh, quad_degree, f);
    errors(level) = std::abs(approx - exact_integral);

    lf::refinement::MeshHierarchy mh(
        mesh, std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2));
    mh.RefineRegular();
    mesh = mh.getMesh(1);
  }

  std::cout << "measured errors for each level: " << std::endl;
  std::cout << errors << std::endl << std::endl;

  // estimate the rate of convergence:
  Eigen::MatrixXd A(max_level + 1, 2);
  A.col(0).setOnes();
  A.col(1) = (-Eigen::ArrayXd::LinSpaced(max_level + 1, 1, max_level + 1) *
              std::log(2))
                 .matrix();

  Eigen::VectorXd b = errors.transpose().array().log().matrix();

  // consider only the three smallest meshes:
  A = A.bottomRows(3);
  b = b.bottomRows(3);

  Eigen::VectorXd x =
      A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
  std::cout << "estimated order of convergence = " << x(1) << std::endl;
}
