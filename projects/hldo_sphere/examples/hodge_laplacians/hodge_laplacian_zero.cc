/**
 * @file hodge_laplacians.cc
 * @brief Solve hodge laplacian source problem for fixed function u = 0 and
 * fixed k = 100
 */

#include <hodge_laplacians_source_problems.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/tp_triag_mesh_builder.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>
#include <mesh_function_whitney_one.h>
#include <mesh_function_whitney_two.h>
#include <mesh_function_whitney_zero.h>
#include <norms.h>
#include <sphere_triag_mesh_builder.h>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using lf::uscalfe::operator-;

/**
 * @brief Concatenate objects defining an operator<<(std::ostream&)
 * @param args A variadic pack of objects implementing
 * `operator<<(std::ostream&)`
 * @returns A string with the objects concatenated
 */
template <typename... Args>
static std::string concat(Args &&...args) {
  std::ostringstream ss;
  (ss << ... << args);
  return ss.str();
}

/**
 * @brief stores information to recover convergence properties
 */
struct ProblemSolution {
  std::shared_ptr<const lf::mesh::Mesh> mesh;
  Eigen::VectorXd mu_zero;
  Eigen::VectorXd mu_one;
  Eigen::VectorXd mu_two;
};

/**
 * @brief Prints the L2 norm errors of the experiment
 */
int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " max_refinement_level" << std::endl;
    exit(1);
  }

  const double k = 100;
  const double r = 2;
  const unsigned refinement_level = atoi(argv[1]);
  std::cout << "max_refinement_level : " << refinement_level << std::endl;

  // Read the mesh from the gmsh file
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

  projects::hldo_sphere::mesh::SphereTriagMeshBuilder sphere =
      projects::hldo_sphere::mesh::SphereTriagMeshBuilder(std::move(factory));

  // we always take the same radius
  sphere.setRadius(1);

  // mathematica function output requries the following helpers
  auto Power = [](double a, double b) -> double { return std::pow(a, b); };
  auto Sin = [](double a) -> double { return std::sin(a); };
  auto Cos = [](double a) -> double { return std::cos(a); };

  // righthandside for the zero and two form
  auto f_zero = [&](const Eigen::Vector3d &x_vec) -> double { return 0; };

  // righthandside for the one form
  auto f_one = [&](const Eigen::Vector3d &x_vec) -> Eigen::VectorXd {
    Eigen::VectorXd ret(3);
    ret << 0, 0, 0;
    return ret;
  };

  // Solve the problem for each mesh in the hierarchy
  std::vector<ProblemSolution> solutions(refinement_level + 1);

  projects::hldo_sphere::discretization::HodgeLaplaciansSourceProblems
      lse_builder;
  lse_builder.SetLoadFunctions(f_zero, f_one, f_zero);

  for (lf::base::size_type lvl = 0; lvl < refinement_level + 1; ++lvl) {
    // get mesh
    sphere.setRefinementLevel(lvl);
    const std::shared_ptr<lf::mesh::Mesh> mesh = sphere.Build();

    auto &sol = solutions[lvl];
    sol.mesh = mesh;

    // setup the problem
    lse_builder.SetMesh(mesh);
    lse_builder.SetK(k);

    // compute the system
    lse_builder.Compute();

    // solve the system
    lse_builder.Solve();

    // store solutions
    sol.mu_zero = lse_builder.GetMuZero();
    sol.mu_one = std::get<0>(lse_builder.GetMuOne());
    sol.mu_two = std::get<1>(lse_builder.GetMuTwo());
  }

  // Compute the analytic solution of the problem
  auto u_zero = [&](const Eigen::Vector3d &x_vec) -> double { return 0; };

  // Compute the analytic solution of the problem
  auto u_one = [&](const Eigen::Vector3d &x_vec) -> Eigen::Vector3d {
    Eigen::VectorXd ret(3);

    ret << 0, 0, 0;
    return ret;
  };

  for (lf::base::size_type lvl = 0; lvl <= refinement_level; ++lvl) {
    // create mesh Functions for solutions
    auto &sol = solutions[lvl];
    projects::hldo_sphere::post_processing::MeshFunctionWhitneyZero mf_zero(
        sol.mu_zero, sol.mesh);
    projects::hldo_sphere::post_processing::MeshFunctionWhitneyOne mf_one(
        sol.mu_one, sol.mesh);
    projects::hldo_sphere::post_processing::MeshFunctionWhitneyTwo mf_two(
        sol.mu_two, sol.mesh);

    lf::mesh::utils::MeshFunctionGlobal<decltype(u_zero)> mf_zero_ana(u_zero);

    lf::mesh::utils::MeshFunctionGlobal<decltype(u_one)> mf_one_ana(u_one);

    // Compute the error of the solutions
    auto mf_diff_zero = mf_zero - mf_zero_ana;
    auto mf_diff_one = mf_one - mf_one_ana;
    auto mf_diff_two = mf_two - mf_zero_ana;

    auto square_scalar = [](double a) -> double { return a * a; };
    auto square_vector = [](Eigen::VectorXd a) -> double {
      return a.squaredNorm();
    };

    lf::quad::QuadRule qr =
        lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 6);

    // Perform post processing on the data
    lf::io::VtkWriter writer(sol.mesh, concat("result_zero_", lvl, ".vtk"));

    // get error for zero form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> L2_zero =
        projects::hldo_sphere::post_processing::L2norm(sol.mesh, mf_diff_zero,
                                                       square_scalar, qr);
    const double error_zero = std::get<0>(L2_zero);
    const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_zero =
        std::get<1>(L2_zero);
    writer.WriteCellData(concat("mf_zero_diff_", lvl), data_set_error_zero);

    // get error for one form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> L2_one =
        projects::hldo_sphere::post_processing::L2norm(sol.mesh, mf_diff_one,
                                                       square_vector, qr);
    const double error_one = std::get<0>(L2_one);
    const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_one =
        std::get<1>(L2_one);
    writer.WriteCellData(concat("mf_one_diff_", lvl), data_set_error_one);

    // get error for two form
    const std::pair<double, lf::mesh::utils::CodimMeshDataSet<double>> L2_two =
        projects::hldo_sphere::post_processing::L2norm(sol.mesh, mf_diff_two,
                                                       square_scalar, qr);
    const double error_two = std::get<0>(L2_two);
    const lf::mesh::utils::CodimMeshDataSet<double> data_set_error_two =
        std::get<1>(L2_two);
    writer.WriteCellData(concat("mf_two_diff_", lvl), data_set_error_two);

    // print all the errors
    std::cout << lvl << ' ' << solutions[lvl].mesh->NumEntities(2) << ' '
              << error_zero << ' ' << error_one << ' ' << error_two
              << std::endl;
  }

  return 0;
}
