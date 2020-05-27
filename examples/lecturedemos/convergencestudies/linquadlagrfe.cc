/**
 * @file linquadlagrfe.cc
 * @brief Creates convergence plots for experiment 3.2.3.7 and 3.2.3.8
 * @author Tobias Rohner
 * @date April 2020
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/fe/fe.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

namespace po = boost::program_options;

/**
 * @brief Solves the Poisson problem on a given mesh
 * @param mesh The mesh on which to solve the PDE
 * @param fe_space The Finite Element Space to use
 * @returns The dof vector corresponding to the solution of the PDE
 *
 * The Dirichlet boundary conditions are set to zero and the load is given by
 * \f[
        f(x) = 2\pi^2\sin(\pi x_1)\sin\pi x_2)
   \f]
 */
Eigen::VectorXd solvePoisson(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>>
        &fe_space) {
  // Define the load function of the manufactured solution
  const auto load = [](const Eigen::Vector2d &x) -> double {
    return 2 * M_PI * M_PI * std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
  };
  const lf::mesh::utils::MeshFunctionGlobal mf_load(load);

  // Intialize the matrix and vector providers
  const lf::mesh::utils::MeshFunctionConstant<double> mf_alpha(1);
  const lf::mesh::utils::MeshFunctionConstant<double> mf_gamma(0);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider element_matrix_provider(
      fe_space, mf_alpha, mf_gamma);
  lf::uscalfe::ScalarLoadElementVectorProvider element_vector_provider(
      fe_space, mf_load,
      {{lf::base::RefEl::kTria(),
        lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 6)},
       {lf::base::RefEl::kQuad(),
        lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 6)}});

  // Assemble the system matrix and right hand side
  const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
  lf::assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());
  std::cout << "\t\t> Assembling System Matrix" << std::endl;
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, element_matrix_provider,
                                      A_COO);
  std::cout << "\t\t> Assembling right Hand Side" << std::endl;
  lf::assemble::AssembleVectorLocally(0, dofh, element_vector_provider, rhs);

  // Enforce zero dirichlet boundary conditions
  std::cout << "\t\t> Enforcing Boundary Conditions" << std::endl;
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);
  const auto selector = [&](unsigned int idx) -> std::pair<bool, double> {
    return {boundary(dofh.Entity(idx)), 0};
  };
  lf::assemble::FixFlaggedSolutionComponents(selector, A_COO, rhs);

  // Solve the LSE using the cholesky decomposition
  std::cout << "\t\t> Solving LSE" << std::endl;
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
  const Eigen::VectorXd solution = solver.solve(rhs);

  // Return the resulting dof vector
  return solution;
}

int main(int argc, char *argv[]) {
  const int num_meshes = 7;

  po::options_description desc("allowed options");
  desc.add_options()("output,o", po::value<std::string>(),
                     "Name of the output file");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  if (vm.count("output") == 0) {
    std::cout << desc << std::endl;
    exit(1);
  }
  const std::string output_file = vm["output"].as<std::string>();

  // The analytic solution
  const auto u = [](const Eigen::VectorXd &x) -> double {
    return std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
  };
  const lf::mesh::utils::MeshFunctionGlobal mf_u(u);
  // The gradient of the analytic solution
  const auto u_grad = [](const Eigen::Vector2d &x) -> Eigen::Vector2d {
    Eigen::Vector2d grad;
    grad[0] = M_PI * std::cos(M_PI * x[0]) * std::sin(M_PI * x[1]);
    grad[1] = M_PI * std::sin(M_PI * x[0]) * std::cos(M_PI * x[1]);
    return grad;
  };
  const lf::mesh::utils::MeshFunctionGlobal mf_u_grad(u_grad);

  const boost::filesystem::path here = __FILE__;
  const boost::filesystem::path mesh_folder = here.parent_path() / "meshes";
  Eigen::MatrixXd results(num_meshes, 7);
  for (int mesh_idx = 0; mesh_idx < num_meshes; ++mesh_idx) {
    std::cout << "> Mesh Nr. " << mesh_idx << std::endl;

    // Load the mesh
    std::cout << "\t> Loading Mesh" << std::endl;
    const std::string mesh_name =
        "unitsquare" + std::to_string(mesh_idx) + ".msh";
    const boost::filesystem::path mesh_file = mesh_folder / mesh_name;
    auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    const lf::io::GmshReader reader(std::move(factory), mesh_file.string());
    const auto mesh = reader.mesh();

    // Solve the problem with linear finite elements
    std::cout << "\t> Linear Lagrangian FE";
    const auto fe_space_o1 =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
    std::cout << " (" << fe_space_o1->LocGlobMap().NumDofs() << " DOFs)"
              << std::endl;
    const Eigen::VectorXd solution_o1 = solvePoisson(mesh, fe_space_o1);
    const lf::fe::MeshFunctionFE<double, double> mf_o1(fe_space_o1,
                                                            solution_o1);
    const lf::fe::MeshFunctionGradFE<double, double> mf_grad_o1(
        fe_space_o1, solution_o1);

    // Solve the problem with quadratic finite elements
    std::cout << "\t> Quadratic Lagrangian FE";
    const auto fe_space_o2 =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh);
    std::cout << " (" << fe_space_o2->LocGlobMap().NumDofs() << " DOFs)"
              << std::endl;
    const Eigen::VectorXd solution_o2 = solvePoisson(mesh, fe_space_o2);
    const lf::fe::MeshFunctionFE<double, double> mf_o2(fe_space_o2,
                                                            solution_o2);
    const lf::fe::MeshFunctionGradFE<double, double> mf_grad_o2(
        fe_space_o2, solution_o2);

    // Compute the H1 and L2 errors
    std::cout << "\t> Computing Error Norms" << std::endl;
    const auto quadrule_provider = [](const lf::mesh::Entity &entity) {
      return lf::quad::make_QuadRule(entity.RefEl(), 6);
    };
    const double H1_err_o1 = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh, lf::mesh::utils::squaredNorm(mf_grad_o1 - mf_u_grad),
        quadrule_provider));
    const double H1_err_o2 = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh, lf::mesh::utils::squaredNorm(mf_grad_o2 - mf_u_grad),
        quadrule_provider));
    const double L2_err_o1 = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh, lf::mesh::utils::squaredNorm(mf_o1 - mf_u), quadrule_provider));
    const double L2_err_o2 = std::sqrt(lf::fe::IntegrateMeshFunction(
        *mesh, lf::mesh::utils::squaredNorm(mf_o2 - mf_u), quadrule_provider));

    // Store the mesh width, the number of DOFs and the errors in the results
    // matrix
    results(mesh_idx, 0) = std::sqrt(1. / mesh->NumEntities(0));
    results(mesh_idx, 1) = fe_space_o1->LocGlobMap().NumDofs();
    results(mesh_idx, 2) = fe_space_o2->LocGlobMap().NumDofs();
    results(mesh_idx, 3) = H1_err_o1;
    results(mesh_idx, 4) = H1_err_o2;
    results(mesh_idx, 5) = L2_err_o1;
    results(mesh_idx, 6) = L2_err_o2;
  }

  // Output the resulting errors to a file
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream file;
  file.open(output_file);
  file << results.format(CSVFormat);
  file.close();

  return 0;
}
