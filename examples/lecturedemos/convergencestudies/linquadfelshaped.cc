/**
 * @file linquadfelshaped.cc
 * @brief Creates convergence plots for experiment 3.2.3.10
 * @author Tobias Rohner
 * @date April 2020
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <memory>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <lf/mesh/mesh.h>
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>


namespace po = boost::program_options;



Eigen::VectorXd solvePoisson(const std::shared_ptr<const lf::mesh::Mesh> &mesh, const std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> &fe_space) {
    // Define the boundary values
    const auto u_bd = [](const Eigen::Vector2d &x) -> double {
	const double r = x.norm();
	double phi = std::atan2(x[1], x[0]);
	if (phi < 0) {
	    phi += 2*M_PI;
	}
	return std::pow(r, 2./3) * std::sin(2./3 * phi);
    };

    // Initialize the matrix provider
    const lf::mesh::utils::MeshFunctionConstant<double> mf_alpha(1);
    const lf::mesh::utils::MeshFunctionConstant<double> mf_gamma(0);
    lf::uscalfe::ReactionDiffusionElementMatrixProvider element_matrix_provider(fe_space, mf_alpha, mf_gamma);

    // Assemble the system matrix (RHS is zero because we have no load)
    const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
    lf::assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());
    std::cout << "\t\t> Assembling System Matrix" << std::endl;
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, element_matrix_provider, A_COO);

    // Enforce the dirichlet boundary conditions
    std::cout << "\t\t> Enforcing Boundary Conditions" << std::endl;
    const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);
    const auto selector = [&](unsigned int idx) -> std::pair<bool, double> {
	const lf::mesh::Entity &entity = dofh.Entity(idx);
	if (!boundary(entity)) {
	    return {false, 0};
	}
	const lf::geometry::Geometry *geom = entity.Geometry();
	// Find out where the evaluation node corresponding to the DOF is
	const auto shape_function_layout = fe_space->ShapeFunctionLayout(entity.RefEl());
	const auto num_int_dofs = dofh.NumInteriorDofs(entity);
	const auto int_glob_dof_idxs = dofh.InteriorGlobalDofIndices(entity);
	int int_dof_idx;
	for (int_dof_idx = 0 ; int_dof_idx < num_int_dofs ; ++int_dof_idx) {
	    if (int_glob_dof_idxs[int_dof_idx] == idx) {
		break;
	    }
	}
	const Eigen::VectorXd eval_node = shape_function_layout->EvaluationNodes().col(int_dof_idx);
	const Eigen::Vector2d pos = geom->Global(eval_node);
	// Return the value on the boundary at the position of the evaluation node corresponding to the dof
	return {true, u_bd(pos)};
    };
    lf::assemble::FixFlaggedSolutionComponents(selector, A_COO, rhs);

    // Solve the LSE using Cholesky decomposition
    std::cout << "\t\t> Solving LSE" << std::endl;
    Eigen::SparseMatrix<double> A = A_COO.makeSparse();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
    const Eigen::VectorXd solution = solver.solve(rhs);

    // Return the resulting solution vector
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
    const auto u = [](const Eigen::Vector2d &x) -> double {
	const double r = x.norm();
	double phi = std::atan2(x[1], x[0]);
	if (phi < 0) {
	    phi += 2*M_PI;
	}
	return std::pow(r, 2./3) * std::sin(2./3 * phi);
    };
    lf::mesh::utils::MeshFunctionGlobal mf_u(u);
    // The gradient of the analytic solution
    const auto u_grad = [](const Eigen::Vector2d &x) -> Eigen::Vector2d {
	const double r = x.norm();
	double phi = std::atan2(x[1], x[0]);
	if (phi < 0) {
	    phi += 2*M_PI;
	}
	Eigen::Vector2d grad;
	grad[0] = 2./3 * std::pow(r, -4./3) * (x[0]*std::sin(2./3*phi) - x[1]*std::cos(2./3*phi));
	grad[1] = 2./3 * std::pow(r, -4./3) * (x[1]*std::sin(2./3*phi) + x[0]*std::cos(2./3*phi));
	return grad;
    };
    lf::mesh::utils::MeshFunctionGlobal mf_u_grad(u_grad);

    const boost::filesystem::path here = __FILE__;
    const boost::filesystem::path mesh_folder = here.parent_path() / "meshes";
    Eigen::MatrixXd results(num_meshes, 5);
    for (int mesh_idx = 0 ; mesh_idx < num_meshes ; ++mesh_idx) {
	std::cout << "> Mesh Nr. " << mesh_idx << std::endl;

	// Load the mesh
	std::cout << "\t> Loading Mesh" << std::endl;
	const std::string mesh_name = "L" + std::to_string(mesh_idx) + ".msh";
	const boost::filesystem::path mesh_file = mesh_folder / mesh_name;
	auto factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
	const lf::io::GmshReader reader(std::move(factory), mesh_file.string());
	const auto mesh = reader.mesh();

	// Solve the problem with linear finite elements
	std::cout << "\t> Linear Lagrangian FE";
	const auto fe_space_o1 = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
	std::cout << " (" << fe_space_o1->LocGlobMap().NumDofs() << " DOFs)" << std::endl;
	const Eigen::VectorXd solution_o1 = solvePoisson(mesh, fe_space_o1);
	const lf::uscalfe::MeshFunctionGradFE<double, double> mf_grad_o1(fe_space_o1, solution_o1);

	// Solve the problem with quadratic finite elements
	std::cout << "\t> Quadratic Lagrangian FE";
	const auto fe_space_o2 = std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(mesh);
	std::cout << " (" << fe_space_o2->LocGlobMap().NumDofs() << " DOFs)" << std::endl;
	const Eigen::VectorXd solution_o2 = solvePoisson(mesh, fe_space_o2);
	const lf::uscalfe::MeshFunctionGradFE<double, double> mf_grad_o2(fe_space_o2, solution_o2);

	// Compute the errors
	std::cout << "\t> Computing Error Norms" << std::endl;
	const auto quadrule_provider = [](const lf::mesh::Entity &entity) {
	    return lf::quad::make_QuadRule(entity.RefEl(), 10);
	};
	const double H1_err_o1 = std::sqrt(lf::uscalfe::IntegrateMeshFunction(*mesh, lf::mesh::utils::squaredNorm(mf_grad_o1 - mf_u_grad), quadrule_provider));
	const double H1_err_o2 = std::sqrt(lf::uscalfe::IntegrateMeshFunction(*mesh, lf::mesh::utils::squaredNorm(mf_grad_o2 - mf_u_grad), quadrule_provider));

	// Store the mesh width, the number of DOFs and the errors
	results(mesh_idx, 0) = std::sqrt(1. / mesh->NumEntities(0));
	results(mesh_idx, 1) = fe_space_o1->LocGlobMap().NumDofs();
	results(mesh_idx, 2) = fe_space_o2->LocGlobMap().NumDofs();
	results(mesh_idx, 3) = H1_err_o1;
	results(mesh_idx, 4) = H1_err_o2;
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
