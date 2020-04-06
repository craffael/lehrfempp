/**
 * @file linquadlagrfe.cc
 * @brief Creates convergence plots for experiment 3.2.3.7
 * @author Tobias Rohner
 * @date April 2020
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <lf/base/base.h>
#include <lf/uscalfe/uscalfe.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/assemble/assemble.h>
#include <lf/quad/quad.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/io/io.h>
#include <boost/filesystem.hpp>



Eigen::VectorXd solvePoisson(const std::shared_ptr<const lf::mesh::Mesh> &mesh, const std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> &fe_space) {
    // Define the load function of the manufactured solution
    const auto load = [](const Eigen::Vector2d &x) -> double {
	return 2*M_PI*M_PI * std::sin(M_PI*x[0]) * std::sin(M_PI*x[1]);
    };
    const lf::mesh::utils::MeshFunctionGlobal mf_load(load);

    // Intialize the matrix and vector providers
    const lf::mesh::utils::MeshFunctionConstant<double> mf_alpha(1);
    const lf::mesh::utils::MeshFunctionConstant<double> mf_gamma(0);
    lf::uscalfe::ReactionDiffusionElementMatrixProvider element_matrix_provider(fe_space, mf_alpha, mf_gamma);
    lf::uscalfe::ScalarLoadElementVectorProvider element_vector_provider(fe_space, mf_load, {{lf::base::RefEl::kTria(), lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 6)}, {lf::base::RefEl::kQuad(), lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 6)}});

    // Assemble the system matrix and right hand side
    const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
    lf::assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());
    std::cout << "\t\t> Assembling System Matrix" << std::endl;
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, element_matrix_provider, A_COO);
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

    const boost::filesystem::path here = __FILE__;
    const boost::filesystem::path mesh_folder = here.parent_path() / "meshes";
    for (int mesh_idx = 0 ; mesh_idx < num_meshes ; ++mesh_idx) {
	std::cout << "> Mesh Nr. " << mesh_idx << std::endl;

	// Load the mesh
	std::cout << "\t> Loading Mesh" << std::endl;
	const std::string mesh_name = "unitsquare" + std::to_string(mesh_idx) + ".msh";
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
	const lf::uscalfe::MeshFunctionGradFE<double, double> mf_grad_o2(fe_space_o1, solution_o2);
    }

    return 0;
}
