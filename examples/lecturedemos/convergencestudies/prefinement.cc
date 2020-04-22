/**
 * @file prefinement.cc	
 * @brief Creates convergence plots for experiment 3.2.3.11
 * @author Tobias Rohner
 * @date April 2020
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES

#include "fespacelagrangeon.h"
#include <iostream>
#include <lf/uscalfe/uscalfe.h>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <string>
#include <cmath>
#include <tuple>
#include <lf/mesh/mesh.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <cstdlib>
#include <vector>



namespace po = boost::program_options;



std::shared_ptr<lf::mesh::Mesh> getSquareDomain() {
    lf::mesh::hybrid2d::MeshFactory factory(2);
    // Add the vertices
    std::vector<lf::mesh::MeshFactory::size_type> vertices;
    Eigen::Vector2d vertex_coord;
    vertex_coord << 0, 0;
    vertices.push_back(factory.AddPoint(vertex_coord));
    vertex_coord << 1, 0;
    vertices.push_back(factory.AddPoint(vertex_coord));
    vertex_coord << 1, 1;
    vertices.push_back(factory.AddPoint(vertex_coord));
    vertex_coord << 0, 1;
    vertices.push_back(factory.AddPoint(vertex_coord));
    // Add the triangles
    Eigen::Matrix<double, Eigen::Dynamic, 3> coords(2, 3);
    lf::mesh::MeshFactory::size_type nodes[3];
    coords << 0, 1, 0,
	      0, 0, 1;
    nodes[0] = vertices[0];
    nodes[1] = vertices[1];
    nodes[2] = vertices[3];
    auto geom_tria1 = std::make_unique<lf::geometry::TriaO1>(coords);
    factory.AddEntity(lf::base::RefEl::kTria(), nodes, std::move(geom_tria1));
    coords << 1, 1, 0,
	      0, 1, 1;
    nodes[0] = vertices[1];
    nodes[1] = vertices[2];
    nodes[2] = vertices[3];
    auto geom_tria2 = std::make_unique<lf::geometry::TriaO1>(coords);
    factory.AddEntity(lf::base::RefEl::kTria(), nodes, std::move(geom_tria2));
    // Build the mesh
    return factory.Build();

}



std::tuple<double, double> computeErrorsSquareDomain(const std::shared_ptr<const lf::mesh::Mesh> &mesh, const std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> & fe_space) {
    // Get the degree of the fe space
    const auto degree = fe_space->ShapeFunctionLayout(lf::base::RefEl::kTria())->Degree();
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

    // Define the load function of the manufactured solution
    const auto load = [](const Eigen::Vector2d &x) -> double {
      return 2 * M_PI * M_PI * std::sin(M_PI * x[0]) * std::sin(M_PI * x[1]);
    };
    const lf::mesh::utils::MeshFunctionGlobal mf_load(load);
  
    // Assemble the system matrix and right hand side
    const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
    lf::assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());
    std::cout << "\t\t> Assembling System Matrix" << std::endl;
    const lf::mesh::utils::MeshFunctionConstant<double> mf_alpha(1);
    const lf::mesh::utils::MeshFunctionConstant<double> mf_gamma(0);
    lf::uscalfe::ReactionDiffusionElementMatrixProvider element_matrix_provider(
        fe_space, mf_alpha, mf_gamma);
    lf::assemble::AssembleMatrixLocally(0, dofh, dofh, element_matrix_provider,
                                        A_COO);
    std::cout << "\t\t> Assembling right Hand Side" << std::endl;
    lf::uscalfe::ScalarLoadElementVectorProvider element_vector_provider(
        fe_space, mf_load,
        {{lf::base::RefEl::kTria(),
          lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2*degree-1)},
         {lf::base::RefEl::kQuad(),
          lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 2*degree-1)}});
    lf::assemble::AssembleVectorLocally(0, dofh, element_vector_provider, rhs);
  
    // Enforce zero dirichlet boundary conditions
    std::cout << "\t\t> Enforcing Boundary Conditions" << std::endl;
    const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);
    const auto selector = [&](unsigned int idx) -> std::pair<bool, double> {
      const auto& entity = dofh.Entity(idx);
      return {entity.Codim() > 0 && boundary(entity), 0};
    };
    lf::assemble::FixFlaggedSolutionComponents(selector, A_COO, rhs);
  
    // Solve the LSE using the cholesky decomposition
    std::cout << "\t\t> Solving LSE" << std::endl;
    Eigen::SparseMatrix<double> A = A_COO.makeSparse();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
    const Eigen::VectorXd solution = solver.solve(rhs);
    const lf::uscalfe::MeshFunctionFE<double, double> mf_numeric(fe_space, solution);
    const lf::uscalfe::MeshFunctionGradFE<double, double> mf_numeric_grad(fe_space, solution);

    // Compute the H1 and L2 errors
    std::cout << "\t\t> Computing Error Norms" << std::endl;
    const auto qr_segment = lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), 2*degree-1);
    const auto qr_tria = lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2*degree-1);
    const auto quadrule_provider = [&](const lf::mesh::Entity &entity) {
	const lf::base::RefEl refel = entity.RefEl();
	if (refel == lf::base::RefEl::kTria()) {
	    return qr_tria;
	}
	else if (refel == lf::base::RefEl::kSegment()) {
	    return qr_segment;
	}
	else {
	    return lf::quad::make_QuadRule(refel, 2*degree-1);
	}
    };
    const double H1_err = std::sqrt(lf::uscalfe::IntegrateMeshFunction(*mesh, lf::mesh::utils::squaredNorm(mf_u_grad - mf_numeric_grad), quadrule_provider));
    const double L2_err = std::sqrt(lf::uscalfe::IntegrateMeshFunction(*mesh, lf::mesh::utils::squaredNorm(mf_u - mf_numeric), quadrule_provider));

    // Return the errors
    return {H1_err, L2_err};
}



std::tuple<double, double> computeErrorsLDomain(const std::shared_ptr<const lf::mesh::Mesh> &mesh_p, const std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> & fe_space_p) {
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
}



int main(int argc, char *argv[]) {
    po::options_description desc("Allowed options");
    desc.add_options()
    ("output,o", po::value<std::string>(), "Name of the output file")
    ("max_p,p", po::value<unsigned>(), "Maximum polynomial degree");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("output") == 0 || vm.count("max_p") == 0) {
	std::cout << desc << std::endl;
	exit(1);
    }
    const std::string output_file = vm["output"].as<std::string>();
    const unsigned max_p = vm["max_p"].as<unsigned>();

    const boost::filesystem::path here = __FILE__;
    // Load the unit square mesh
    const auto square_mesh = getSquareDomain();
    // Load the L-shaped domain mesh
    const boost::filesystem::path L_mesh_file = here.parent_path() / "meshes" / "L0.msh";
    auto mesh_factory_L = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    const lf::io::GmshReader reader_L(std::move(mesh_factory_L), L_mesh_file.string());
    const auto L_mesh = reader_L.mesh();

    // Compute the errors for different polynomial degrees
    Eigen::MatrixXd results(max_p, 7);
    for (unsigned p = 1 ; p <= max_p ; ++p) {
	std::cout << "> Polynomial Degree: " << p << std::endl;

	// Solve the problem on the unit square domain
	std::cout << "\t> Unit Square Domain";
	const auto fe_space_square = std::make_shared<FeSpaceLagrangeON<double>>(square_mesh, p);
	std::cout << " (" << fe_space_square->LocGlobMap().NumDofs() << " DOFs)" << std::endl;
	const auto [H1_square, L2_square] = computeErrorsSquareDomain(square_mesh, fe_space_square);
	
	// Solve the problem on the L-shaped domain
	std::cout << "\t> L-shaped Domain";
	const auto fe_space_L = std::make_shared<FeSpaceLagrangeON<double>>(L_mesh, p);
	std::cout << " (" << fe_space_L->LocGlobMap().NumDofs() << " DOFs)" << std::endl;
	const auto [H1_L, L2_L] = computeErrorsLDomain(L_mesh, fe_space_L);

	// Store the computed quantities in the results matrix
	results(p-1, 0) = p;
	results(p-1, 1) = fe_space_square->LocGlobMap().NumDofs();
	results(p-1, 2) = fe_space_L->LocGlobMap().NumDofs();
	results(p-1, 3) = H1_square;
	results(p-1, 4) = L2_square;
	results(p-1, 5) = H1_L;
	results(p-1, 6) = L2_L;
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
