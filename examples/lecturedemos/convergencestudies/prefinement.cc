/**
 * @file prefinement.cc
 * @brief Creates convergence plots for experiment 3.2.3.11
 * @author Tobias Rohner
 * @date April 2020
 * @copyright MIT License
 */

#define _USE_MATH_DEFINES

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "fespacelagrangeon.h"

namespace po = boost::program_options;

/**
 * @brief Builds a mesh on [0, 1]^2 from two triangles
 * @returns A shared pointer to a mesh
 */
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
  coords << 0, 1, 0, 0, 0, 1;
  nodes[0] = vertices[0];
  nodes[1] = vertices[1];
  nodes[2] = vertices[3];
  auto geom_tria1 = std::make_unique<lf::geometry::TriaO1>(coords);
  factory.AddEntity(lf::base::RefEl::kTria(), nodes, std::move(geom_tria1));
  coords << 1, 1, 0, 0, 1, 1;
  nodes[0] = vertices[1];
  nodes[1] = vertices[2];
  nodes[2] = vertices[3];
  auto geom_tria2 = std::make_unique<lf::geometry::TriaO1>(coords);
  factory.AddEntity(lf::base::RefEl::kTria(), nodes, std::move(geom_tria2));
  // Build the mesh
  return factory.Build();
}

/**
 * @brief Builds a mesh on [-1, 1]^2 \ (]0, 1[x]-1, 0[) from four triangles
 * @returns A shared pointer to a mesh
 */
std::shared_ptr<lf::mesh::Mesh> getLDomain() {
  lf::mesh::hybrid2d::MeshFactory factory(2);
  // Add the vertices
  std::vector<lf::mesh::MeshFactory::size_type> vertices;
  Eigen::Vector2d vertex_coord;
  vertex_coord << -1, -1;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << 0, -1;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << 0, 0;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << 1, 0;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << 1, 1;
  vertices.push_back(factory.AddPoint(vertex_coord));
  vertex_coord << -1, 1;
  vertices.push_back(factory.AddPoint(vertex_coord));
  // Add the triangles
  Eigen::Matrix<double, Eigen::Dynamic, 3> coords(2, 3);
  lf::mesh::MeshFactory::size_type nodes[3];
  coords << -1, 0, -1, -1, 0, 1;
  nodes[0] = vertices[0];
  nodes[1] = vertices[2];
  nodes[2] = vertices[5];
  auto geom_tria1 = std::make_unique<lf::geometry::TriaO1>(coords);
  factory.AddEntity(lf::base::RefEl::kTria(), nodes, std::move(geom_tria1));
  coords << 0, 1, -1, 0, 1, 1;
  nodes[0] = vertices[2];
  nodes[1] = vertices[4];
  nodes[2] = vertices[5];
  auto geom_tria2 = std::make_unique<lf::geometry::TriaO1>(coords);
  factory.AddEntity(lf::base::RefEl::kTria(), nodes, std::move(geom_tria2));
  coords << -1, 0, 0, -1, -1, 0;
  nodes[0] = vertices[0];
  nodes[1] = vertices[1];
  nodes[2] = vertices[2];
  auto geom_tria3 = std::make_unique<lf::geometry::TriaO1>(coords);
  factory.AddEntity(lf::base::RefEl::kTria(), nodes, std::move(geom_tria3));
  coords << 0, 1, 1, 0, 0, 1;
  nodes[0] = vertices[2];
  nodes[1] = vertices[3];
  nodes[2] = vertices[4];
  auto geom_tria4 = std::make_unique<lf::geometry::TriaO1>(coords);
  factory.AddEntity(lf::base::RefEl::kTria(), nodes, std::move(geom_tria4));
  // Build the mesh
  return factory.Build();
}

/**
 * @brief Get the error of a test problem for the given mesh and FE space
 * @param mesh The mesh on which to solve the PDE
 * @param fe_space The Finite Element Space to use for the computation
 *
 * The test problem is formulated on [0, 1]^2 with Dirichlet boundary conditions
 and has the analytic solution
 * \f[
        u(x) = \sin(\pi x_1)\sin(\pi x_2)
   \f]
 */
std::tuple<double, double> computeErrorsSquareDomain(
    const std::shared_ptr<lf::mesh::Mesh> &mesh,
    const std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>>
        &fe_space) {
  // Get the degree of the fe space
  const auto degree =
      fe_space->ShapeFunctionLayout(lf::base::RefEl::kTria())->Degree();
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
        lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2 * degree - 1)},
       {lf::base::RefEl::kQuad(),
        lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 2 * degree - 1)}});
  lf::assemble::AssembleVectorLocally(0, dofh, element_vector_provider, rhs);

  // Enforce zero dirichlet boundary conditions
  std::cout << "\t\t> Enforcing Boundary Conditions" << std::endl;
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);
  const auto selector = [&](unsigned int idx) -> std::pair<bool, double> {
    const auto &entity = dofh.Entity(idx);
    return {entity.Codim() > 0 && boundary(entity), 0};
  };
  lf::assemble::FixFlaggedSolutionComponents(selector, A_COO, rhs);

  // Solve the LSE using the cholesky decomposition
  std::cout << "\t\t> Solving LSE" << std::endl;
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
  const Eigen::VectorXd solution = solver.solve(rhs);
  const lf::uscalfe::MeshFunctionFE<double, double> mf_numeric(fe_space,
                                                               solution);
  const lf::uscalfe::MeshFunctionGradFE<double, double> mf_numeric_grad(
      fe_space, solution);

  // Compute the H1 and L2 errors
  std::cout << "\t\t> Computing Error Norms" << std::endl;
  const auto qr_segment =
      lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), 2 * degree - 1);
  const auto qr_tria =
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2 * degree - 1);
  const auto qr_quad =
      lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 2 * degree - 1);
  const auto quadrule_provider = [&](const lf::mesh::Entity &entity) {
    const lf::base::RefEl refel = entity.RefEl();
    if (refel == lf::base::RefEl::kTria()) {
      return qr_tria;
    } else if (refel == lf::base::RefEl::kSegment()) {
      return qr_segment;
    } else if (refel == lf::base::RefEl::kQuad()) {
      return qr_quad;
    } else {
      return lf::quad::make_QuadRule(refel, 2 * degree - 1);
    }
  };
  const double H1_err = std::sqrt(lf::uscalfe::IntegrateMeshFunction(
      *mesh, lf::mesh::utils::squaredNorm(mf_u_grad - mf_numeric_grad),
      quadrule_provider));
  const double L2_err = std::sqrt(lf::uscalfe::IntegrateMeshFunction(
      *mesh, lf::mesh::utils::squaredNorm(mf_u - mf_numeric),
      quadrule_provider));

  // Return the errors
  return {H1_err, L2_err};
}

/**
 * @brief Get the error of a test problem for the given mesh and FE space
 * @param mesh The mesh on which to solve the PDE
 * @param fe_space The Finite Element Space to use for the computation
 *
 * The test problem is formulated on [-1, 1]^2 \ (]0, 1[x]-1, 0[) with Dirichlet
 boundary conditions and has the analytic solution
 * \f[
        u(r, \phi) = r^{\frac{2}{3}}\sin(\frac{2}{3}\phi)
   \f]
 */
std::tuple<double, double> computeErrorsLDomain(
    const std::shared_ptr<lf::mesh::Mesh> &mesh,
    const std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>>
        &fe_space) {
  // The analytic solution
  const auto u = [](const Eigen::Vector2d &x) -> double {
    const double r = x.norm();
    double phi = std::atan2(x[1], x[0]);
    if (phi < 0) {
      phi += 2 * M_PI;
    }
    return std::pow(r, 2. / 3) * std::sin(2. / 3 * phi);
  };
  lf::mesh::utils::MeshFunctionGlobal mf_u(u);
  // The gradient of the analytic solution
  const auto u_grad = [](const Eigen::Vector2d &x) -> Eigen::Vector2d {
    const double r = x.norm();
    double phi = std::atan2(x[1], x[0]);
    if (phi < 0) {
      phi += 2 * M_PI;
    }
    Eigen::Vector2d grad;
    grad[0] = 2. / 3 * std::pow(r, -4. / 3) *
              (x[0] * std::sin(2. / 3 * phi) - x[1] * std::cos(2. / 3 * phi));
    grad[1] = 2. / 3 * std::pow(r, -4. / 3) *
              (x[1] * std::sin(2. / 3 * phi) + x[0] * std::cos(2. / 3 * phi));
    return grad;
  };
  lf::mesh::utils::MeshFunctionGlobal mf_u_grad(u_grad);

  // Get a few useful variables
  const lf::assemble::DofHandler &dofh = fe_space->LocGlobMap();
  const unsigned degree =
      fe_space->ShapeFunctionLayout(lf::base::RefEl::kTria())->Degree();

  // Assemble the system matrix
  std::cout << "\t\t> Assembling System Matrix" << std::endl;
  const lf::mesh::utils::MeshFunctionConstant<double> mf_alpha(1);
  const lf::mesh::utils::MeshFunctionConstant<double> mf_gamma(0);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider element_matrix_provider(
      fe_space, mf_alpha, mf_gamma);
  lf::assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, element_matrix_provider,
                                      A_COO);

  // The right hand side is zero because we have no load
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());

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
    const auto shape_function_layout =
        fe_space->ShapeFunctionLayout(entity.RefEl());
    const auto num_dofs = dofh.NumLocalDofs(entity);
    const auto glob_dof_idxs = dofh.GlobalDofIndices(entity);
    int dof_idx;
    for (dof_idx = 0; dof_idx < num_dofs; ++dof_idx) {
      if (glob_dof_idxs[dof_idx] == idx) {
        break;
      }
    }
    const Eigen::VectorXd eval_node =
        shape_function_layout->EvaluationNodes().col(dof_idx);
    const Eigen::Vector2d pos = geom->Global(eval_node);
    // Return the value on the boundary at the position of the evaluation node
    // corresponding to the dof
    return {true, u(pos)};
  };
  lf::assemble::FixFlaggedSolutionComponents(selector, A_COO, rhs);

  // Solve the LSE using Cholesky decomposition
  std::cout << "\t\t> Solving LSE" << std::endl;
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(A);
  const Eigen::VectorXd solution = solver.solve(rhs);
  const lf::uscalfe::MeshFunctionFE<double, double> mf_numeric(fe_space,
                                                               solution);
  const lf::uscalfe::MeshFunctionGradFE<double, double> mf_numeric_grad(
      fe_space, solution);

  // Compute the H1 and L2 errors
  std::cout << "\t\t> Computing Error Norms" << std::endl;
  const auto qr_segment =
      lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), 2 * degree - 1);
  const auto qr_tria =
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), 2 * degree - 1);
  const auto qr_quad =
      lf::quad::make_QuadRule(lf::base::RefEl::kQuad(), 2 * degree - 1);
  const auto quadrule_provider = [&](const lf::mesh::Entity &entity) {
    const lf::base::RefEl refel = entity.RefEl();
    if (refel == lf::base::RefEl::kTria()) {
      return qr_tria;
    } else if (refel == lf::base::RefEl::kSegment()) {
      return qr_segment;
    } else if (refel == lf::base::RefEl::kQuad()) {
      return qr_quad;
    } else {
      return lf::quad::make_QuadRule(refel, 2 * degree - 1);
    }
  };
  const double H1_err = std::sqrt(lf::uscalfe::IntegrateMeshFunction(
      *mesh, lf::mesh::utils::squaredNorm(mf_u_grad - mf_numeric_grad),
      quadrule_provider));
  const double L2_err = std::sqrt(lf::uscalfe::IntegrateMeshFunction(
      *mesh, lf::mesh::utils::squaredNorm(mf_u - mf_numeric),
      quadrule_provider));

  // Return the errors
  return {H1_err, L2_err};
}

int main(int argc, char *argv[]) {
  po::options_description desc("Allowed options");
  desc.add_options()("output,o", po::value<std::string>(),
                     "Name of the output file")(
      "max_p,p", po::value<unsigned>(), "Maximum polynomial degree");
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
  const auto L_mesh = getLDomain();

  // Compute the errors for different polynomial degrees
  Eigen::MatrixXd results(max_p, 7);
  for (unsigned p = 1; p <= max_p; ++p) {
    std::cout << "> Polynomial Degree: " << p << std::endl;

    // Solve the problem on the unit square domain
    std::cout << "\t> Unit Square Domain";
    const auto fe_space_square =
        std::make_shared<FeSpaceLagrangeON<double>>(square_mesh, p);
    std::cout << " (" << fe_space_square->LocGlobMap().NumDofs() << " DOFs)"
              << std::endl;
    const auto [H1_square, L2_square] =
        computeErrorsSquareDomain(square_mesh, fe_space_square);

    // Solve the problem on the L-shaped domain
    std::cout << "\t> L-shaped Domain";
    const auto fe_space_L =
        std::make_shared<FeSpaceLagrangeON<double>>(L_mesh, p);
    std::cout << " (" << fe_space_L->LocGlobMap().NumDofs() << " DOFs)"
              << std::endl;
    const auto [H1_L, L2_L] = computeErrorsLDomain(L_mesh, fe_space_L);

    // Store the computed quantities in the results matrix
    results(p - 1, 0) = p;
    results(p - 1, 1) = fe_space_square->LocGlobMap().NumDofs();
    results(p - 1, 2) = fe_space_L->LocGlobMap().NumDofs();
    results(p - 1, 3) = H1_square;
    results(p - 1, 4) = L2_square;
    results(p - 1, 5) = H1_L;
    results(p - 1, 6) = L2_L;
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
