/**
 * @file nested_cylinders.cc
 * @brief Compute the convergence of the nested cylinders experiment
 */

#include <algorithm>
#include <boost/filesystem.hpp>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>

#include <lf/assemble/dofhandler.h>
#include <lf/io/gmsh_reader.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/lambda_mesh_data_set.h>
#include <lf/quad/quad.h>

#include <boost/program_options.hpp>

#include <annulus_triag_mesh_builder.h>
#include <build_system_matrix.h>
#include <mesh_function_interpolation.h>
#include <mesh_function_velocity.h>
#include <norms.h>
#include <piecewise_const_element_matrix_provider.h>
#include <piecewise_const_element_vector_provider.h>
#include <solution_to_mesh_data_set.h>

using lf::uscalfe::operator-;

/**
 * @brief Solve the nested cylinders problem with zero potential at the boundary
 * @param mesh A shared pointer to the mesh on which to solve the PDE
 * @param dofh The dofhandler used for the simulation
 * @param r The radius of the inner cylinder
 * @param R The radius of the outer cylinder
 * @param omega1 The angular velocity of the inner cylinder
 * @param omega2 The angular velocity of the outer cylinder
 * @param modified_penalty If true, use the modified penalty term instead of the
 * original one
 * @returns A vector of basis function coefficients for the solution
 */
Eigen::VectorXd solveNestedCylindersZeroBC(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh, double r, double R, double omega1,
    double omega2, bool modified_penalty) {
  // The volume forces are equal to zero everywhere
  auto f = [](const Eigen::Vector2d & /*unused*/) {
    return Eigen::Vector2d::Zero();
  };
  // Drive the inner and outer cylinder with omega1 and omega2
  const double eps = 1e-10;
  auto dirichlet_funct = [&](const lf::mesh::Entity &edge) -> Eigen::Vector2d {
    const auto geom = edge.Geometry();
    const auto vertices = geom->Global(edge.RefEl().NodeCoords());
    if (vertices.col(0).norm() <= R + eps &&
        vertices.col(0).norm() >= R - eps &&
        vertices.col(1).norm() <= R + eps &&
        vertices.col(1).norm() >= R - eps) {
      return omega2 * R * (vertices.col(1) - vertices.col(0)).normalized();
    }
    if (vertices.col(0).norm() <= r + eps &&
        vertices.col(0).norm() >= r - eps &&
        vertices.col(1).norm() <= r + eps &&
        vertices.col(1).norm() >= r - eps) {
      return omega1 * r * (vertices.col(0) - vertices.col(1)).normalized();
    }
    return Eigen::Vector2d::Zero();
  };
  const auto dirichlet =
      *lf::mesh::utils::make_LambdaMeshDataSet(dirichlet_funct);

  // Asemble the LSE
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);
  // Assemble the Matrix
  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());
  const projects::ipdg_stokes::assemble::PiecewiseConstElementMatrixProvider
      elem_mat_provider(100, boundary, modified_penalty);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elem_mat_provider, A);
  // Assemble the right hand side
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());
  const projects::ipdg_stokes::assemble::PiecewiseConstElementVectorProvider
      elem_vec_provider(100, f, lf::quad::make_TriaQR_MidpointRule(), boundary,
                        dirichlet);
  lf::assemble::AssembleVectorLocally(0, dofh, elem_vec_provider, rhs);

  // Enforce the no-flow boundary conditions
  auto selector = [&](lf::base::size_type idx) -> std::pair<bool, double> {
    const auto &entity = dofh.Entity(idx);
    if (entity.RefEl() == lf::base::RefElType::kPoint && boundary(entity)) {
      return {true, 0};
    }
    return {false, 0};
  };
  lf::assemble::FixFlaggedSolutionComponents(selector, A, rhs);

  // Solve the LSE using sparse LU
  Eigen::SparseMatrix<double> As = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(As);
  return solver.solve(rhs);
}

/**
 * @brief Solve the nested cylinders problem with constant potential at the
 * boundary
 * @param mesh A shared pointer to the mesh on which to solve the PDE
 * @param dofh The dofhandler used for the simulation
 * @param r The radius of the inner cylinder
 * @param R The radius of the outer cylinder
 * @param omega1 The angular velocity of the inner cylinder
 * @param omega2 The angular velocity of the outer cylinder
 * @param modified_penalty If true, use the modified penalty term instead of the
 * original one
 * @returns A vector of basis function coefficients for the solution
 */
Eigen::VectorXd solveNestedCylindersNonzeroBC(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh, double r, double R, double omega1,
    double omega2, bool modified_penalty) {
  // The volume forces are equal to zero everywhere
  auto f = [](const Eigen::Vector2d & /*unused*/) {
    return Eigen::Vector2d::Zero();
  };
  // Drive the inner and outer cylinder with omega1 and omega2
  const double eps = 1e-10;
  auto dirichlet_funct = [&](const lf::mesh::Entity &edge) -> Eigen::Vector2d {
    const auto geom = edge.Geometry();
    const auto vertices = geom->Global(edge.RefEl().NodeCoords());
    if (vertices.col(0).norm() <= R + eps &&
        vertices.col(0).norm() >= R - eps &&
        vertices.col(1).norm() <= R + eps &&
        vertices.col(1).norm() >= R - eps) {
      return omega2 * R * (vertices.col(1) - vertices.col(0)).normalized();
    }
    if (vertices.col(0).norm() <= r + eps &&
        vertices.col(0).norm() >= r - eps &&
        vertices.col(1).norm() <= r + eps &&
        vertices.col(1).norm() >= r - eps) {
      return omega1 * r * (vertices.col(0) - vertices.col(1)).normalized();
    }
    return Eigen::Vector2d::Zero();
  };
  const auto dirichlet =
      *lf::mesh::utils::make_LambdaMeshDataSet(dirichlet_funct);

  // Asemble the LSE
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);
  // Assemble the Matrix
  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());
  const projects::ipdg_stokes::assemble::PiecewiseConstElementMatrixProvider
      elem_mat_provider(100, boundary, modified_penalty);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elem_mat_provider, A);
  // Assemble the right hand side
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofh.NumDofs());
  const projects::ipdg_stokes::assemble::PiecewiseConstElementVectorProvider
      elem_vec_provider(100, f, lf::quad::make_TriaQR_MidpointRule(), boundary,
                        dirichlet);
  lf::assemble::AssembleVectorLocally(0, dofh, elem_vec_provider, rhs);

  // Combine the basis functions on the inside boundary to a single one
  const lf::base::size_type M_outer = 0;
  const lf::base::size_type M_inner = 1;
  // Create a mapping of DOF indices to remove the afterwards unused DOFs
  std::vector<lf::base::size_type> dofmap(dofh.NumDofs());
  lf::base::size_type idx = 2;
  for (lf::base::size_type dof = 0; dof < dofh.NumDofs(); ++dof) {
    const auto &entity = dofh.Entity(dof);
    const auto geom = entity.Geometry();
    if (entity.RefEl() == lf::base::RefElType::kPoint && boundary(entity)) {
      if (geom->Global(entity.RefEl().NodeCoords()).norm() > R - eps) {
        dofmap[dof] = M_outer;
      } else {
        dofmap[dof] = M_inner;
      }
    } else {
      dofmap[dof] = idx++;
    }
  }
  // Apply this mapping to the triplets of the matrix
  std::for_each(A.triplets().begin(), A.triplets().end(),
                [&](Eigen::Triplet<double> &trip) {
                  trip = Eigen::Triplet<double>(
                      dofmap[trip.row()], dofmap[trip.col()], trip.value());
                });
  // Apply the mapping to the right hand side vector
  Eigen::VectorXd rhs_mapped = Eigen::VectorXd::Zero(idx);
  for (lf::base::size_type dof = 0; dof < dofh.NumDofs(); ++dof) {
    rhs_mapped[dofmap[dof]] += rhs[dof];
  }

  // Set the potential on the outer boundary to zero
  auto selector = [&](lf::base::size_type idx) -> std::pair<bool, double> {
    return {idx == M_outer, 0};
  };
  lf::assemble::FixFlaggedSolutionComponents(selector, A, rhs);

  // Solve the LSE using sparse LU
  Eigen::SparseMatrix<double> As_mapped = A.makeSparse().block(0, 0, idx, idx);
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(As_mapped);
  const Eigen::VectorXd sol_mapped = solver.solve(rhs_mapped);

  // Apply the inverse mapping to recovr the basis function coefficients for the
  // original basis functions
  Eigen::VectorXd sol = Eigen::VectorXd::Zero(dofh.NumDofs());
  for (lf::base::size_type dof = 0; dof < dofh.NumDofs(); ++dof) {
    sol[dof] = sol_mapped[dofmap[dof]];
  }

  return sol;
}

/**
 * @brief Concatenate objects defining an operator<<(std::ostream&)
 * @param args A variadic pack of objects implementing
 * `operator<<(std::ostream&)`
 * @returns A string with the objects concatenated
 */
template <typename... Args>
static std::string concat(Args &&... args) {
  std::ostringstream ss;
  (ss << ... << args);
  return ss.str();
}

/**
 * @brief outputs the L2 and DG norm errors for the nested cylinders experiment
 *
 * Three different command line arguments can be provided:
 *   - builder
 *   - files
 *   - irregular
 *
 * Providing builder as a command line argument will build the meshes with
 * #AnnulusTriagMeshBuilder. Providing files as a command line argument will use
 * uniform meshes generated by GMSH. Providing irregular as a command line
 * argument will use meshes with a sudden jump in mesh resolution.
 */
int main(int argc, char *argv[]) {
  const double r = 0.25;
  const double R = 1;
  const double omega1 = 0;
  const double omega2 = 1;

  // Parse the command line options
  std::string mesh_selection;
  boost::program_options::options_description desc{"Options"};
  desc.add_options()("help,h", "Help Screen")(
      "type", boost::program_options::value<std::string>(&mesh_selection),
      "Type of mesh to use. Either 'builder', 'files' or 'irregular'");
  boost::program_options::positional_options_description pos_desc;
  pos_desc.add("type", 1);
  boost::program_options::command_line_parser parser{argc, argv};
  parser.options(desc).positional(pos_desc).allow_unregistered();
  boost::program_options::parsed_options po = parser.run();
  boost::program_options::variables_map vm;
  boost::program_options::store(po, vm);
  boost::program_options::notify(vm);
  if (vm.count("help") != 0U) {
    std::cout << desc << std::endl;
  }

  std::vector<std::shared_ptr<lf::mesh::Mesh>> meshes;
  if (mesh_selection == "files") {
    // Read the mesh from the gmsh file
    boost::filesystem::path meshpath = __FILE__;
    meshpath = meshpath.parent_path();
    for (int i = 0; i <= 4; ++i) {
      const auto meshfile = meshpath / concat("annulus", std::setw(2),
                                              std::setfill('0'), i, ".msh");
      std::unique_ptr<lf::mesh::MeshFactory> factory =
          std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
      lf::io::GmshReader reader(std::move(factory), meshfile.string());
      meshes.push_back(reader.mesh());
    }
  } else if (mesh_selection == "builder") {
    // Build a sequence of meshes
    for (unsigned i = 0U; i < 8U; ++i) {
      const unsigned nx = 4U << i;
      const double dx = 2 * M_PI * (r + R) / 2 / nx;
      const unsigned ny = std::max(static_cast<unsigned>((R - r) / dx), 1U);

      // Build the mesh
      std::unique_ptr<lf::mesh::MeshFactory> factory =
          std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
      projects::ipdg_stokes::mesh::AnnulusTriagMeshBuilder builder(
          std::move(factory));
      builder.setBottomLeftCorner(0, r);
      builder.setTopRightCorner(1, R);
      builder.setNumXCells(nx);
      builder.setNumYCells(ny);
      meshes.push_back(builder.Build());
    }
  } else if (mesh_selection == "irregular") {
    boost::filesystem::path meshpath = __FILE__;
    const auto mesh_irregular_path =
        meshpath.parent_path() / "annulus_irregular.msh";
    const auto mesh_irregular_inverted_path =
        meshpath.parent_path() / "annulus_irregular_inverted.msh";
    std::unique_ptr<lf::mesh::MeshFactory> factory_irregular =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    std::unique_ptr<lf::mesh::MeshFactory> factory_irregular_inverted =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader_irregular(std::move(factory_irregular),
                                        mesh_irregular_path.string());
    lf::io::GmshReader reader_irregular_inverted(
        std::move(factory_irregular_inverted),
        mesh_irregular_inverted_path.string());
    meshes.push_back(reader_irregular.mesh());
    meshes.push_back(reader_irregular_inverted.mesh());
  } else {
    std::cout << desc << std::endl;
    exit(1);
  }

  // Compute the analytic solution of the problem
  const double C1 = 2 * (omega1 * r * r - omega2 * R * R) / (r * r - R * R);
  const double C2 = ((omega1 - omega2) * r * r * R * R) / (R * R - r * r);
  auto analytic_velocity = [&](const Eigen::Vector2d &x) -> Eigen::Vector2d {
    const double radius = x.norm();
    Eigen::Vector2d vec;
    vec << x[1], -x[0];
    vec.normalize();
    return -(0.5 * C1 * radius + C2 / radius) * vec;
  };
  auto analytic_gradient = [&](const Eigen::Vector2d &x) -> Eigen::Matrix2d {
    const double r2 = x.squaredNorm();
    Eigen::Matrix2d g;
    g << 2 * C2 * x[0] * x[1] / r2 / r2,
        -C1 / 2 - (C2 * r2 - 2 * C2 * x[1] * x[1]) / r2 / r2,
        C1 / 2 + (C2 * r2 - 2 * C2 * x[0] * x[0]) / r2 / r2,
        -2 * C2 * x[0] * x[1] / r2 / r2;
    return g;
  };

  // Solve the problem on each mesh and compute the error
  for (const auto &mesh : meshes) {
    lf::assemble::UniformFEDofHandler dofh(
        mesh,
        {{lf::base::RefEl::kPoint(), 1}, {lf::base::RefEl::kSegment(), 1}});
    const auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
    const Eigen::VectorXd solution_zero =
        solveNestedCylindersZeroBC(mesh, dofh, r, R, omega1, omega2, false);
    const Eigen::VectorXd solution_nonzero =
        solveNestedCylindersNonzeroBC(mesh, dofh, r, R, omega1, omega2, false);
    const Eigen::VectorXd solution_zero_modified =
        solveNestedCylindersZeroBC(mesh, dofh, r, R, omega1, omega2, true);
    const Eigen::VectorXd solution_nonzero_modified =
        solveNestedCylindersNonzeroBC(mesh, dofh, r, R, omega1, omega2, true);
    // Create mesh functions for the analytic and numerical solutions
    const auto velocity_exact =
        lf::uscalfe::MeshFunctionGlobal(analytic_velocity);
    const auto grad_exact = lf::uscalfe::MeshFunctionGlobal(analytic_gradient);
    const auto velocity_zero =
        projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double,
                                                                     double>(
            fe_space, solution_zero);
    const auto velocity_nonzero =
        projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double,
                                                                     double>(
            fe_space, solution_nonzero);
    const auto velocity_zero_modified =
        projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double,
                                                                     double>(
            fe_space, solution_zero_modified);
    const auto velocity_nonzero_modified =
        projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double,
                                                                     double>(
            fe_space, solution_nonzero_modified);
    // Store the solution
    lf::io::VtkWriter writer_zero(
        mesh, concat("result_zero_", dofh.NumDofs(), ".vtk"));
    lf::io::VtkWriter writer_nonzero(
        mesh, concat("result_nonzero_", dofh.NumDofs(), ".vtk"));
    writer_zero.WriteCellData("velocity", velocity_zero);
    writer_zero.WriteCellData("velocity_modified", velocity_zero_modified);
    writer_nonzero.WriteCellData("velocity", velocity_nonzero);
    writer_nonzero.WriteCellData("velocity_modified", velocity_nonzero_modified);
    writer_zero.WriteCellData("analytic", lf::uscalfe::MeshFunctionGlobal(analytic_velocity));
    writer_nonzero.WriteCellData("analytic", lf::uscalfe::MeshFunctionGlobal(analytic_velocity));
    // Compute the difference between the numerical and the analytical solution
    auto diff_velocity_zero = velocity_zero - velocity_exact;
    auto diff_velocity_zero_modified = velocity_zero_modified - velocity_exact;
    auto diff_velocity_nonzero = velocity_nonzero - velocity_exact;
    auto diff_velocity_nonzero_modified = velocity_nonzero_modified - velocity_exact;
    auto diff_gradient_zero = -grad_exact;
    auto diff_gradient_zero_modified = -grad_exact;
    auto diff_gradient_nonzero = -grad_exact;
    auto diff_gradient_nonzero_modified = -grad_exact;
    const auto qr_provider = [](const lf::mesh::Entity &e) {
      return lf::quad::make_QuadRule(e.RefEl(), 0);
    };
    const double L2_zero = projects::ipdg_stokes::post_processing::L2norm(
        mesh, diff_velocity_zero, qr_provider);
    const double L2_nonzero = projects::ipdg_stokes::post_processing::L2norm(
        mesh, diff_velocity_nonzero, qr_provider);
    ;
    const double DG_zero = projects::ipdg_stokes::post_processing::DGnorm(
        mesh, diff_velocity_zero, diff_gradient_zero, qr_provider);
    const double DG_nonzero = projects::ipdg_stokes::post_processing::DGnorm(
        mesh, diff_velocity_nonzero, diff_gradient_nonzero, qr_provider);
    const double L2_zero_modified =
        projects::ipdg_stokes::post_processing::L2norm(
            mesh, diff_velocity_zero_modified, qr_provider);
    const double L2_nonzero_modified =
        projects::ipdg_stokes::post_processing::L2norm(
            mesh, diff_velocity_nonzero_modified, qr_provider);
    const double DG_zero_modified =
        projects::ipdg_stokes::post_processing::DGnorm(
            mesh, diff_velocity_zero_modified, diff_gradient_zero_modified,
            qr_provider);
    const double DG_nonzero_modified =
        projects::ipdg_stokes::post_processing::DGnorm(
            mesh, diff_velocity_nonzero_modified,
            diff_gradient_nonzero_modified, qr_provider);
    std::cout << mesh->NumEntities(2) << ' ' << L2_zero << ' ' << DG_zero << ' '
              << L2_nonzero << ' ' << DG_nonzero << ' ' << L2_zero_modified
              << ' ' << DG_zero_modified << ' ' << L2_nonzero_modified << ' '
              << DG_nonzero_modified << std::endl;
  }

  return 0;
}
