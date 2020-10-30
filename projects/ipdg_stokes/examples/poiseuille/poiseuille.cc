/**
 * @file poiseuille.cc
 * @brief Solve for the poiseuille velocity profile
 */

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <lf/assemble/assemble.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/base/base.h>
#include <lf/io/io.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/hybrid2d/mesh.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/tp_triag_mesh_builder.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>

#include <build_system_matrix.h>
#include <mesh_function_velocity.h>
#include <norms.h>
#include <piecewise_const_element_matrix_provider.h>
#include <piecewise_const_element_vector_provider.h>
#include <solution_to_mesh_data_set.h>

using lf::uscalfe::operator-;

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
 * @brief Compute the analytic flow velocity of the poiseuille flow
 * @param h Half the distance between the two plates
 * @param flowrate The flowrate
 * @param y The y position at which to get the flow velocity. It must be in the
 * interval [-h, h]
 * @returns The flow velocity at point y
 */
double poiseuilleVelocity(double h, double flowrate, double y) {
  const double u_max = flowrate * 3 / 4 / h;
  return u_max * (1 - y * y / h / h);
}

/**
 * @brief stores information to recover convergence properties
 */
struct ProblemSolution {
  std::shared_ptr<const lf::mesh::Mesh> mesh;
  std::shared_ptr<const lf::assemble::DofHandler> dofh;
  Eigen::SparseMatrix<double> A;
  Eigen::SparseMatrix<double> A_modified;
  Eigen::VectorXd rhs;
  Eigen::VectorXd solution;
  Eigen::VectorXd solution_modified;
};

/**
 * @brief Prints the L2 and DG norm errors of the poiseuille flow experiment
 */
int main(int argc, char *argv[]) {
  const unsigned refinement_level = 6;
  const double flowrate = 1;
  const double h = 0.25;

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << "[regular,stretched]" << std::endl;
    std::cerr << "\tregular  : Builds a mesh where the resolution in the x- "
                 "and y-direction are equal"
              << std::endl;
    std::cerr << "\tirregular: Builds a mesh where the resolution in the "
                 "x-direction is 10 times less than in the y-direction"
              << std::endl;
    exit(1);
  }

  // Read the mesh from the gmsh file
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::utils::TPTriagMeshBuilder builder(std::move(factory));
  builder.setBottomLeftCorner(0, 0);
  builder.setTopRightCorner(2, 2 * h);
  if (strcmp(argv[1], "irregular") == 0) {
    builder.setNumXCells(std::max(2, 2 * static_cast<int>(1.0 / h) / 10));
  } else {
    builder.setNumXCells(2 * static_cast<int>(1.0 / h));
  }
  builder.setNumYCells(2);
  auto mesh0 = builder.Build();

  // Generate a mesh hierarchy by regular refinement
  const auto mesh_hierarchy =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh0,
                                                              refinement_level);

  // Solve the problem for each mesh in the hierarchy
  std::vector<ProblemSolution> solutions(refinement_level + 1);
  for (lf::base::size_type lvl = 0; lvl < refinement_level + 1; ++lvl) {
    const auto &mesh = mesh_hierarchy->getMesh(lvl);
    auto &sol = solutions[lvl];
    sol.mesh = mesh;
    sol.dofh = std::shared_ptr<const lf::assemble::DofHandler>(
        new lf::assemble::UniformFEDofHandler(
            mesh, {{lf::base::RefEl::kPoint(), 1},
                   {lf::base::RefEl::kSegment(), 1}}));
    // No volume forces are present
    auto f = [](const Eigen::Vector2d & /*unused*/) -> Eigen::Vector2d {
      return Eigen::Vector2d::Zero();
    };
    // Enforce a Poiseuille in- and outflow and no-slip boundary conditions at
    // the tube boundaries
    auto dirichlet_funct =
        [&](const lf::mesh::Entity &edge) -> Eigen::Vector2d {
      static constexpr double eps = 1e-10;
      const auto geom = edge.Geometry();
      const auto vertices = geom->Global(edge.RefEl().NodeCoords());
      const Eigen::Vector2d midpoint = vertices.rowwise().sum() / 2;
      Eigen::Vector2d v = Eigen::Vector2d::Zero();
      if (vertices(0, 0) >= -eps && vertices(0, 0) <= eps &&
          vertices(0, 1) >= -eps && vertices(0, 1) <= eps) {
        // The edge is part of the inflow boundary
        v << poiseuilleVelocity(h, flowrate, midpoint[1] - h), 0;
      }
      if (vertices(0, 0) >= 2 - eps && vertices(0, 0) <= 2 + eps &&
          vertices(0, 1) >= 2 - eps && vertices(0, 1) <= 2 + eps) {
        // The edge is part of the outflow boundary
        v << poiseuilleVelocity(h, flowrate, midpoint[1] - h), 0;
      }
      return v;
    };
    const auto [A, rhs, offset_function] =
        projects::ipdg_stokes::assemble::buildSystemMatrixInOutFlow(
            sol.mesh, *(sol.dofh), f, dirichlet_funct, 100,
            lf::quad::make_TriaQR_MidpointRule(), false);
    sol.A = A.makeSparse();
    sol.rhs = rhs;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(sol.A);
    sol.solution = solver.solve(rhs) + offset_function;
    const auto [A_modified, rhs_modified, offset_function_modified] =
        projects::ipdg_stokes::assemble::buildSystemMatrixInOutFlow(
            sol.mesh, *(sol.dofh), f, dirichlet_funct, 100,
            lf::quad::make_TriaQR_MidpointRule(), true);
    sol.A_modified = A_modified.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_modified(
        sol.A_modified);
    sol.solution_modified =
        solver_modified.solve(rhs_modified) + offset_function_modified;
  }

  // Compute the analytic solution of the problem
  auto analytic_velocity = [&](const Eigen::Vector2d &x) -> Eigen::Vector2d {
    const double u_max = flowrate * 3 / 4 / h;
    Eigen::Vector2d v;
    v << u_max * (1 - (x[1] - h) * (x[1] - h) / h / h), 0;
    return v;
  };
  auto analytic_gradient = [&](const Eigen::Vector2d &x) -> Eigen::Matrix2d {
    const double u_max = flowrate * 3 / 4 / h;
    Eigen::Matrix2d g;
    g << 0, 0, u_max * 2 * (x[1] - h) / h / h, 0;
    return g;
  };

  // Perform post processing on the data
  lf::io::VtkWriter writer(solutions.back().mesh, "result.vtk");
  for (lf::base::size_type lvl = 0; lvl <= refinement_level; ++lvl) {
    const auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(
            solutions[lvl].mesh);
    const auto velocity_exact =
        lf::mesh::utils::MeshFunctionGlobal(analytic_velocity);
    const auto gradient_exact =
        lf::mesh::utils::MeshFunctionGlobal(analytic_gradient);
    const auto velocity =
        projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double,
                                                                     double>(
            fe_space, solutions[lvl].solution);
    const auto velocity_modified =
        projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double,
                                                                     double>(
            fe_space, solutions[lvl].solution_modified);
    // Store the result on the finest mesh to vtk
    if (lvl == refinement_level) {
      const auto v = *lf::mesh::utils::make_LambdaMeshDataSet(
          [&](const lf::mesh::Entity &entity) -> Eigen::Vector2d {
            const Eigen::Vector2d center =
                entity.Geometry()
                    ->Global(entity.RefEl().NodeCoords())
                    .rowwise()
                    .sum() /
                entity.RefEl().NumNodes();
            return analytic_velocity(center);
          });
      writer.WriteCellData(concat("v_", solutions[lvl].mesh->NumEntities(2)),
                           velocity);
      writer.WriteCellData(
          concat("v__modified_", solutions[lvl].mesh->NumEntities(2)),
          *lf::mesh::utils::make_LambdaMeshDataSet(
              [&](const lf::mesh::Entity &e) {
                return velocity_modified(e,
                                         Eigen::Vector2d::Constant(1. / 3))[0];
              }));
      writer.WriteCellData("analytic", v);
    }
    // Compute the error in the velocity
    auto diff_v = velocity - velocity_exact;
    auto diff_v_modified = velocity_modified - velocity_exact;
    // Compute the error in the gradient of the velocity
    auto diff_g = -gradient_exact;
    auto diff_g_modified = -gradient_exact;
    const auto qr_provider = [](const lf::mesh::Entity &e) {
      return lf::quad::make_QuadRule(e.RefEl(), 0);
    };
    const double L2 = projects::ipdg_stokes::post_processing::L2norm(
        solutions[lvl].mesh, diff_v, qr_provider);
    const double DG = projects::ipdg_stokes::post_processing::DGnorm(
        solutions[lvl].mesh, diff_v, diff_g, qr_provider);
    const double L2_modified = projects::ipdg_stokes::post_processing::L2norm(
        solutions[lvl].mesh, diff_v_modified, qr_provider);
    const double DG_modified = projects::ipdg_stokes::post_processing::DGnorm(
        solutions[lvl].mesh, diff_v_modified, diff_g_modified, qr_provider);
    std::cout << lvl << ' ' << solutions[lvl].mesh->NumEntities(2) << ' ' << L2
              << ' ' << DG << ' ' << L2_modified << ' ' << DG_modified
              << std::endl;
  }

  return 0;
}
