/**
 * @file convergence.cc
 * @brief Produces h-convergence results for the lid driven cavity experiment
 */

#define _USE_MATH_DEFINES
#include <build_system_matrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/io/gmsh_reader.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/tp_triag_mesh_builder.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/refinement/mesh_function_transfer.h>
#include <lf/refinement/refinement.h>
#include <mesh_function_velocity.h>
#include <norms.h>
#include <piecewise_const_element_matrix_provider.h>
#include <piecewise_const_element_vector_provider.h>
#include <solution_to_mesh_data_set.h>

#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

using lf::uscalfe::operator-;
using lf::uscalfe::operator*;

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
 * @brief Outputs L2 and DG norm errors for the lid driven cavity experiment
 *
 * The lid driven cavity problem is solved on a mesh hierarchy of
 * `refinement_level+1` tensor product meshes. The convergence data is printed
 * to stdout in the following order: level vertex_count L2 L2_weighted  DG
 * L2_modified L2_weighted_modified DG_modified
 */
int main() {
  const unsigned refinement_level = 7;

  // Build the mesh
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::utils::TPTriagMeshBuilder builder(std::move(factory));
  builder.setBottomLeftCorner(0, 0);
  builder.setTopRightCorner(1, 1);
  builder.setNumXCells(2);
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
    // No volume forces are present in this experiment
    auto f = [](const Eigen::Vector2d & /*unused*/) -> Eigen::Vector2d {
      return Eigen::Vector2d::Zero();
    };
    // The top lid is driven with velocity 1
    auto dirichlet_funct = [](const lf::mesh::Entity &edge) -> Eigen::Vector2d {
      static constexpr double eps = 1e-10;
      const auto *const geom = edge.Geometry();
      const auto vertices = geom->Global(edge.RefEl().NodeCoords());
      Eigen::Vector2d v;
      v << 1, 0;
      if (vertices(1, 0) <= 1 + eps && vertices(1, 0) >= 1 - eps &&
          vertices(1, 1) <= 1 + eps && vertices(1, 1) >= 1 - eps) {
        return -v;
      }
      return Eigen::Vector2d::Zero();
    };
    const auto [A, rhs] =
        projects::ipdg_stokes::assemble::buildSystemMatrixNoFlow(
            sol.mesh, *(sol.dofh), f, dirichlet_funct, 100,
            lf::quad::make_TriaQR_MidpointRule(), false);
    sol.A = A.makeSparse();
    sol.rhs = rhs;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(sol.A);
    sol.solution = solver.solve(rhs);
    const auto [A_modified, rhs_modified] =
        projects::ipdg_stokes::assemble::buildSystemMatrixNoFlow(
            sol.mesh, *(sol.dofh), f, dirichlet_funct, 100,
            lf::quad::make_TriaQR_MidpointRule(), true);
    sol.A_modified = A_modified.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_modified(
        sol.A_modified);
    sol.solution_modified = solver_modified.solve(rhs_modified);
  }

  // Perform post processing on the data
  const auto fe_space_fine =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(
          solutions.back().mesh);
  projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double, double>
      velocity_exact(fe_space_fine, solutions.back().solution);
  projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double, double>
      velocity_exact_modified(fe_space_fine,
                              solutions.back().solution_modified);
  lf::io::VtkWriter writer(solutions.back().mesh, "result.vtk");
  for (lf::base::size_type lvl = 0; lvl < refinement_level; ++lvl) {
    const auto fe_space_lvl =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(
            solutions[lvl].mesh);
    projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double, double>
        velocity_lvl(fe_space_lvl, solutions[lvl].solution);
    projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double, double>
        velocity_lvl_modified(fe_space_lvl, solutions[lvl].solution_modified);
    lf::refinement::MeshFunctionTransfer velocity_fine(
        *mesh_hierarchy, velocity_lvl, lvl, mesh_hierarchy->NumLevels() - 1);
    lf::refinement::MeshFunctionTransfer velocity_fine_modified(
        *mesh_hierarchy, velocity_lvl_modified, lvl,
        mesh_hierarchy->NumLevels() - 1);
    writer.WriteCellData(concat("v_", solutions[lvl].mesh->NumEntities(2)),
                         velocity_fine);
    writer.WriteCellData(
        concat("v_modified", solutions[lvl].mesh->NumEntities(2)),
        velocity_fine_modified);
    // We need to implement our own binary mesh function for  multiplication
    const auto qr_provider = [](const lf::mesh::Entity &e) {
      return lf::quad::make_QuadRule(e.RefEl(), 0);
    };
    const auto weight =
        lf::mesh::utils::MeshFunctionGlobal([](const Eigen::Vector2d &x) {
          return x[1] >= 0.9 ? (1 - std::cos(M_PI / 0.1 * (1 - x[1]))) / 2 : 1.;
        });
    const auto diff = velocity_fine - velocity_exact;
    const auto diff_weighted = weight * diff;
    const auto diff_modified = velocity_fine_modified - velocity_exact_modified;
    const auto diff_weighted_modified = weight * diff_modified;
    const double L2 = projects::ipdg_stokes::post_processing::L2norm(
        solutions.back().mesh, diff, qr_provider);
    const double L2w = projects::ipdg_stokes::post_processing::L2norm(
        solutions.back().mesh, diff_weighted, qr_provider);
    const double DG = projects::ipdg_stokes::post_processing::DGnorm(
        solutions.back().mesh, diff,
        lf::mesh::utils::MeshFunctionConstant<Eigen::Matrix2d>(
            Eigen::Matrix2d::Zero()),
        qr_provider);
    const double L2_modified = projects::ipdg_stokes::post_processing::L2norm(
        solutions.back().mesh, diff_modified, qr_provider);
    const double L2w_modified = projects::ipdg_stokes::post_processing::L2norm(
        solutions.back().mesh, diff_weighted_modified, qr_provider);
    const double DG_modified = projects::ipdg_stokes::post_processing::DGnorm(
        solutions.back().mesh, diff_modified,
        lf::mesh::utils::MeshFunctionConstant<Eigen::Matrix2d>(
            Eigen::Matrix2d::Zero()),
        qr_provider);
    std::cout << lvl << ' ' << solutions[lvl].mesh->NumEntities(2) << ' ' << L2
              << ' ' << L2w << ' ' << DG << ' ' << L2_modified << ' '
              << L2w_modified << ' ' << DG_modified << std::endl;
  }

  return 0;
}
