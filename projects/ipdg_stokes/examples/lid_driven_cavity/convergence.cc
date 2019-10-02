/**
 * @file convergence.cc
 * @brief Produces h-convergence results for the lid driven cavity experiment
 */

#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

#include <lf/assemble/dofhandler.h>
#include <lf/io/gmsh_reader.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/tp_triag_mesh_builder.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>

#include <build_system_matrix.h>
#include <mesh_hierarchy_function.h>
#include <norms.h>
#include <piecewise_const_element_matrix_provider.h>
#include <piecewise_const_element_vector_provider.h>
#include <solution_to_mesh_data_set.h>

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
  std::shared_ptr<lf::mesh::MeshFactory> factory =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::mesh::hybrid2d::TPTriagMeshBuilder builder(factory);
  builder.setBottomLeftCorner(0, 0);
  builder.setTopRightCorner(1, 1);
  builder.setNoXCells(2);
  builder.setNoYCells(2);
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
      const auto geom = edge.Geometry();
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

  // Bring the solutions on the meshes to the finest mesh in order to store them
  // to vtk
  std::map<lf::base::size_type,
           std::function<Eigen::Vector2d(const lf::mesh::Entity &,
                                         const Eigen::Vector2d &)>>
      velocity;
  std::map<lf::base::size_type,
           std::function<Eigen::Matrix2d(const lf::mesh::Entity &,
                                         const Eigen::Vector2d &)>>
      gradient;
  std::map<lf::base::size_type,
           std::function<Eigen::Vector2d(const lf::mesh::Entity &,
                                         const Eigen::Vector2d &)>>
      velocity_modified;
  std::map<lf::base::size_type,
           std::function<Eigen::Matrix2d(const lf::mesh::Entity &,
                                         const Eigen::Vector2d &)>>
      gradient_modified;
  for (lf::base::size_type lvl = 0; lvl < refinement_level + 1; ++lvl) {
    const auto v = projects::ipdg_stokes::post_processing::extractVelocity(
        solutions[lvl].mesh, *(solutions[lvl].dofh), solutions[lvl].solution);
    const auto v_modified =
        projects::ipdg_stokes::post_processing::extractVelocity(
            solutions[lvl].mesh, *(solutions[lvl].dofh),
            solutions[lvl].solution_modified);
    velocity[lvl] = [v](const lf::mesh::Entity &entity, const Eigen::Vector2d &
                        /*unused*/) -> Eigen::Vector2d { return v(entity); };
    gradient[lvl] =
        [](const lf::mesh::Entity & /*unused*/, const Eigen::Vector2d &
           /*unused*/) -> Eigen::Matrix2d { return Eigen::Matrix2d::Zero(); };
    velocity_modified[lvl] =
        [v_modified](
            const lf::mesh::Entity &entity, const Eigen::Vector2d &
            /*unused*/) -> Eigen::Vector2d { return v_modified(entity); };
    gradient_modified[lvl] =
        [](const lf::mesh::Entity & /*unused*/, const Eigen::Vector2d &
           /*unused*/) -> Eigen::Matrix2d { return Eigen::Matrix2d::Zero(); };
  }
  auto fine_velocity =
      projects::ipdg_stokes::post_processing::bringToFinestMesh(*mesh_hierarchy,
                                                                velocity);
  auto fine_gradient =
      projects::ipdg_stokes::post_processing::bringToFinestMesh(*mesh_hierarchy,
                                                                gradient);
  auto fine_velocity_modified =
      projects::ipdg_stokes::post_processing::bringToFinestMesh(
          *mesh_hierarchy, velocity_modified);
  auto fine_gradient_modified =
      projects::ipdg_stokes::post_processing::bringToFinestMesh(
          *mesh_hierarchy, gradient_modified);

  // Perform post processing on the data
  lf::io::VtkWriter writer(solutions.back().mesh, "result.vtk");
  for (lf::base::size_type lvl = 0; lvl < refinement_level; ++lvl) {
    writer.WriteCellData(concat("v_", solutions[lvl].mesh->NumEntities(2)),
                         *lf::mesh::utils::make_LambdaMeshDataSet(
                             [&](const lf::mesh::Entity &e) {
                               return fine_velocity[lvl](
                                   e, Eigen::Vector2d::Zero());
                             }));
    writer.WriteCellData(
        concat("v_modified_", solutions[lvl].mesh->NumEntities(2)),
        *lf::mesh::utils::make_LambdaMeshDataSet(
            [&](const lf::mesh::Entity &e) {
              return fine_velocity_modified[lvl](e, Eigen::Vector2d::Zero());
            }));
    // The error in the velocity
    auto diff_v = [&](const lf::mesh::Entity &entity,
                      const Eigen::Vector2d &x) -> Eigen::Vector2d {
      return fine_velocity[lvl](entity, x) -
             fine_velocity[refinement_level](entity, x);
    };
    auto diff_v_modified = [&](const lf::mesh::Entity &entity,
                               const Eigen::Vector2d &x) -> Eigen::Vector2d {
      return fine_velocity_modified[lvl](entity, x) -
             fine_velocity_modified[refinement_level](entity, x);
    };
    // The weighted error in the velocity
    auto diff_v_weighted = [&](const lf::mesh::Entity &entity,
                               const Eigen::Vector2d &x) -> Eigen::Vector2d {
      static constexpr double delta = 0.1;
      const double weight = x[1] >= 1 - delta
                                ? (1 - std::cos(M_PI / delta * (1 - x[1]))) / 2
                                : 1.;
      return weight * diff_v(entity, x);
    };
    auto diff_v_weighted_modified =
        [&](const lf::mesh::Entity &entity,
            const Eigen::Vector2d &x) -> Eigen::Vector2d {
      static constexpr double delta = 0.1;
      const double weight = x[1] >= 1 - delta
                                ? (1 - std::cos(M_PI / delta * (1 - x[1]))) / 2
                                : 1.;
      return weight * diff_v_modified(entity, x);
    };
    // The difference in the gradient of the velocity
    auto diff_g = [&](const lf::mesh::Entity &entity,
                      const Eigen::Vector2d &x) -> Eigen::Matrix2d {
      return fine_gradient[lvl](entity, x) -
             fine_gradient[refinement_level](entity, x);
    };
    auto diff_g_modified = [&](const lf::mesh::Entity &entity,
                               const Eigen::Vector2d &x) -> Eigen::Matrix2d {
      return fine_gradient_modified[lvl](entity, x) -
             fine_gradient_modified[refinement_level](entity, x);
    };
    const double L2 = projects::ipdg_stokes::post_processing::L2norm(
        solutions[refinement_level].mesh, diff_v, 0);
    const double L2w = projects::ipdg_stokes::post_processing::L2norm(
        solutions[refinement_level].mesh, diff_v_weighted, 10);
    const double DG = projects::ipdg_stokes::post_processing::DGnorm(
        solutions[refinement_level].mesh, diff_v, diff_g, 0);
    const double L2_modified = projects::ipdg_stokes::post_processing::L2norm(
        solutions[refinement_level].mesh, diff_v_modified, 0);
    const double L2w_modified = projects::ipdg_stokes::post_processing::L2norm(
        solutions[refinement_level].mesh, diff_v_weighted_modified, 10);
    const double DG_modified = projects::ipdg_stokes::post_processing::DGnorm(
        solutions[refinement_level].mesh, diff_v_modified, diff_g_modified, 0);
    std::cout << lvl << ' ' << solutions[lvl].mesh->NumEntities(2) << ' ' << L2
              << ' ' << L2w << ' ' << DG << ' ' << L2_modified << ' '
              << L2w_modified << ' ' << DG_modified << std::endl;
  }

  return 0;
}
