/**
 * @file manufactured_solution.cc
 * @brief Determines the h-convergence of the manufactured solution
 */

#define _USE_MATH_DEFINES

#include <boost/filesystem.hpp>
#include <cmath>
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
using lf::uscalfe::operator*;

/**
 * @brief Compute the analytic flow velocity
 * @param n Number of rings around the origin
 * @param x The position at which to get the velocity
 * @returns The velocity at position x
 */
Eigen::Vector2d computeU(int n, Eigen::Vector2d x) {
  const double r = x.norm();
  const double u0 = 1 - std::cos(2 * M_PI * n * r);
  Eigen::Vector2d u;
  u << -x[1], x[0];
  return u0 * u.normalized();
}

/**
 * @brief Compute the analytic gradient of the flow velocity
 * @param n Number of rings around the origin
 * @param x The position at which to get the gradient
 * @returns The gradient of the velocity at position x
 */
Eigen::Matrix2d computeUGrad(int n, Eigen::Vector2d x) {
  const double r = x.norm();
  const double r2 = r * r;
  const double r3 = r2 * r;
  const double phi = 2 * M_PI * n * r;
  const double cphi = std::cos(phi);
  const double sphi = std::sin(phi);
  Eigen::Matrix2d grad;
  grad << x[0] * x[1] / r3 * (1 - cphi) -
              2 * M_PI * n * x[0] * x[1] / r2 * sphi,
      (1. / r - x[0] * x[0] / r3) * (1 - cphi) +
          2 * M_PI * n * x[0] * x[0] / r2 * sphi,
      (x[0] * x[1] / r3 - 1. / r) * (1 - cphi) +
          2 * M_PI * n * x[1] * x[1] / r2 * sphi,
      -x[0] * x[1] / r3 * (1 - cphi) + 2 * M_PI * n * x[0] * x[1] / r2 * sphi;
  return grad;
}

/**
 * @brief Compute the analytic volumetric forces
 * @param n Number of rings around the origin
 * @param x The position at which to get the force
 * @returns The volumetric forces at position x
 */
Eigen::Vector2d computeF(int n, Eigen::Vector2d x) {
  const double r = x.norm();
  const double f0 = r == 0 ? -6 * n * n * M_PI * M_PI
                           : (1 -
                              (1 + 4 * n * n * M_PI * M_PI * r * r) *
                                  std::cos(2 * M_PI * n * r) -
                              2 * M_PI * n * r * std::sin(2 * M_PI * n * r)) /
                                 (r * r);
  Eigen::Vector2d f;
  f << -x[1], x[0];
  return f0 * f.normalized();
}

/**
 * @brief Concatenate objects defining an operator<<(std::ostream&)
 * @param args A variadic pack of objects implementing
 * `operator<<(std::ostream&)`
 * @returns A string with the objects concatenated
 */
template <typename... Args>
static std::string concat(Args&&... args) {
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
 * @brief Outputs the L2 and DG norm errors for the manufactured solution
 *
 * The manufactured solution is solved on a sequence of `mesh_count` meshes
 * loaded from the files disk??.msh. The convergence of the solution is printed
 * to stdout. Furthermore the convergence of the scaled solution to minimize the
 * L2-norm is also printed to stdout.
 */
int main() {
  const unsigned mesh_count = 6;
  const int n = 3;
  const unsigned quadrule_degree = 10;

  // Compute the solution for each mesh
  std::vector<ProblemSolution> solutions(mesh_count);
  for (lf::base::size_type lvl = 0; lvl < mesh_count; ++lvl) {
    boost::filesystem::path meshpath = __FILE__;
    meshpath = meshpath.parent_path() / concat("disk", lvl, ".msh");
    std::unique_ptr<lf::mesh::MeshFactory> factory =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
    lf::io::GmshReader reader(std::move(factory), meshpath.string());
    auto mesh = reader.mesh();
    auto& sol = solutions[lvl];
    sol.mesh = mesh;
    sol.dofh = std::shared_ptr<const lf::assemble::DofHandler>(
        new lf::assemble::UniformFEDofHandler(
            mesh, {{lf::base::RefEl::kPoint(), 1},
                   {lf::base::RefEl::kSegment(), 1}}));
    // Use the volume forces returned by computeF
    auto f = [n](const Eigen::Vector2d& x) -> Eigen::Vector2d {
      return computeF(n, x);
    };
    // The velocity on the boundary is zero everywhere
    auto dirichlet_funct = [](const lf::mesh::Entity &
                              /*unused*/) -> Eigen::Vector2d {
      return Eigen::Vector2d::Zero();
    };
    const auto [A, rhs] =
        projects::ipdg_stokes::assemble::buildSystemMatrixNoFlow(
            sol.mesh, *(sol.dofh), f, dirichlet_funct, 1,
            lf::quad::make_QuadRule(lf::base::RefEl::kTria(), quadrule_degree),
            false);
    sol.A = A.makeSparse();
    sol.rhs = rhs;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(sol.A);
    sol.solution = solver.solve(rhs);
    const auto [A_modified, rhs_modified] =
        projects::ipdg_stokes::assemble::buildSystemMatrixNoFlow(
            sol.mesh, *(sol.dofh), f, dirichlet_funct, 1,
            lf::quad::make_QuadRule(lf::base::RefEl::kTria(), quadrule_degree),
            true);
    sol.A_modified = A_modified.makeSparse();
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_modified(
        sol.A_modified);
    sol.solution_modified = solver_modified.solve(rhs_modified);
  }

  // Perform post processing on the data
  const auto analyticU = lf::mesh::utils::MeshFunctionGlobal(
      [&](const Eigen::Vector2d& x) -> Eigen::Vector2d {
        return computeU(n, x);
      });
  const auto analyticF = lf::mesh::utils::MeshFunctionGlobal(
      [&](const Eigen::Vector2d& x) -> Eigen::Vector2d {
        return computeF(n, x);
      });
  for (lf::base::size_type lvl = 0; lvl < mesh_count; ++lvl) {
    const auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(
            solutions[lvl].mesh);
    const auto velocity_exact = lf::mesh::utils::MeshFunctionGlobal(
        [n](const Eigen::Vector2d& x) { return computeU(n, x); });
    const auto grad_exact = lf::mesh::utils::MeshFunctionGlobal(
        [n](const Eigen::Vector2d& x) -> Eigen::Matrix2d {
          return computeUGrad(n, x);
        });
    projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double, double>
        velocity(fe_space, solutions[lvl].solution);
    projects::ipdg_stokes::post_processing::MeshFunctionVelocity<double, double>
        velocity_modified(fe_space, solutions[lvl].solution_modified);
    lf::io::VtkWriter writer(solutions[lvl].mesh,
                             concat("result", lvl, ".vtk"));
    writer.WriteCellData("analyticU", analyticU);
    writer.WriteCellData("analyticF", analyticF);
    writer.WriteCellData("U", velocity);
    writer.WriteCellData("U_modified", velocity_modified);
    // The error in the velocity
    const auto diff_v = velocity - velocity_exact;
    const auto diff_v_modified = velocity_modified - velocity_exact;
    // The error in the gradient of the velocity
    const auto diff_grad = -grad_exact;
    const auto diff_grad_modified = -grad_exact;
    // Approximately compute the factor the numerical solution is off by
    const auto qr_provider = [](const lf::mesh::Entity& e) {
      return lf::quad::make_QuadRule(e.RefEl(), 0);
    };
    const double factor = lf::uscalfe::IntegrateMeshFunction(
                              *(solutions[lvl].mesh),
                              lf::uscalfe::transpose(velocity) * velocity_exact,
                              qr_provider)(0, 0) /
                          lf::uscalfe::IntegrateMeshFunction(
                              *(solutions[lvl].mesh),
                              lf::uscalfe::squaredNorm(velocity), qr_provider);
    const double factor_modified =
        lf::uscalfe::IntegrateMeshFunction(
            *(solutions[lvl].mesh),
            lf::uscalfe::transpose(velocity_modified) * velocity_exact,
            qr_provider)(0, 0) /
        lf::uscalfe::IntegrateMeshFunction(
            *(solutions[lvl].mesh), lf::uscalfe::squaredNorm(velocity_modified),
            qr_provider);
    std::cout << factor << std::endl;
    // The error in the corrected velocity
    const auto velocity_scaled =
        lf::mesh::utils::MeshFunctionConstant<double>(factor) * velocity;
    const auto velocity_scaled_modified =
        lf::mesh::utils::MeshFunctionConstant(factor_modified) *
        velocity_modified;
    const auto diff_v_fac = velocity_scaled - velocity_exact;
    const auto diff_v_fac_modified = velocity_scaled_modified - velocity_exact;
    // The error in the gradient of the corrected velocty
    const auto& diff_g_fac = diff_grad;
    const auto& diff_g_fac_modified = diff_grad_modified;
    const double L2 = projects::ipdg_stokes::post_processing::L2norm(
        solutions[lvl].mesh, diff_v, qr_provider);
    const double DG = projects::ipdg_stokes::post_processing::DGnorm(
        solutions[lvl].mesh, diff_v, diff_grad, qr_provider);
    const double L2f = projects::ipdg_stokes::post_processing::L2norm(
        solutions[lvl].mesh, diff_v_fac, qr_provider);
    const double DGf = projects::ipdg_stokes::post_processing::DGnorm(
        solutions[lvl].mesh, diff_v_fac, diff_g_fac, qr_provider);
    const double L2_modified = projects::ipdg_stokes::post_processing::L2norm(
        solutions[lvl].mesh, diff_v_modified, qr_provider);
    const double DG_modified = projects::ipdg_stokes::post_processing::DGnorm(
        solutions[lvl].mesh, diff_v_modified, diff_grad_modified, qr_provider);
    const double L2f_modified = projects::ipdg_stokes::post_processing::L2norm(
        solutions[lvl].mesh, diff_v_fac_modified, qr_provider);
    const double DGf_modified = projects::ipdg_stokes::post_processing::DGnorm(
        solutions[lvl].mesh, diff_v_fac_modified, diff_g_fac_modified,
        qr_provider);
    std::cout << lvl << ' ' << solutions[lvl].mesh->NumEntities(2) << ' ' << L2
              << ' ' << DG << ' ' << L2f << ' ' << DGf << ' ' << L2_modified
              << ' ' << DG_modified << ' ' << L2f_modified << ' '
              << DGf_modified << std::endl;
  }
}
