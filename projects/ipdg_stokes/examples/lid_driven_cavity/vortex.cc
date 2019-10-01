/**
 * @file vortex.cc
 * @brief Solve the lid driven cavity experiment on a nonuniform mesh
 */

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <boost/filesystem.hpp>
#include <cassert>

#include <lf/assemble/dofhandler.h>
#include <lf/io/gmsh_reader.h>
#include <lf/io/vtk_writer.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/utils/tp_triag_mesh_builder.h>
#include <lf/quad/quad.h>
#include <lf/refinement/refinement.h>

#include <build_system_matrix.h>
#include <piecewise_const_element_matrix_provider.h>
#include <piecewise_const_element_vector_provider.h>
#include <solution_to_mesh_data_set.h>

/**
 * @brief Solves the lid driven cavity problem on a domain [0,100]x[0,100]
 * @param mesh A shared pointer to the mesh on which to solve the PDE
 * @param dofh The dofhandler to use for the simulation
 * @param modified If true, use the modified penalty term
 * otherwise use the original one
 * @returns A vector containing the basis function coefficients of the solution
 */
Eigen::VectorXd solveLidDrivenCavity(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh, bool modified = false) {
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
    v << 1. / 100, 0;
    if (vertices(1, 0) <= 100 + eps && vertices(1, 0) >= 100 - eps &&
        vertices(1, 1) <= 100 + eps && vertices(1, 1) >= 100 - eps) {
      return v;
    }
    return Eigen::Vector2d::Zero();
  };

  // Solve the LSE using sparse cholesky
  const auto [A, rhs] =
      projects::ipdg_stokes::assemble::buildSystemMatrixNoFlow(
          mesh, dofh, f, dirichlet_funct, 1,
          lf::quad::make_TriaQR_MidpointRule(), modified);
  Eigen::SparseMatrix<double> As = A.makeSparse();
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(As);
  return solver.solve(rhs);
}

/**
 * @brief stores the solution of the lid driven cavity experiment to a vtk file
 */
int main() {
  const double mu = 1;
  const double sigma = 1;
  const double rho = 1;

  // Load the mesh
  boost::filesystem::path meshpath = __FILE__;
  meshpath = meshpath.parent_path() / "mesh.msh";
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(factory), meshpath.c_str());
  auto mesh = reader.mesh();

  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh);

  lf::assemble::UniformFEDofHandler dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1}, {lf::base::RefEl::kSegment(), 1}});
  std::cout << "solving original" << std::endl;
  const Eigen::VectorXd solution_original =
      solveLidDrivenCavity(mesh, dofh, false);
  std::cout << "extracting original" << std::endl;
  const auto c_original =
      projects::ipdg_stokes::post_processing::extractBasisFunctionCoefficients(
          mesh, dofh, solution_original);
  const auto v_original =
      projects::ipdg_stokes::post_processing::extractVelocity(
          mesh, dofh, solution_original);
  std::cout << "solving modified" << std::endl;
  const Eigen::VectorXd solution_modified =
      solveLidDrivenCavity(mesh, dofh, true);
  std::cout << "extracting modified" << std::endl;
  const auto c_modified =
      projects::ipdg_stokes::post_processing::extractBasisFunctionCoefficients(
          mesh, dofh, solution_modified);
  const auto v_modified =
      projects::ipdg_stokes::post_processing::extractVelocity(
          mesh, dofh, solution_modified);

  std::cout << "writing" << std::endl;
  lf::io::VtkWriter writer(mesh, "vortex.vtk");
  writer.WritePointData("coefficients_original", c_original);
  writer.WritePointData("coefficients_modified", c_modified);
  writer.WriteCellData("velocity_original", v_original);
  writer.WriteCellData("velocity_modified", v_modified);

  return 0;
}
