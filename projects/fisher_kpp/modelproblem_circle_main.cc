/** @file
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 13.05.20
 *  @copyright ETH Zurich
 */

#include <lf/io/io.h>

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "strangsplitting.h"

using FisherKPP::StrangSplit;

int main(int /*argc*/, char ** /*argv*/) {
  /* Obtain mesh */
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const std::filesystem::path here = __FILE__;

  /* In case of one circle: Mesh 1. */
  auto mesh_file = (here.parent_path() / "/meshes/circle.msh").string();
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  const std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  /* Finite Element Space */
  std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  /* Dofhandler */
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  std::cout << "N_dofs :" << N_dofs << '\n';

  /* Initial Population density */
  Eigen::VectorXd u0(N_dofs);
  u0.setZero();

  /* In case of the circle: Mesh 1. For one source only: */
  u0(0) = 0.001;
  /* For two sources on the circle: */
  /* u0(0) = 0.0005;
   * u0(53) = 0.0005;
   */

  std::cout << "norm u0 " << u0.norm() << '\n';

  /* Diffusion Coefficient */
  auto c = [](const Eigen::Vector2d & /*x*/) -> double {
    /* In case of the circle: Mesh 1. */
    return 0.007;
  };

  /* Growth Factor */
  /* In case of the circle: Mesh 1. */
  const double lambda = 0.008;

  /* Boundary Conditions. */

  /* In case of the circle: Mesh 1. */
  auto h = [](const Eigen::Vector2d & /*x*/) -> double { return 0.0; };
  Eigen::MatrixXd L(N_dofs, N_dofs);
  L.setZero();

  /* Strang Splitting Method
   * SDIRK-2 evolution of linear parabolic term
   * Exact Evolution for nonlinear reaction term
   */

  /* Total number of timesteps */
  const unsigned int m = 100;
  const double T = 1.;  // the timestepsize tau will equal T/m = 0.01

  /* First we assemble the carrying capacity maps */

  Eigen::VectorXd cap(N_dofs);
  cap.setZero();
  cap = 0.8 * Eigen::VectorXd::Ones(N_dofs);

  /* Now we may compute the solution */
  StrangSplit StrangSplitter(fe_space, T, m, lambda, c, h, L);

  Eigen::MatrixXd sol(N_dofs, 20);

  sol.col(0) = StrangSplitter.Evolution(cap, u0);
  std::cout << "sol1" << '\n';
  sol.col(1) = StrangSplitter.Evolution(cap, sol.col(0));
  std::cout << "sol2" << '\n';
  sol.col(2) = StrangSplitter.Evolution(cap, sol.col(1));
  std::cout << "sol3" << '\n';
  sol.col(3) = StrangSplitter.Evolution(cap, sol.col(2));
  std::cout << "sol4" << '\n';
  sol.col(4) = StrangSplitter.Evolution(cap, sol.col(3));
  std::cout << "sol5" << '\n';

  /* Use VTK-Writer for Visualization of solution */

  for (int k = 1; k < 6; k++) {
    std::stringstream filename;
    filename << "sol" << k << ".vtk";

    lf::io::VtkWriter vtk_writer(mesh_p, filename.str());
    auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2);
    for (int global_idx = 0; global_idx < N_dofs; global_idx++) {
      nodal_data->operator()(dofh.Entity(global_idx)) =
          sol.col(k - 1)[global_idx];
    }

    vtk_writer.WritePointData("sol", *nodal_data);
  }

  return 0;
}
