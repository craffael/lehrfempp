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

using FisherKPP::localQuadFunction;
using FisherKPP::StrangSplit;

int main(int /*argc*/, char** /*argv*/) {
  /* Obtain mesh */
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  std::filesystem::path here = __FILE__;

  /* In case of three circles: Mesh 2. */
  auto mesh_file = (here.parent_path() / "/meshes/threecircles.msh").string();

  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = reader.mesh();
  /* Finite Element Space */
  std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  /* Dofhandler */
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  std::cout << "N_dofs :" << N_dofs << std::endl;

  /* Initial Population density */
  Eigen::VectorXd u0(N_dofs);
  u0.setZero();
  /* In case of three circles: Mesh 2. */
  u0(0) = 0.00008;

  /* Diffusion Coefficient */
  auto c = [](const Eigen::Vector2d& /*x*/) -> double {
    /* In case of three circles: Mesh 2. */
    return 0.004;
  };

  /* Growth Factor */
  /* In case of three circles: Mesh 2. */
  double lambda = 0.0042;

  /* Boundary Conditions. */

  /* In case of homogeneous Neumann boundary conditions for Mesh 2. */
  /* auto h = [] (const Eigen::Vector2d& x) -> double { return 0.0;};
   * Eigen::MatrixXd L(N_dofs, N_dofs); L.setZero();
   */

  /* Else in case of non-local boundary condition for Mesh 2. */

  /* This predicate returns true for nodes on the boundary */
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  auto boundary_nodes = [&bd_flags, &dofh](unsigned int idx) -> bool {
    return bd_flags(dofh.Entity(idx));
  };

  /* This predicate returns true for edges on the boundary */
  auto bd_flags_edge{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  auto edge_pred = [&bd_flags_edge](const lf::mesh::Entity& edge) -> bool {
    return bd_flags_edge(edge);
  };

  auto h = [fe_space, mesh_p, edge_pred](const Eigen::Vector2d& x) -> double {
    double res = 0.0;

    /* The kernel g_1: No bounds on \| x - y \|. */
    auto g_1 = [&x](const Eigen::Vector2d& y) -> double {
      double tmp_res = 0.0;
      tmp_res = (1.0 / (1.0 + (x - y).squaredNorm()));
      return tmp_res;
    };

    /* The kernel g_2: Only an upper bound. */
    auto g_2 = [&x](const Eigen::Vector2d& y) -> double {
      double tmp_res = 0.0;
      if ((x - y).norm() <= 1.75) {
        tmp_res = (1.0 / (1.0 + (x - y).squaredNorm()));
      }
      return tmp_res;
    };

    /* The kernel g_3: Lower and Upper bound. */
    auto g_3 = [&x](const Eigen::Vector2d& y) -> double {
      double tmp_res = 0.0;
      if (1.25 <= (x - y).norm() && (x - y).norm() <= 1.75) {
        tmp_res = (1.0 / (1.0 + (x - y).squaredNorm()));
      }
      return tmp_res;
    };

    const auto* fe = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
    res = localQuadFunction(
        *mesh_p,
        {{lf::base::RefEl::kSegment(),
          lf::quad::make_QuadRule(lf::base::RefEl::kSegment(),
                                  2 * fe->Degree())},
         {lf::base::RefEl::kTria(), lf::quad::make_TriaQR_EdgeMidpointRule()}},
        g_1, 1, edge_pred);

    return res;
  };

  Eigen::MatrixXd L(N_dofs, N_dofs);
  for (int j = 0; j < N_dofs; j++) {
    if (boundary_nodes(j)) {
      auto L_j = [fe_space, &dofh, N_dofs, edge_pred,
                  j](const Eigen::Vector2d& x) -> double {
        /* The kernel g_1: No bounds on \| x - y \|. */
        auto g_1_j = [&x](const Eigen::Vector2d& y) -> double {
          double tmp_res = 0.0;
          tmp_res = (1.0 / (1.0 + (x - y).squaredNorm()));
          return tmp_res;
        };

        /* The kernel g_2: Only an upper bound. */
        auto g_2_j = [&x](const Eigen::Vector2d& y) -> double {
          double tmp_res = 0.0;
          if ((x - y).norm() <= 1.75) {
            tmp_res = (1.0 / (1.0 + (x - y).squaredNorm()));
          }
          return tmp_res;
        };

        /* The kernel g_3: Lower and Upper bound. */
        auto g_3_j = [&x](const Eigen::Vector2d& y) -> double {
          double tmp_res = 0.0;
          if (1.25 <= (x - y).norm() && (x - y).norm() <= 1.75) {
            tmp_res = (1.0 / (1.0 + (x - y).squaredNorm()));
          }
          return tmp_res;
        };

        Eigen::VectorXd L1(N_dofs);
        L1.setZero();

        lf::mesh::utils::MeshFunctionGlobal mf_g{g_1_j};
        lf::uscalfe::ScalarLoadEdgeVectorProvider<double, decltype(mf_g),
                                                  decltype(edge_pred)>
            edgeVec_y(fe_space, mf_g, edge_pred);
        lf::assemble::AssembleVectorLocally(1, dofh, edgeVec_y, L1);

        return L1(j);
      };

      Eigen::VectorXd L2(N_dofs);
      L2.setZero();

      lf::mesh::utils::MeshFunctionGlobal mf_L{L_j};
      lf::uscalfe::ScalarLoadEdgeVectorProvider<double, decltype(mf_L),
                                                decltype(edge_pred)>
          edgeVec_x(fe_space, mf_L, edge_pred);
      lf::assemble::AssembleVectorLocally(1, dofh, edgeVec_x, L2);

      L.col(j) = L2;
    }

    else {
      L.col(j) = Eigen::VectorXd::Zero(N_dofs);
    }
  }

  /* Strang Splitting Method
   * SDIRK-2 evolution of linear parabolic term
   * Exact Evolution for nonlinear reaction term
   */

  /* Total number of timesteps */
  unsigned int m = 100;
  double T = 1.;  // the timestepsize tau will equal T/m = 0.01

  /* First we assemble the carrying capacity maps */

  Eigen::VectorXd cap(N_dofs);
  cap.setZero();
  cap = 0.8 * Eigen::VectorXd::Ones(N_dofs);

  /* Now we may compute the solution */
  StrangSplit StrangSplitter(fe_space, T, m, lambda, c, h, L);

  Eigen::MatrixXd sol(N_dofs, 20);

  sol.col(0) = StrangSplitter.Evolution(cap, u0);
  std::cout << "sol1" << std::endl;
  sol.col(1) = StrangSplitter.Evolution(cap, sol.col(0));
  std::cout << "sol2" << std::endl;
  sol.col(2) = StrangSplitter.Evolution(cap, sol.col(1));
  std::cout << "sol3" << std::endl;
  sol.col(3) = StrangSplitter.Evolution(cap, sol.col(2));
  std::cout << "sol4" << std::endl;
  sol.col(4) = StrangSplitter.Evolution(cap, sol.col(3));
  std::cout << "sol5" << std::endl;

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
