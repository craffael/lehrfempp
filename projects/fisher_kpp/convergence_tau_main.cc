/** @file
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 13.05.20
 *  @copyright ETH Zurich
 */

#include <lf/io/io.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "strangsplitting.h"

using FisherKPP::localQuadFunction;
using FisherKPP::StrangSplit;

double getMeshSize(const std::shared_ptr<const lf::mesh::Mesh> &mesh_p) {
  double mesh_size = 0.0;
  /* Find maximal edge length */
  double edge_length;
  for (const lf::mesh::Entity *edge : mesh_p->Entities(1)) {
    /* Compute the length of the edge */
    auto endpoints = lf::geometry::Corners(*(edge->Geometry()));
    edge_length = (endpoints.col(0) - endpoints.col(1)).norm();
    if (mesh_size < edge_length) {
      mesh_size = edge_length;
    }
  }
  return mesh_size;
}

int main(int /*argc*/, char ** /*argv*/) {
  /* Diffusion Coefficient */
  auto c = [](const Eigen::Vector2d & /*x*/) -> double { return 1.2; };
  /* Growth Factor */
  double lambda = 2.1;

  /* Obtain mesh */
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  std::filesystem::path here = __FILE__;
  auto mesh_file = (here.parent_path() / "/meshes/test4.msh").string();
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();

  /* Finite Element Space */
  std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  /* Dofhandler */
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  std::cout << "Num of dofs " << N_dofs << std::endl;

  /* Initial Population density */
  Eigen::VectorXd u0(N_dofs);
  u0.setZero();
  u0(104) = 0.3;

  /* Non Local Boundary Conditions */

  /* This predicate returns true for nodes on the boundary */
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  auto boundary_nodes = [&bd_flags, &dofh](unsigned int idx) -> bool {
    return bd_flags(dofh.Entity(idx));
  };

  /* This predicate returns true for edges on the boundary */
  auto bd_flags_edge{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  auto edge_pred = [&bd_flags_edge](const lf::mesh::Entity &edge) -> bool {
    return bd_flags_edge(edge);
  };

  /* This h is to be used as a function handle for the gain term.
   * Use it in the MassEdgeMatrix Provider.
   */
  auto h = [fe_space, mesh_p, edge_pred](const Eigen::Vector2d &x) -> double {
    double res = 0.0;

    auto g = [&x](const Eigen::Vector2d &y) -> double {
      double tmp_res = 0.0;
      if ((x - y).norm() >= 15 && (x - y).norm() <= 35) {
        tmp_res = (1.0 / (1.0 + (x - y).squaredNorm()));
      }
      return tmp_res;
    };

    const auto *fe = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
    res = localQuadFunction(
        *mesh_p,
        {{lf::base::RefEl::kSegment(),
          lf::quad::make_QuadRule(lf::base::RefEl::kSegment(),
                                  2 * fe->Degree())},
         {lf::base::RefEl::kTria(), lf::quad::make_TriaQR_EdgeMidpointRule()}},
        g, 1, edge_pred);
    return res;
  };

  /* In what follows, the loss term is assembled. */
  Eigen::MatrixXd L(N_dofs, N_dofs);
  for (int j = 0; j < N_dofs; j++) {
    if (boundary_nodes(j)) {
      auto L_j = [fe_space, &dofh, N_dofs, edge_pred,
                  j](const Eigen::Vector2d &x) -> double {
        auto g_j = [&x](const Eigen::Vector2d &y) -> double {
          double tmp_res = 0.0;
          if ((x - y).norm() >= 15 && (x - y).norm() <= 35) {
            tmp_res = (1.0 / (1.0 + (x - y).squaredNorm()));
          }
          return tmp_res;
        };

        Eigen::VectorXd L1(N_dofs);
        L1.setZero();

        lf::mesh::utils::MeshFunctionGlobal mf_g{g_j};
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
  Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> numSteps(6);
  numSteps.setZero();
  Eigen::VectorXd tau(6);
  tau.setZero();
  double T = 1.;

  for (int i = 0; i < 6; i++) {
    numSteps(i) = static_cast<unsigned int>(80 * std::pow(2, i));
    tau(i) = T / numSteps(i);
  }

  /* Carrying Capacity. */
  Eigen::VectorXd cap(N_dofs);
  cap.setZero();
  cap = 0.8 * Eigen::VectorXd::Ones(N_dofs);

  StrangSplit StrangSplitter1(fe_space, T, numSteps(0), lambda, c, h, L);
  StrangSplit StrangSplitter2(fe_space, T, numSteps(1), lambda, c, h, L);
  StrangSplit StrangSplitter3(fe_space, T, numSteps(2), lambda, c, h, L);
  StrangSplit StrangSplitter4(fe_space, T, numSteps(3), lambda, c, h, L);
  StrangSplit StrangSplitter5(fe_space, T, numSteps(4), lambda, c, h, L);
  StrangSplit StrangSplitter6(fe_space, T, numSteps(5), lambda, c, h, L);

  Eigen::VectorXd sol1(N_dofs);
  sol1.setZero();
  sol1 = StrangSplitter1.Evolution(cap, u0);
  std::cout << "sol1 " << std::endl;

  Eigen::VectorXd sol2(N_dofs);
  sol2.setZero();
  sol2 = StrangSplitter2.Evolution(cap, u0);
  std::cout << "sol2 " << std::endl;

  Eigen::VectorXd sol3(N_dofs);
  sol3.setZero();
  sol3 = StrangSplitter3.Evolution(cap, u0);
  std::cout << "sol3 " << std::endl;

  Eigen::VectorXd sol4(N_dofs);
  sol4.setZero();
  sol4 = StrangSplitter4.Evolution(cap, u0);
  std::cout << "sol4 " << std::endl;

  Eigen::VectorXd sol5(N_dofs);
  sol5.setZero();
  sol5 = StrangSplitter5.Evolution(cap, u0);
  std::cout << "sol5 " << std::endl;

  Eigen::VectorXd sol6(N_dofs);
  sol6.setZero();
  sol6 = StrangSplitter6.Evolution(cap, u0);
  std::cout << "sol6 " << std::endl;

  Eigen::VectorXd eL2(5);
  eL2.setZero();
  eL2(0) = (sol6 - sol1).lpNorm<2>();
  eL2(1) = (sol6 - sol2).lpNorm<2>();
  eL2(2) = (sol6 - sol3).lpNorm<2>();
  eL2(3) = (sol6 - sol4).lpNorm<2>();
  eL2(4) = (sol6 - sol5).lpNorm<2>();
  std::cout << "errors computed" << std::endl;

  for (int l = 0; l < 6; l++) {
    if (l < 5) {
      std::cout << "errorL2 solcol0 for l = " << l << " : " << eL2(l)
                << std::endl;
    }

    std::cout << "number of timesteps for l = " << l << " : " << numSteps(l)
              << std::endl;
    std::cout << "tau for l = " << l << " : " << tau(l) << std::endl;
  }

  /* Define output file format */
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");

  std::ofstream file;

  file.open("eL2.csv");
  file << eL2.format(CSVFormat);
  file.close();

  file.open("tau.csv");
  file << tau.format(CSVFormat);
  file.close();

  file.open("numSteps.csv");
  file << numSteps.format(CSVFormat);
  file.close();

  file.open("sol1.csv");
  file << sol1.format(CSVFormat);
  file.close();

  file.open("sol2.csv");
  file << sol2.format(CSVFormat);
  file.close();

  file.open("sol3.csv");
  file << sol3.format(CSVFormat);
  file.close();

  file.open("sol4.csv");
  file << sol4.format(CSVFormat);
  file.close();

  file.open("sol5.csv");
  file << sol5.format(CSVFormat);
  file.close();

  file.open("sol6.csv");
  file << sol6.format(CSVFormat);
  file.close();

  return 0;
}
