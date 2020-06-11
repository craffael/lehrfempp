/** @file
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 27.05.20
 *  @copyright ETH Zurich
 */

#include <lf/assemble/dofhandler.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/mesh_function_transfer.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <boost/filesystem.hpp>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "strangsplitting.cc"

using namespace FisherKPP;

double getMeshSize(const std::shared_ptr<const lf::mesh::Mesh>& mesh_p) {
  double mesh_size = 0.0;
  /* Find maximal edge length */
  double edge_length;
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    /* Compute the length of the edge */
    auto endpoints = lf::geometry::Corners(*(edge->Geometry()));
    edge_length = (endpoints.col(0) - endpoints.col(1)).norm();
    if (mesh_size < edge_length) {
      mesh_size = edge_length;
    }
  }
  return mesh_size;
}

int main(int /*argc*/, char** /*argv*/) {
  /* Obtain mesh */
  std::unique_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  boost::filesystem::path here = __FILE__;
  auto mesh_file = (here.parent_path() / "/meshes/test1.msh").string();
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();
  /* Mesh hierarchy */
  std::unique_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory2 =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::refinement::MeshHierarchy hierarchy(mesh_p, std::move(mesh_factory2));
  /* Finite element space on finest mesh */
  std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_fine;

  Eigen::VectorXi numDofs(5);
  numDofs.setZero();
  Eigen::VectorXd meshsizes(5);
  meshsizes.setZero();
  Eigen::VectorXi m(5);
  m.setZero();
  Eigen::VectorXd tau(5);
  tau.setZero();

  /* Final time */
  double T = 1.;

  m(0) = 100;
  tau(0) = T / m(0);
  meshsizes(0) = 3.52544;

  Eigen::VectorXd mu1(39);
  mu1.setZero();
  Eigen::VectorXd mu2(125);
  mu2.setZero();
  Eigen::VectorXd mu3(444);
  mu3.setZero();
  Eigen::VectorXd mu4(1670);
  mu4.setZero();
  Eigen::VectorXd mu5(6474);
  mu5.setZero();

  Eigen::VectorXd mu1_ipol(6474);
  mu1_ipol.setZero();
  Eigen::VectorXd mu2_ipol(6474);
  mu2_ipol.setZero();
  Eigen::VectorXd mu3_ipol(6474);
  mu3_ipol.setZero();
  Eigen::VectorXd mu4_ipol(6474);
  mu4_ipol.setZero();

  /* Diffusion Coefficient */
  auto c = [](const Eigen::Vector2d & /*x*/) -> double { return 1.2; };
  /* Growth Factor */
  double lambda = 2.1;

  hierarchy.RefineRegular();
  hierarchy.RefineRegular();
  hierarchy.RefineRegular();
  hierarchy.RefineRegular();

  for (int l = 4; l >= 0; l--) {
    mesh_p = hierarchy.getMesh(l);
    /* Finite Element Space */
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    /* Dofhandler */
    const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
    const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

    numDofs(l) = N_dofs;
    meshsizes(l) = getMeshSize(mesh_p);

    std::cout << "Num of Dofs " << N_dofs << std::endl;
    std::cout << "Meshsizes " << meshsizes(l) << std::endl;

    /* Initial Population density */
    Eigen::VectorXd u0(N_dofs);
    u0.setZero();

    if (l == 0) {
      u0(13) = 0.3;
    } else if (l == 1) {
      u0(26) = 0.3;
    } else if (l == 2) {
      u0(52) = 0.3;
    } else if (l == 3) {
      u0(104) = 0.3;
    } else if (l == 4) {
      u0(208) = 0.3;
    }

    /* Non Local Boundary Conditions */

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

    /* This h is to be used as a function handle for the gain term.
     * Use it in the MassEdgeMatrix Provider.
     */
    auto h = [fe_space, mesh_p, edge_pred](const Eigen::Vector2d& x) -> double {
      double res = 0.0;

      auto g = [&x](const Eigen::Vector2d& y) -> double {
        double tmp_res = 0.0;
        if ((x - y).norm() >= 15 && (x - y).norm() <= 35) {
          tmp_res = (1.0 / (1.0 + (x - y).squaredNorm()));
        }
        return tmp_res;
      };

      auto fe = fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());
      res = localQuadFunction(
          *mesh_p,
          {{lf::base::RefEl::kSegment(),
            lf::quad::make_QuadRule(lf::base::RefEl::kSegment(),
                                    2 * fe->Degree())},
           {lf::base::RefEl::kTria(),
            lf::quad::make_TriaQR_EdgeMidpointRule()}},
          g, 1, edge_pred);
      return res;
    };

    /* In what follows, the loss term is assembled. */
    Eigen::MatrixXd L(N_dofs, N_dofs);
    for (int j = 0; j < N_dofs; j++) {
      if (boundary_nodes(j)) {
        auto L_j = [fe_space, &dofh, N_dofs, edge_pred,
                    j](const Eigen::Vector2d& x) -> double {
          auto g_j = [&x](const Eigen::Vector2d& y) -> double {
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
    std::cout << "L norm " << L.norm() << std::endl;

    /* Strang Splitting Method
     * SDIRK-2 evolution of linear parabolic term
     * Exact Evolution for nonlinear reaction term
     */

    if (l > 0) {
      tau(l) = (tau(0) / meshsizes(0)) * meshsizes(l);
      m(l) = std::round(T / tau(l));
    }

    std::cout << "tau " << tau(l) << std::endl;
    std::cout << "M " << m(l) << std::endl;

    /* Carrying Capacity */
    Eigen::VectorXd cap(N_dofs);
    cap.setZero();
    cap = 0.8 * Eigen::VectorXd::Ones(N_dofs);

    /* Now we may compute the solution */
    if (l == 0) {
      StrangSplit StrangSplitter1(fe_space, T, m(0), lambda, c, h, L);
      mu1 = StrangSplitter1.Evolution(cap, u0);
      const lf::uscalfe::MeshFunctionFE mf_coarse(fe_space, mu1);
      const lf::refinement::MeshFunctionTransfer mf_fine(hierarchy, mf_coarse,
                                                         0, 4);
      mu1_ipol = lf::uscalfe::NodalProjection(*fe_space_fine, mf_fine);
    }

    if (l == 1) {
      StrangSplit StrangSplitter2(fe_space, T, m(1), lambda, c, h, L);
      mu2 = StrangSplitter2.Evolution(cap, u0);
      const lf::uscalfe::MeshFunctionFE mf_coarse(fe_space, mu2);
      const lf::refinement::MeshFunctionTransfer mf_fine(hierarchy, mf_coarse,
                                                         1, 4);
      mu2_ipol = lf::uscalfe::NodalProjection(*fe_space_fine, mf_fine);
    }
    if (l == 2) {
      StrangSplit StrangSplitter3(fe_space, T, m(2), lambda, c, h, L);
      mu3 = StrangSplitter3.Evolution(cap, u0);
      const lf::uscalfe::MeshFunctionFE mf_coarse(fe_space, mu3);
      const lf::refinement::MeshFunctionTransfer mf_fine(hierarchy, mf_coarse,
                                                         2, 4);
      mu3_ipol = lf::uscalfe::NodalProjection(*fe_space_fine, mf_fine);
    }

    if (l == 3) {
      StrangSplit StrangSplitter4(fe_space, T, m(3), lambda, c, h, L);
      mu4 = StrangSplitter4.Evolution(cap, u0);
      const lf::uscalfe::MeshFunctionFE mf_coarse(fe_space, mu4);
      const lf::refinement::MeshFunctionTransfer mf_fine(hierarchy, mf_coarse,
                                                         3, 4);
      mu4_ipol = lf::uscalfe::NodalProjection(*fe_space_fine, mf_fine);
    }

    if (l == 4) {
      StrangSplit StrangSplitter5(fe_space, T, m(4), lambda, c, h, L);
      mu5 = StrangSplitter5.Evolution(cap, u0);
      fe_space_fine = fe_space;
    }
  }

  /* Compare the solution to the Reference solution on the finest mesh. */
  Eigen::Vector4d eL2;
  eL2.setZero();

  eL2(0) = (mu1_ipol - mu5).lpNorm<2>();
  eL2(1) = (mu2_ipol - mu5).lpNorm<2>();
  eL2(2) = (mu3_ipol - mu5).lpNorm<2>();
  eL2(3) = (mu4_ipol - mu5).lpNorm<2>();

  for (int l = 0; l < 4; l++) {
    std::cout << "numdofs for l = " << l << " : " << numDofs(l) << std::endl;
    std::cout << "meshsize for l = " << l << " : " << meshsizes(l) << std::endl;
    std::cout << "for l = " << l << "number of timesteps m : " << m(l)
              << std::endl;
    std::cout << "for l = " << l << "timestep size : " << tau(l) << std::endl;

    std::cout << "L2 error of solution with respect to reference solution on "
                 "finest mesh, for l = "
              << l << " : " << eL2(l) << std::endl;
  }

  /* Define output file format */
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream file;

  file.open("eL2.csv");
  file << eL2.format(CSVFormat);
  file.close();

  file.open("m.csv");
  file << m.format(CSVFormat);
  file.close();

  file.open("tau.csv");
  file << tau.format(CSVFormat);
  file.close();

  file.open("NumDofs.csv");
  file << numDofs.format(CSVFormat);
  file.close();

  file.open("meshsizes.csv");
  file << meshsizes.format(CSVFormat);
  file.close();

  /* SOLUTION */
  file.open("mu1.csv");
  file << mu1.format(CSVFormat);
  file.close();

  file.open("mu2.csv");
  file << mu2.format(CSVFormat);
  file.close();

  file.open("mu3.csv");
  file << mu3.format(CSVFormat);
  file.close();

  file.open("mu4.csv");
  file << mu4.format(CSVFormat);
  file.close();

  file.open("mu5.csv");
  file << mu5.format(CSVFormat);
  file.close();

  /* Interpolated Solution */
  file.open("mu1_ipol.csv");
  file << mu1_ipol.format(CSVFormat);
  file.close();

  file.open("mu2_ipol.csv");
  file << mu2_ipol.format(CSVFormat);
  file.close();

  file.open("mu3_ipol.csv");
  file << mu3_ipol.format(CSVFormat);
  file.close();

  file.open("mu4_ipol.csv");
  file << mu4_ipol.format(CSVFormat);
  file.close();

  return 0;
}
