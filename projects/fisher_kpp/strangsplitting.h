/** @file
 *  @brief Bachelor Thesis Fisher/KPP
 *  @author Amélie Justine Loher
 *  @date 01.04.20
 *  @copyright ETH Zurich
 */

#ifndef STRANG_SPLIT_H
#define STRANG_SPLIT_H

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>
#include <stdio.h>
#include <stdlib.h>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <memory>
#include <utility>

namespace FisherKPP {

class StrangSplit {
 public:
  /* Disabled constructors */
  StrangSplit() = delete;
  StrangSplit(const StrangSplit &) = delete;
  StrangSplit(StrangSplit &&) = delete;
  StrangSplit &operator=(const StrangSplit &) = delete;
  StrangSplit &operator=(const StrangSplit &&) = delete;
  /* Main constructor */
  template <typename DIFF_COEFF, typename NONLOC_BC>
  explicit StrangSplit(
      std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> &fe_space,
      double T, unsigned int m, double lambda, DIFF_COEFF &&c, NONLOC_BC &&h,
      const Eigen::MatrixXd &L);
  /* Destructor */
  virtual ~StrangSplit() = default;

  /* Member Functions */
  Eigen::VectorXd diffusionEvolutionOperator(bool firstcall,
                                             const Eigen::VectorXd &mu);
  Eigen::VectorXd Evolution(const Eigen::VectorXd &cap,
                            const Eigen::VectorXd &mu);

 private:
  /* Finite Element Space */
  const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> fe_space_;
  /* Final Time */
  double T_;
  /* Number of Timesteps */
  unsigned int m_;
  /* Time step size */
  double tau_;
  /* Growth Factor */
  double lambda_;
  /* coefficient for SDIRK-2 Butcher Tableau */
  double kappa_;
  /* Galerkin matrix corresponding to the negative Laplacian with Robin Boundary
   * Conditions */
  Eigen::SparseMatrix<double> A_;
  /* Galerkin matrix for the Mass Matrix */
  Eigen::SparseMatrix<double> M_;
  /* Precompute LU decomposition needed for time stepping */
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver1;
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver2;
};

template <typename FUNCTOR>
auto localQuadFunction(
    const lf::mesh::Mesh &mesh,
    std::map<lf::base::RefEl, lf::quad::QuadRule> quadrules, FUNCTOR &&f,
    unsigned int codim,
    const std::function<bool(const lf::mesh::Entity &)> &pred) {
  LF_ASSERT_MSG(mesh.DimMesh() >= codim, "Illegal codim = " << codim);
  /* Variable for summing the result */
  using value_t = std::invoke_result_t<FUNCTOR, Eigen::VectorXd>;
  value_t sum_var{};
  /* Loop over entities of co-dimension codim */
  for (const lf::mesh::Entity *entity : mesh.Entities(codim)) {
    if (pred(*entity)) {
      /* Obtain geometry information for entity */
      const lf::geometry::Geometry &geo{*entity->Geometry()};
      /* obtain quadrature rule suitable for entity type */
      auto tmp = quadrules.find(entity->RefEl());
      if (tmp != quadrules.end()) {
        /* A quadrature rule has been found */
        const lf::quad::QuadRule &qr{tmp->second};
        /* Number of quadrature points */
        const lf::base::size_type P = qr.NumPoints();
        /* Quadrature points */
        const Eigen::MatrixXd &zeta_ref{qr.Points()};
        /* Map quadrature points to physical/world coordinates */
        const Eigen::MatrixXd zeta{geo.Global(zeta_ref)};
        /* Quadrature weights */
        const Eigen::VectorXd &w_ref{qr.Weights()};
        /* Gramian determinants */
        const Eigen::VectorXd gram_dets{geo.IntegrationElement(zeta_ref)};
        /* Iterate over the quadrature points */
        for (int l = 0; l < P; ++l) {
          sum_var += w_ref[l] * f(zeta.col(l)) * gram_dets[l];
        }
      } else {
        LF_VERIFY_MSG(false, "Missing quadrature rule for " << entity->RefEl());
      }
    }
  }
  return sum_var;
}

/* Function for the assembly of both Galerkin Matrices, the Mass matrix and the
 * Stiffness matrix */
template <typename DIFF_COEFF, typename NONLOC_BC>
std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
assembleGalerkinMatrices(const lf::assemble::DofHandler &dofh, DIFF_COEFF &&c,
                         NONLOC_BC &&h, const Eigen::MatrixXd &L) {
  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> A_M;

  std::shared_ptr<const lf::mesh::Mesh> mesh = dofh.Mesh();
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1)};
  auto edges_predicate = [&bd_flags](const lf::mesh::Entity &edge) -> bool {
    return bd_flags(edge);
  };

  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  lf::assemble::COOMatrix<double> M_COO(N_dofs, N_dofs);
  lf::assemble::COOMatrix<double> G_COO(N_dofs, N_dofs);

  auto one = [](const Eigen::Vector2d & /*x*/) -> double { return 1.0; };
  auto zero = [](const Eigen::Vector2d & /*x*/) -> double { return 0.0; };

  lf::mesh::utils::MeshFunctionGlobal mf_c{c};
  lf::mesh::utils::MeshFunctionGlobal mf_h{h};
  lf::mesh::utils::MeshFunctionGlobal mf_one{one};
  lf::mesh::utils::MeshFunctionGlobal mf_zero{zero};

  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_c),
                                                      decltype(mf_zero)>
      elMat_Stiff(fe_space, mf_c, mf_zero);
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(mf_zero),
                                                      decltype(mf_one)>
      elMat_Mass(fe_space, mf_zero, mf_one);
  lf::uscalfe::MassEdgeMatrixProvider<double, decltype(mf_h),
                                      decltype(edges_predicate)>
      elMat_Edge(fe_space, mf_h, edges_predicate);

  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elMat_Stiff, A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elMat_Mass, M_COO);
  lf::assemble::AssembleMatrixLocally(1, dofh, dofh, elMat_Edge, A_COO);

  Eigen::SparseMatrix<double> A_sps = A_COO.makeSparse();
  Eigen::SparseMatrix<double> M(N_dofs, N_dofs);
  M = M_COO.makeSparse();

  Eigen::MatrixXd A_ds(A_sps);
  Eigen::MatrixXd A_ds1 = A_ds - L;

  Eigen::SparseMatrix<double> A(N_dofs, N_dofs);
  A = A_ds1.sparseView();

  A_M = std::make_pair(A, M);

  return A_M;
}

/* Constructor for StrangSplit */
template <typename DIFF_COEFF, typename NONLOC_BC>
StrangSplit::StrangSplit(
    std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> &fe_space,
    double T, unsigned m, double lambda, DIFF_COEFF &&c, NONLOC_BC &&h,
    const Eigen::MatrixXd &L)
    : fe_space_(fe_space), T_(T), m_(m), lambda_(lambda) {
  const lf::assemble::DofHandler &dofh{fe_space_->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
      galerkinpair = assembleGalerkinMatrices(dofh, c, h, L);
  A_ = galerkinpair.first;
  M_ = galerkinpair.second;

  /* Butcher Tableau Coefficient for SDIRK-2 */
  kappa_ = 1.0 - 0.5 * sqrt(2.0);

  tau_ = T_ / m_;
  solver1.compute(M_ + 0.5 * tau_ * kappa_ * A_);
  LF_VERIFY_MSG(solver1.info() == Eigen::Success, "LU decomposition failed");
  solver2.compute(M_ + tau_ * kappa_ * A_);
  LF_VERIFY_MSG(solver2.info() == Eigen::Success, "LU decomposition failed");
}

} /* namespace FisherKPP */

#endif
