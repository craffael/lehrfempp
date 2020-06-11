/** @file
 *	@brief Bachelor Thesis Fisher/KPP
 *	@author Amélie Justine Loher
 *	@date 01.04.20
 *  @copyright ETH Zurich
 */

#include "strangsplitting.h"

namespace FisherKPP {

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
    : fe_space_(fe_space), T_(T), m_(m), lambda_(lambda){

  const lf::assemble::DofHandler &dofh{fe_space_->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>
      galerkinpair = assembleGalerkinMatrices(dofh, c, h, L); 
  A_ = galerkinpair.first;
  M_ = galerkinpair.second;

  /* Butcher Tableau Coefficient for SDIRK-2 */
  kappa_ = 1.0 - 0.5 * sqrt(2.0);
  
  tau_ = T_/m_;
  solver1.compute(M_ + 0.5 * tau_ * kappa_ * A_);
  LF_VERIFY_MSG(solver1.info() == Eigen::Success, "LU decomposition failed");
  solver2.compute(M_ + tau_ * kappa_ * A_);
  LF_VERIFY_MSG(solver2.info() == Eigen::Success, "LU decomposition failed");
}

/* Member Function StrangSplit
 * Computes the Evolution Operator for the linear parabolic diffusion term
 */
Eigen::VectorXd StrangSplit::diffusionEvolutionOperator(bool firstcall, const Eigen::VectorXd &mu) {
  Eigen::VectorXd evol_op;
  
  /* Precomputation */
/*  solver.compute(M_ + tau * kappa_ * A_);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
 */
  Eigen::VectorXd rhs = -A_ * mu;

  if(firstcall) {
	/* First stage SDIRK-2 */
  	Eigen::VectorXd k1 = solver1.solve(rhs);
  	LF_VERIFY_MSG(solver1.info() == Eigen::Success, "LU decomposition failed");
  	/* Second stage SDIRK-2 */
  	Eigen::VectorXd k2 = solver1.solve(rhs - 0.5 * tau_ * (1 - kappa_) * A_ * k1);
  	LF_VERIFY_MSG(solver1.info() == Eigen::Success, "LU decomposition failed");

  	evol_op = mu + 0.5 * tau_ * (1 - kappa_) * k1 + 0.5 * tau_ * kappa_ * k2;
  
  } else {
	 /* First stage SDIRK-2 */
    Eigen::VectorXd k1 = solver2.solve(rhs);
    LF_VERIFY_MSG(solver2.info() == Eigen::Success, "LU decomposition failed");
    /* Second stage SDIRK-2 */
    Eigen::VectorXd k2 = solver2.solve(rhs - tau_ * (1 - kappa_) * A_ * k1);
    LF_VERIFY_MSG(solver2.info() == Eigen::Success, "LU decomposition failed");

    evol_op = mu + tau_ * (1 - kappa_) * k1 + tau_ * kappa_ * k2;
  }
  
  return evol_op;
}

/* Member Function StrangSplit
 * Computes the Evolution for m_ timesteps
 */
Eigen::VectorXd StrangSplit::Evolution(const Eigen::VectorXd &cap,
                                       const Eigen::VectorXd &mu) {
  const lf::assemble::DofHandler &dofh{fe_space_->LocGlobMap()};
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  Eigen::VectorXd sol(N_dofs);
  Eigen::VectorXd sol_cur(N_dofs);
  Eigen::VectorXd sol_next(N_dofs);

  Eigen::VectorXd ones = Eigen::VectorXd::Ones(N_dofs);

  /* Time Steps */
  //double tau = T_ / m_;

  /* Strang Splitting Method
   * First half time step [0, tau/2]: Diffusion
   */
  bool firstcall = true;
  sol_cur = StrangSplit::diffusionEvolutionOperator(firstcall, mu);
  firstcall = false;
  /* loop through time steps */
  for (unsigned int i = 1; i < m_ - 1; i++) {
    /* Reaction for next full time step */
    sol_next = cap.cwiseQuotient(ones + (cap.cwiseQuotient(sol_cur) - ones) *
                                            std::exp(-lambda_ * tau_));
    sol_cur = sol_next;
    /* Diffusion for next full time step */
    sol_next = StrangSplit::diffusionEvolutionOperator(firstcall, sol_cur);
    sol_cur = sol_next;
  }

  /* Last half time step: Reaction */
  sol_next = cap.cwiseQuotient(ones + (cap.cwiseQuotient(sol_cur) - ones) *
                                          std::exp(-lambda_ * tau_ / 2.));
  sol_cur = sol_next;

  sol = sol_cur;
  return sol;
}

} /* namespace FisherKPP */
