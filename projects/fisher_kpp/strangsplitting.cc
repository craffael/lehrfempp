/** @file
 *	@brief Bachelor Thesis Fisher/KPP
 *	@author Amélie Justine Loher
 *	@date 01.04.20
 *  @copyright ETH Zurich
 */

#include "strangsplitting.h"

namespace FisherKPP {

/* Member Function StrangSplit
 * Computes the Evolution Operator for the linear parabolic diffusion term
 */
Eigen::VectorXd StrangSplit::diffusionEvolutionOperator(
    bool firstcall, const Eigen::VectorXd &mu) {
  Eigen::VectorXd evol_op;

  /* Precomputation */
  /*  solver.compute(M_ + tau * kappa_ * A_);
    LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
   */
  const Eigen::VectorXd rhs = -A_ * mu;

  if (firstcall) {
    /* First stage SDIRK-2 */
    const Eigen::VectorXd k1 = solver1.solve(rhs);
    LF_VERIFY_MSG(solver1.info() == Eigen::Success, "LU decomposition failed");
    /* Second stage SDIRK-2 */
    const Eigen::VectorXd k2 =
        solver1.solve(rhs - 0.5 * tau_ * (1 - kappa_) * A_ * k1);
    LF_VERIFY_MSG(solver1.info() == Eigen::Success, "LU decomposition failed");

    evol_op = mu + 0.5 * tau_ * (1 - kappa_) * k1 + 0.5 * tau_ * kappa_ * k2;

  } else {
    /* First stage SDIRK-2 */
    const Eigen::VectorXd k1 = solver2.solve(rhs);
    LF_VERIFY_MSG(solver2.info() == Eigen::Success, "LU decomposition failed");
    /* Second stage SDIRK-2 */
    const Eigen::VectorXd k2 =
        solver2.solve(rhs - tau_ * (1 - kappa_) * A_ * k1);
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

  const Eigen::VectorXd ones = Eigen::VectorXd::Ones(N_dofs);

  /* Time Steps */
  // double tau = T_ / m_;

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
