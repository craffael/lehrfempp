#include "whitney_zero_hodge_laplacian.h"

#include <laplace_element_matrix_provider.h>
#include <lf/assemble/assembler.h>
#include <load_element_vector_provider.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::discretization {

void WhitneyZeroHodgeLaplace::WhitneyZeroHodgeLaplace::Compute() {
  // create element matrix provider
  projects::hldo_sphere::assemble::LaplaceElementMatrixProvider matrix_provider;

  const lf::assemble::DofHandler &dof_handler =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kPoint(), 1}});

  // create COO Matrix
  const lf::assemble::size_type n_dofs(dof_handler.NumDofs());

  lf::assemble::COOMatrix<double> coo_mat(n_dofs, n_dofs);
  coo_mat.setZero();

  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dof_handler, dof_handler, matrix_provider, coo_mat);

  // create element vector provider
  projects::hldo_sphere::assemble::LoadElementVectorProvider vector_provider{
      f_};

  // create load vector
  Eigen::VectorXd phi(dof_handler.NumDofs());
  phi.setZero();

  // assemble the global vector over entities with codim=0:
  AssembleVectorLocally(0, dof_handler, vector_provider, phi);

  // set the class attributes
  phi_ = phi;
  coo_matrix_ = coo_mat;
}

}  // namespace projects::hldo_sphere::discretization
