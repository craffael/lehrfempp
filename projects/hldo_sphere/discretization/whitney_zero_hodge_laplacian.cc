#include "whitney_zero_hodge_laplacian.h"

#include <lf/uscalfe/uscalfe.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::discretization {

void WhitneyZeroHodgeLaplace::Compute() {
  // create element matrix provider
  lf::uscalfe::LinearFELaplaceElementMatrix matrix_provider;

  // create fe_space for the ScalarLoadElementVectorProvider
  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p_);

  // create COO Matrix
  const lf::assemble::size_type n_dofs(dof_handler_p_->NumDofs());
  lf::assemble::COOMatrix<double> coo_mat(n_dofs, n_dofs);
  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, *dof_handler_p_, *dof_handler_p_, matrix_provider, coo_mat);

  // create element vector provider
  lf::mesh::utils::MeshFunctionGlobal mf_f(f_);
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      vector_provider{fe_space_p, mf_f};

  // create load vector
  Eigen::VectorXd phi(dof_handler_p_->NumDofs());

  // assemble the global vector over entities with codim=0:
  AssembleVectorLocally(0, *dof_handler_p_, vector_provider, phi);

  // set the class attributes
  phi_ = phi;
  coo_matrix_p_ = std::shared_ptr<lf::assemble::COOMatrix<double>>{&coo_mat};
}

}  // namespace projects::hldo_sphere::discretization
