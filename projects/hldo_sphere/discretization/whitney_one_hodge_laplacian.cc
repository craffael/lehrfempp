#include "whitney_one_hodge_laplacian.h"

#include <curl_element_vector_provider.h>
#include <lf/assemble/assembler.h>
#include <lf/quad/quad.h>
#include <load_element_vector_provider.h>
#include <mass_element_matrix_provider.h>
#include <whitney_one_form_curl_element_matrix_provider.h>
#include <whitney_one_form_grad_element_matrix_provider.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::discretization {

void WhitneyOneHodgeLaplace::WhitneyOneHodgeLaplace::Compute() {
  // create element matrix provider
  projects::hldo_sphere::assemble::WhitneyOneFormCurlElementMatrixProvider
      matrix_curl_provider;
  projects::hldo_sphere::assemble::WhitneyOneFormGradElementMatrixProvider
      matrix_grad_provider;
  projects::hldo_sphere::assemble::MassElementMatrixProvider
      matrix_mass_provider;

  const lf::assemble::DofHandler &dof_handler_curl =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kSegment(), 1}});
  const lf::assemble::DofHandler &dof_handler_bary =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kPoint(), 1}});

  // create COO Matrix
  const lf::assemble::size_type n_dofs_curl(dof_handler_curl.NumDofs());
  const lf::assemble::size_type n_dofs_bary(dof_handler_bary.NumDofs());

  lf::assemble::COOMatrix<double> coo_A_11(n_dofs_curl, n_dofs_curl);
  coo_A_11.setZero();
  lf::assemble::COOMatrix<double> coo_A_12(n_dofs_curl, n_dofs_bary);
  coo_A_12.setZero();
  lf::assemble::COOMatrix<double> coo_A_22(n_dofs_bary, n_dofs_bary);
  coo_A_22.setZero();

  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dof_handler_curl, dof_handler_curl, matrix_curl_provider, coo_A_11);

  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dof_handler_bary, dof_handler_curl, matrix_grad_provider, coo_A_12);

  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dof_handler_bary, dof_handler_bary, matrix_mass_provider, coo_A_22);

  // create full matrix
  lf::assemble::COOMatrix<double> full_matrix(n_dofs_curl + n_dofs_bary,
                                              n_dofs_curl + n_dofs_bary);
  full_matrix.setZero();

  // iterate over all triplets of the previously computed matrice and add up
  // entries
  const std::vector<Eigen::Triplet<double>> triplets_A_11 = coo_A_11.triplets();
  const std::vector<Eigen::Triplet<double>> triplets_A_12 = coo_A_12.triplets();
  const std::vector<Eigen::Triplet<double>> triplets_A_22 = coo_A_22.triplets();

  // Add A_11
  for (Eigen::Triplet<double> triplet : triplets_A_11) {
    int col = triplet.col();
    int row = triplet.row();
    double val = triplet.value();
    full_matrix.AddToEntry(row, col, val);
  }

  // Add A_12 and A_21
  for (Eigen::Triplet<double> triplet : triplets_A_12) {
    int col = triplet.col() + n_dofs_curl;
    int row = triplet.row();
    double val = triplet.value();

    // Add A_12
    full_matrix.AddToEntry(row, col, -val);

    // Add A_21
    full_matrix.AddToEntry(col, row, val);
  }

  // Add A_22
  for (Eigen::Triplet<double> triplet : triplets_A_22) {
    int col = triplet.col() + n_dofs_curl;
    int row = triplet.row() + n_dofs_curl;
    double val = triplet.value();
    full_matrix.AddToEntry(row, col, val);
  }

  coo_matrix_ = full_matrix;

  // define quad rule with sufficiantly high degree since the
  // baricentric coordinate functions and whitney one form basis functions have
  // degree 1
  lf::quad::QuadRule quadrule{lf::quad::make_TriaQR_EdgeMidpointRule()};

  // create element vector provider
  projects::hldo_sphere::assemble::CurlElementVectorProvider vector_provider(
      f_, quadrule);

  // create load vector
  Eigen::VectorXd phi(n_dofs_curl);
  phi.setZero();

  // assemble the global vector over entities with codim=0:
  AssembleVectorLocally(0, dof_handler_curl, vector_provider, phi);

  // create full vector
  Eigen::VectorXd full_vec(n_dofs_curl + n_dofs_bary);
  full_vec.setZero();
  full_vec.head(n_dofs_curl) = phi;

  // set the class attributes
  phi_ = full_vec;
}

}  // namespace projects::hldo_sphere::discretization
