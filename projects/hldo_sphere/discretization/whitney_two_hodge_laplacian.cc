#include "whitney_two_hodge_laplacian.h"

#include <curl_element_vector_provider.h>
#include <lf/assemble/assembler.h>
#include <lf/quad/quad.h>
#include <load_element_vector_provider.h>
#include <mass_element_matrix_provider.h>
#include <rot_w_one_form_div_element_matrix_provider.h>
#include <rot_w_one_form_dot_element_matrix_provider.h>
#include <whitney_two_element_vector_provider.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::discretization {

void WhitneyTwoHodgeLaplace::WhitneyTwoHodgeLaplace::Compute() {
  // create element matrix provider
  projects::hldo_sphere::assemble::RotWOneFormDotElementMatrixProvider
      matrix_dot_provider;
  projects::hldo_sphere::assemble::RotWOneFormDivElementMatrixProvider
      matrix_curl_provider;

  const lf::assemble::DofHandler &dof_handler_div =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kSegment(), 1}});
  const lf::assemble::DofHandler &dof_handler_const =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kTria(), 1}});

  // create COO Matrix
  const lf::assemble::size_type n_dofs_div(dof_handler_div.NumDofs());
  const lf::assemble::size_type n_dofs_const(dof_handler_const.NumDofs());

  lf::assemble::COOMatrix<double> coo_A_11(n_dofs_div, n_dofs_div);
  coo_A_11.setZero();
  lf::assemble::COOMatrix<double> coo_A_12(n_dofs_div, n_dofs_const);
  coo_A_12.setZero();

  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dof_handler_div, dof_handler_div, matrix_dot_provider, coo_A_11);

  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dof_handler_const, dof_handler_div, matrix_curl_provider, coo_A_12);

  // create full matrix
  lf::assemble::COOMatrix<double> full_matrix(n_dofs_div + n_dofs_const,
                                              n_dofs_div + n_dofs_const);
  full_matrix.setZero();

  // iterate over all triplets of the previously computed matrice and add up
  // entries
  const std::vector<Eigen::Triplet<double>> triplets_A_11 = coo_A_11.triplets();
  const std::vector<Eigen::Triplet<double>> triplets_A_12 = coo_A_12.triplets();

  // Add A_11
  for (Eigen::Triplet<double> triplet : triplets_A_11) {
    int col = triplet.col();
    int row = triplet.row();
    double val = triplet.value();
    full_matrix.AddToEntry(row, col, val);
  }

  // Add A_12 and A_21
  for (Eigen::Triplet<double> triplet : triplets_A_12) {
    int col = triplet.col() + n_dofs_div;
    int row = triplet.row();
    double val = triplet.value();

    // Add A_12
    full_matrix.AddToEntry(row, col, val);

    // Add A_21
    full_matrix.AddToEntry(col, row, val);
  }

  coo_matrix_ = full_matrix;

  // define quad rule with sufficiantly high degree since the
  // baricentric coordinate functions and whitney one form basis functions have
  // degree 1
  lf::quad::QuadRule quadrule{lf::quad::make_TriaQR_EdgeMidpointRule()};

  // create element vector provider
  projects::hldo_sphere::assemble::WhitneyTwoElementVectorProvider
      vector_provider(f_, quadrule);

  // create load vector
  Eigen::VectorXd phi(n_dofs_const);
  phi.setZero();

  // assemble the global vector over entities with codim=0:
  lf::assemble::AssembleVectorLocally(0, dof_handler_const, vector_provider,
                                      phi);

  // create full vector
  Eigen::VectorXd full_vec(n_dofs_div + n_dofs_const);
  full_vec.setZero();
  full_vec.tail(n_dofs_const) = -phi;

  // set the class attributes
  phi_ = full_vec;
}

}  // namespace projects::hldo_sphere::discretization
