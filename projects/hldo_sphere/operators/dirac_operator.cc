#include "dirac_operator.h"

namespace projects::hldo_sphere {

namespace operators {

void DiracOperator::Compute() {
  // create element matrix provider for the matrices A12, A21
  projects::hldo_sphere::assemble::WhitneyOneGradMatrixProvider
      matrix_grad_provider;

  // create matrix provider for  -A23, -A32
  projects::hldo_sphere::assemble::RotWhitneyOneDivMatrixProvider
      matrix_div_provider;

  const lf::assemble::DofHandler &dof_handler_vert =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kPoint(), 1}});
  const lf::assemble::DofHandler &dof_handler_edge =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kSegment(), 1}});
  const lf::assemble::DofHandler &dof_handler_cell =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kTria(), 1}});

  // create COO Matrix
  const lf::assemble::size_type n_dofs_vert(dof_handler_vert.NumDofs());
  const lf::assemble::size_type n_dofs_edge(dof_handler_edge.NumDofs());
  const lf::assemble::size_type n_dofs_cell(dof_handler_cell.NumDofs());

  lf::assemble::COOMatrix<double> coo_A_21(n_dofs_edge, n_dofs_vert);
  coo_A_21.setZero();

  // the m stands for minus
  // note that the RotWOneFormDivElementMatrixProvider returns the
  // negative of the elements we want
  lf::assemble::COOMatrix<double> coo_A_23_m(n_dofs_edge, n_dofs_cell);
  coo_A_23_m.setZero();

  // create s matrix with n_dofs_edge rows and n_dofs_vert columns
  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dof_handler_vert, dof_handler_edge, matrix_grad_provider, coo_A_21);

  // create matrix with n_dof_edge_rows and n_dof_cell columns
  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
      0, dof_handler_cell, dof_handler_edge, matrix_div_provider, coo_A_23_m);

  // create full matrix
  lf::assemble::COOMatrix<complex> full_matrix(
      n_dofs_vert + n_dofs_edge + n_dofs_cell,
      n_dofs_vert + n_dofs_edge + n_dofs_cell);

  full_matrix.setZero();

  // iterate over all triplets of the previously computed matrice and add up
  // entries
  const std::vector<Eigen::Triplet<double>> triplets_A_21 = coo_A_21.triplets();
  const std::vector<Eigen::Triplet<double>> triplets_A_23_m =
      coo_A_23_m.triplets();

  // Add A_12 and A_21
  for (Eigen::Triplet<double> triplet : triplets_A_21) {
    int col = triplet.col();
    int row = triplet.row();
    complex val = complex(triplet.value(), 0.);
    // A_21
    full_matrix.AddToEntry(row + n_dofs_vert, col, val);
    // A_12
    full_matrix.AddToEntry(col, row + n_dofs_vert, val);
  }

  // Add A_23 and A_32
  for (Eigen::Triplet<double> triplet : triplets_A_23_m) {
    int col = triplet.col();
    int row = triplet.row();
    complex val = complex(triplet.value(), 0.);

    // Add A_23
    full_matrix.AddToEntry(row + n_dofs_vert, col + n_dofs_vert + n_dofs_edge,
                           val);

    // Add A_32
    full_matrix.AddToEntry(col + n_dofs_vert + n_dofs_edge, row + n_dofs_vert,
                           val);
  }

  coo_matrix_ = full_matrix;

  // create element vector provider
  projects::hldo_sphere::assemble::LoadVectorProvider<complex>
      vector_provider_0(f0_);

  projects::hldo_sphere::assemble::WhitneyOneVectorProvider<complex>
      vector_provider_1(f1_);

  projects::hldo_sphere::assemble::WhitneyTwoVectorProvider<complex>
      vector_provider_2(f2_);

  // create load vector
  Eigen::Matrix<complex, Eigen::Dynamic, 1> phi0(n_dofs_vert);
  phi0.setZero();
  Eigen::Matrix<complex, Eigen::Dynamic, 1> phi1(n_dofs_edge);
  phi1.setZero();
  Eigen::Matrix<complex, Eigen::Dynamic, 1> phi2(n_dofs_cell);
  phi2.setZero();

  // assemble the global vector over entities with codim=0:
  AssembleVectorLocally(0, dof_handler_vert, vector_provider_0, phi0);
  AssembleVectorLocally(0, dof_handler_edge, vector_provider_1, phi1);
  AssembleVectorLocally(0, dof_handler_cell, vector_provider_2, phi2);

  // create full vector
  Eigen::Matrix<complex, Eigen::Dynamic, 1> full_vec(n_dofs_vert + n_dofs_edge +
                                                     n_dofs_cell);
  full_vec.setZero();
  full_vec.head(n_dofs_vert) = phi0;
  full_vec.segment(n_dofs_vert, n_dofs_edge) = phi1;
  full_vec.tail(n_dofs_cell) = phi2;

  // set the class attributes
  phi_ = full_vec;
}

/**
 * @brief Sets the mesh and creates dof_handler
 * @param mesh_p pointer to the mesh
 *
 * requries all cells in the mesh are triangles
 * requries mesh global dimension to be 3
 */
void DiracOperator::SetMesh(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  // check if cells are triagles
  for (const lf::mesh::Entity *tria : mesh_p->Entities(0)) {
    LF_ASSERT_MSG(
        tria->RefEl() == lf::base::RefEl::kTria(),
        "Mesh must be Triangular, unsupported cell " << tria->RefEl());
  }

  // check if dimension of the mesh is 3
  LF_ASSERT_MSG(mesh_p->DimWorld() == 3,
                "World Dimension must be 3 but is" << mesh_p->DimWorld());

  // set mesh
  mesh_p_ = mesh_p;
}

}  // namespace operators
}  // namespace projects::hldo_sphere
