#include "dirac_operator_source_problem.h"

namespace projects::hldo_sphere {

namespace operators {

DiracOperatorSourceProblem::DiracOperatorSourceProblem()
    : coo_matrix_(lf::assemble::COOMatrix<complex>(1, 1)) {
  // create mesh factory for basic mesh
  std::unique_ptr<lf::mesh::MeshFactory> factory =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

  std::shared_ptr<projects::hldo_sphere::mesh::SphereTriagMeshBuilder> sphere =
      std::make_shared<projects::hldo_sphere::mesh::SphereTriagMeshBuilder>(
          std::move(factory));
  sphere->setRefinementLevel(0);
  sphere->setRadius(1);

  mesh_p_ = sphere->Build();

  k_ = 1.;

  phi_ = Eigen::Matrix<complex, Eigen::Dynamic, 1>(1);

  mu_ = Eigen::Matrix<complex, Eigen::Dynamic, 1>(1);

  // create basic function
  auto f_0 = [](Eigen::Matrix<double, 3, 1> x) -> complex { return 0; };
  auto f_1 = [](Eigen::Matrix<double, 3, 1> x) -> Eigen::Matrix<complex, 3, 1> {
    return Eigen::Matrix<complex, 3, 1>::Zero();
  };
  auto f_2 = [](Eigen::Matrix<double, 3, 1> x) -> complex { return 0; };
  f0_ = f_0;
  f1_ = f_1;
  f2_ = f_2;
}

void DiracOperatorSourceProblem::Compute() {
  // get Dirac Operator Matrix
  projects::hldo_sphere::operators::DiracOperator dirac_operator;
  dirac_operator.SetLoadFunctions(f0_, f1_, f2_);
  dirac_operator.SetMesh(mesh_p_);
  dirac_operator.Compute();

  lf::assemble::COOMatrix<complex> coo_mat = dirac_operator.GetGalerkinMatrix();

  // get righthandside vector
  Eigen::Matrix<complex, Eigen::Dynamic, 1> phi =
      dirac_operator.GetLoadVector();
  phi_ = phi;

  //**********************
  // create mass matrices
  //**********************

  // Zero form mass matrix
  projects::hldo_sphere::assemble::MassMatrixProvider mass_matrix_provider_zero;

  const lf::assemble::DofHandler &dof_handler_zero =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kPoint(), 1}});
  const lf::assemble::size_type n_dofs_zero(dof_handler_zero.NumDofs());

  lf::assemble::COOMatrix<complex> coo_mass_mat_zero(n_dofs_zero, n_dofs_zero);
  coo_mass_mat_zero.setZero();

  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<complex>>(
      0, dof_handler_zero, dof_handler_zero, mass_matrix_provider_zero,
      coo_mass_mat_zero);

  for (Eigen::Triplet<complex> triplet : coo_mass_mat_zero.triplets()) {
    int col = triplet.col();
    int row = triplet.row();
    complex val = std::complex<double>(0., 1.) * k_ * triplet.value();
    coo_mat.AddToEntry(row, col, val);
  };

  // one form mass matrix
  projects::hldo_sphere::assemble::WhitneyOneMassMatrixProvider
      mass_matrix_provider_one;
  const lf::assemble::DofHandler &dof_handler_one =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kSegment(), 1}});
  const lf::assemble::size_type n_dofs_one(dof_handler_one.NumDofs());
  lf::assemble::COOMatrix<complex> coo_mass_mat_one(n_dofs_one, n_dofs_one);
  coo_mass_mat_one.setZero();
  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<complex>>(
      0, dof_handler_one, dof_handler_one, mass_matrix_provider_one,
      coo_mass_mat_one);

  for (Eigen::Triplet<complex> triplet : coo_mass_mat_one.triplets()) {
    int col = triplet.col() + n_dofs_zero;
    int row = triplet.row() + n_dofs_zero;
    complex val = std::complex<double>(0., 1.) * k_ * triplet.value();
    coo_mat.AddToEntry(row, col, val);
  }

  // two form mass matrix
  projects::hldo_sphere::assemble::WhitneyTwoMassMatrixProvider
      mass_matrix_provider_two;
  const lf::assemble::DofHandler &dof_handler_two =
      lf::assemble::UniformFEDofHandler(mesh_p_,
                                        {{lf::base::RefEl::kTria(), 1}});
  const lf::assemble::size_type n_dofs_two(dof_handler_two.NumDofs());
  lf::assemble::COOMatrix<complex> coo_mass_mat_two(n_dofs_two, n_dofs_two);
  coo_mass_mat_two.setZero();
  lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<complex>>(
      0, dof_handler_two, dof_handler_two, mass_matrix_provider_two,
      coo_mass_mat_two);

  for (Eigen::Triplet<complex> triplet : coo_mass_mat_two.triplets()) {
    // n_dofs_one contains the number of edges and hence the
    // dimension of A_{11}
    int col = triplet.col() + n_dofs_zero + n_dofs_one;
    int row = triplet.row() + n_dofs_zero + n_dofs_one;
    complex val = std::complex<double>(0., 1.) * k_ * triplet.value();
    coo_mat.AddToEntry(row, col, val);
  }

  coo_matrix_ = coo_mat;
}

void DiracOperatorSourceProblem::Solve() {
  Eigen::SparseLU<Eigen::SparseMatrix<complex>> solver;
  Eigen::SparseMatrix<complex> sparse_mat = coo_matrix_.makeSparse();
  sparse_mat.makeCompressed();
  solver.analyzePattern(sparse_mat);
  solver.factorize(sparse_mat);
  if (solver.info() != Eigen::Success) {
    throw std::runtime_error("Could not decompose the matrix");
  }

  mu_ = solver.solve(phi_);
}

void DiracOperatorSourceProblem::SetMesh(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
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

Eigen::Matrix<double, Eigen::Dynamic, 1> DiracOperatorSourceProblem::GetMu(
    int index) {
  LF_ASSERT_MSG(index < 3 && index >= 0,
                "Index must be in {0,1,2}, given " << index);

  Eigen::Vector3d n;
  Eigen::Vector3d s;
  n(0) = mesh_p_->NumEntities(2);
  s(0) = 0;
  n(1) = mesh_p_->NumEntities(1);
  s(1) = n(0);
  n(2) = mesh_p_->NumEntities(0);
  s(2) = s(1) + n(1);
  return mu_.segment(s(index), n(index)).real();
}

}  // namespace operators
}  // namespace projects::hldo_sphere
