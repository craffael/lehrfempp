#ifndef THESIS_DISCRETIZATION_DIRAC_OPERATOR_SOURCE_PROBLEM_H
#define THESIS_DISCRETIZATION_DIRAC_OPERATOR_SOURCE_PROBLEM_H

/**
 * @file dirac_operator_source_problem.h
 * @brief Class to discretise the dirac_operator.
 */

#include <dirac_operator.h>
#include <lf/assemble/assembler.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/mesh_interface.h>
#include <lf/quad/quad.h>
#include <load_vector_provider.h>
#include <mass_matrix_provider.h>
#include <sphere_triag_mesh_builder.h>
#include <whitney_one_curl_curl_matrix_provider.h>
#include <whitney_one_grad_matrix_provider.h>
#include <whitney_one_mass_matrix_provider.h>
#include <whitney_two_mass_matrix_provider.h>

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere {

namespace discretization {

/**
 * @brief Computes the Galerkin LSE for the Dirac Operator source problem
 *
 * The computation must be complex valued since the problem definition involves
 * a complex i
 *
 * @note Only triangular meshes are supported
 *
 */
using SCALAR = std::complex<double>;
class DiracOperatorSourceProblem {
 public:
  /**
   * @brief Constructor
   * creates basic mesh (Octaeder with radius 1.0)
   *
   * set k to 1
   *
   * For all three Hodge Laplacians
   * creates zerovalued functions f
   *
   */
  DiracOperatorSourceProblem()
      : coo_matrix_(lf::assemble::COOMatrix<SCALAR>(1, 1)) {
    // create mesh factory for basic mesh
    std::unique_ptr<lf::mesh::MeshFactory> factory =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);

    std::shared_ptr<projects::hldo_sphere::mesh::SphereTriagMeshBuilder>
        sphere = std::make_shared<
            projects::hldo_sphere::mesh::SphereTriagMeshBuilder>(
            std::move(factory));
    sphere->setRefinementLevel(0);
    sphere->setRadius(1);

    mesh_p_ = sphere->Build();

    k_ = 1.;

    phi_ = Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>(1);

    mu_ = Eigen::Matrix<SCALAR, Eigen::Dynamic, 1>(1);

    // create basic function
    auto f_0 = [](Eigen::Matrix<double, 3, 1> x) -> double { return 0; };
    auto f_1 =
        [](Eigen::Matrix<double, 3, 1> x) -> Eigen::Matrix<double, 3, 1> {
      return Eigen::Matrix<double, 3, 1>::Zero();
    };
    auto f_2 = [](Eigen::Matrix<double, 3, 1> x) -> double { return 0; };
    f0_ = f_0;
    f1_ = f_1;
    f2_ = f_2;
  }

  /**
   * @brief Computes the Galerkin LSE for the dirac operator source problem
   *
   * @f[
   *   D \vec{u} + \imath k \vec{u} = \vec{f}
   * @f]
   *
   * The Galerkin matrix will be accessable with `get_galerkin_matrix`
   * The load vector will be accessable with `get_load_vector`
   *
   */
  void Compute() {
    // get Dirac Operator Matrix
    projects::hldo_sphere::discretization::DiracOperator dirac_operator;
    dirac_operator.SetLoadFunctions(f0_, f1_, f2_);
    dirac_operator.SetMesh(mesh_p_);
    dirac_operator.Compute();

    lf::assemble::COOMatrix<SCALAR> coo_mat =
        dirac_operator.GetGalerkinMatrix();

    // get righthandside vector
    Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> phi =
        dirac_operator.GetLoadVector();
    phi_ = phi;

    //**********************
    // create mass matrices
    //**********************

    // Zero form mass matrix
    projects::hldo_sphere::assemble::MassMatrixProvider
        mass_matrix_provider_zero;

    const lf::assemble::DofHandler &dof_handler_zero =
        lf::assemble::UniformFEDofHandler(mesh_p_,
                                          {{lf::base::RefEl::kPoint(), 1}});
    const lf::assemble::size_type n_dofs_zero(dof_handler_zero.NumDofs());
    lf::assemble::COOMatrix<SCALAR> coo_mass_mat_zero(n_dofs_zero, n_dofs_zero);
    coo_mass_mat_zero.setZero();
    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<SCALAR>>(
        0, dof_handler_zero, dof_handler_zero, mass_matrix_provider_zero,
        coo_mass_mat_zero);

    for (Eigen::Triplet<SCALAR> triplet : coo_mass_mat_zero.triplets()) {
      int col = triplet.col();
      int row = triplet.row();
      SCALAR val = std::complex<double>(0., 1.) * k_ * triplet.value();
      coo_mat.AddToEntry(row, col, val);
    };

    // one form mass matrix
    projects::hldo_sphere::assemble::WhitneyOneMassMatrixProvider
        mass_matrix_provider_one;
    const lf::assemble::DofHandler &dof_handler_one =
        lf::assemble::UniformFEDofHandler(mesh_p_,
                                          {{lf::base::RefEl::kSegment(), 1}});
    const lf::assemble::size_type n_dofs_one(dof_handler_one.NumDofs());
    lf::assemble::COOMatrix<SCALAR> coo_mass_mat_one(n_dofs_one, n_dofs_one);
    coo_mass_mat_one.setZero();
    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<SCALAR>>(
        0, dof_handler_one, dof_handler_one, mass_matrix_provider_one,
        coo_mass_mat_one);

    for (Eigen::Triplet<SCALAR> triplet : coo_mass_mat_one.triplets()) {
      int col = triplet.col() + n_dofs_zero;
      int row = triplet.row() + n_dofs_zero;
      SCALAR val = std::complex<double>(0., 1.) * k_ * triplet.value();
      coo_mat.AddToEntry(row, col, val);
    }

    // two form mass matrix
    projects::hldo_sphere::assemble::WhitneyTwoMassMatrixProvider
        mass_matrix_provider_two;
    const lf::assemble::DofHandler &dof_handler_two =
        lf::assemble::UniformFEDofHandler(mesh_p_,
                                          {{lf::base::RefEl::kTria(), 1}});
    const lf::assemble::size_type n_dofs_two(dof_handler_two.NumDofs());
    lf::assemble::COOMatrix<SCALAR> coo_mass_mat_two(n_dofs_two, n_dofs_two);
    coo_mass_mat_two.setZero();
    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<SCALAR>>(
        0, dof_handler_two, dof_handler_two, mass_matrix_provider_two,
        coo_mass_mat_two);

    for (Eigen::Triplet<SCALAR> triplet : coo_mass_mat_two.triplets()) {
      // n_dofs_one contains the number of edges and hence the
      // dimension of A_{11}
      int col = triplet.col() + n_dofs_zero + n_dofs_one;
      int row = triplet.row() + n_dofs_zero + n_dofs_one;
      SCALAR val = std::complex<double>(0., 1.) * k_ * triplet.value();
      coo_mat.AddToEntry(row, col, val);
    }

    coo_matrix_ = coo_mat;
  }

  /**
   *
   * @brief solves the linear systems build in Compute
   *
   * Uses the Matrices stored in `coo_matrix_` and
   * the righthandside vectors stored in `phi_`
   * to compute the basis expansion
   * coefficiants mu_ approximating
   *
   * @note the method `Compute` has to be called prior to
   * calling `Solve`
   *
   * @f[
   *   D \vec{u} + \imath k \vec{u} = \vec{f}
   * @f]
   *
   */
  void Solve() {
    Eigen::SparseLU<Eigen::SparseMatrix<SCALAR>> solver;
    Eigen::SparseMatrix<SCALAR> sparse_mat = coo_matrix_.makeSparse();
    solver.compute(sparse_mat);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("Could not decompose the matrix");
    }

    mu_ = solver.solve(phi_);
  }

  /**
   * @brief Sets the mesh
   * @param mesh_p pointer to the mesh
   *
   * requries all cells in the mesh are triangles
   * requries mesh global dimension to be 3
   */
  void SetMesh(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
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

  /**
   * @brief Sets the load functions
   * @param f0 load function in L^2
   * @param f1 load functions in L^2_t tangential vector field on the sphere
   * @param f2 load functions in L^2
   */
  void SetLoadFunctions(
      std::function<double(const Eigen::Matrix<double, 3, 1> &)> f0,
      std::function<
          Eigen::Matrix<double, 3, 1>(const Eigen::Matrix<double, 3, 1> &)>
          f1,
      std::function<double(const Eigen::Matrix<double, 3, 1> &)> f2) {
    f0_ = f0;
    f1_ = f1;
    f2_ = f2;
  }

  /**
   * @brief Sets k for the souce problems
   *
   * @param k real valued parameter scaling the mass matrix
   */
  void SetK(SCALAR k) {
    // check if real valued
    LF_ASSERT_MSG(std::abs(k.imag()) < 1e-13,
                  "K must be real but has imaginary part" << k.imag());
    k_ = k;
  }

  /**
   * @brief returns the Loadvector
   *
   * This is the righthandside of LSE
   *
   * @note The loadvector must be computed with `Compute` before calling
   * this function
   *
   */
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> GetLoadVector() { return phi_; }

  /**
   * @brief returns the Galerkin Matrix
   *
   * This is the Matrix of the LSE
   *
   * @note The Galerkin matrix must be computed with `Compute` before
   * calling this funciton
   *
   */
  lf::assemble::COOMatrix<SCALAR> GetGalerkinMatrix() { return coo_matrix_; }

  /**
   * @brief retunrs the basis expansion coefficiants for
   * the solution of the source problem
   *
   * The basis expansion coefficiants for the first
   * number of vertices entries are with respect to the
   * barycentric basis functions on the mesh.
   *
   * The basis expansion coefficiants for the next
   * number of edge entries are with respect to the
   * whitney one form basis functions.
   *
   * The last number of cells basis expansion coefficiants
   * are with respect to the cellwise constant basis functions
   *
   * @returns basis expansion coefficiants
   *
   * @note the methods `Compute` and `Solve` must be
   * called before the result is available
   *
   */
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> GetMu() { return mu_; }

 private:
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> f0_;
  std::function<Eigen::Matrix<double, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      f1_;
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> f2_;
  lf::assemble::COOMatrix<SCALAR> coo_matrix_;
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> phi_;
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> mu_;
  SCALAR k_;
};

}  // namespace discretization

}  // namespace projects::hldo_sphere

#endif  // THESIS_DISCRETIZATION_DIRAC_OPERATOR_SOURCE_PROBLEM_H
