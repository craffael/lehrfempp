#ifndef THESIS_DISCRETIZATION_HODGE_LAPLACIANS_SOURCE_PROBLEMS_H
#define THESIS_DISCRETIZATION_HODGE_LAPLACIANS_SOURCE_PROBLEMS_H

/**
 * @file hodge_laplacians_source_problems.h
 * @brief Class to discretise the hodge_laplacian source problems
 */

#include <mass_matrix_provider.h>
#include <whitney_one_hodge_laplacian.h>
#include <whitney_one_mass_matrix_provider.h>
#include <whitney_two_hodge_laplacian.h>
#include <whitney_two_mass_matrix_provider.h>
#include <whitney_zero_hodge_laplacian.h>

#include <Eigen/Dense>
#include <cmath>

namespace projects::hldo_sphere {

namespace discretization {

/**
 * @brief Creates and solves the Discretised Hodge Laplacian source problems
 *
 * @note Only triangular meshes are supported
 *
 */
class HodgeLaplaciansSourceProblems {
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
  HodgeLaplaciansSourceProblems() {
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

    // initialize matrix vector
    coo_matrix_ = std::vector<lf::assemble::COOMatrix<double>>{
        lf::assemble::COOMatrix<double>(1, 1),
        lf::assemble::COOMatrix<double>(1, 1),
        lf::assemble::COOMatrix<double>(1, 1)};

    phi_ = std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>>(3);

    mu_ = std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>>(3);

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
   * @brief Computes the Galerkin LSE for all three source problems
   *
   * @f[
   * - \Delta_l u + k^2 u = f_l \qquad l = 0, 1, 2
   * @f]
   *
   * The Galerkin matrix will be accessable with `get_galerkin_matrix(int
   * index)` The load vector will be accessable with `get_load_vector(int
   * index)`
   *
   */
  void Compute() {
    //*****************
    // Whitney 0 form
    //*****************

    // get the Hodge laplacian for the zero form and take the negative
    projects::hldo_sphere::discretization::WhitneyZeroHodgeLaplace zero_builder;
    zero_builder.SetMesh(mesh_p_);
    zero_builder.SetLoadFunction(f0_);
    zero_builder.Compute();
    lf::assemble::COOMatrix<double> coo_mat_zero_pos =
        zero_builder.GetGalerkinMatrix();
    lf::assemble::COOMatrix<double> coo_mat_zero(coo_mat_zero_pos.rows(),
                                                 coo_mat_zero_pos.cols());
    coo_mat_zero.setZero();
    for (Eigen::Triplet<double> triplet : coo_mat_zero_pos.triplets()) {
      int col = triplet.col();
      int row = triplet.row();
      double val = -triplet.value();
      coo_mat_zero.AddToEntry(row, col, val);
    }

    // create mass Matrix

    projects::hldo_sphere::assemble::MassMatrixProvider
        mass_matrix_provider_zero;

    const lf::assemble::DofHandler &dof_handler_zero =
        lf::assemble::UniformFEDofHandler(mesh_p_,
                                          {{lf::base::RefEl::kPoint(), 1}});
    const lf::assemble::size_type n_dofs_zero(dof_handler_zero.NumDofs());
    lf::assemble::COOMatrix<double> coo_mass_mat_zero(n_dofs_zero, n_dofs_zero);
    coo_mass_mat_zero.setZero();
    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
        0, dof_handler_zero, dof_handler_zero, mass_matrix_provider_zero,
        coo_mass_mat_zero);

    // merge whitney zero for with k^2 * mass Matrix
    for (Eigen::Triplet<double> triplet : coo_mass_mat_zero.triplets()) {
      int col = triplet.col();
      int row = triplet.row();
      double val = k_ * k_ * triplet.value();
      coo_mat_zero.AddToEntry(row, col, val);
    };

    coo_matrix_[0] = coo_mat_zero;

    // get righthandside
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi_zero =
        zero_builder.GetLoadVector();
    phi_[0] = phi_zero;

    //*****************
    // Whitney 1 form
    //*****************

    // get the Hodge laplacian for the one form and take the negative
    projects::hldo_sphere::discretization::WhitneyOneHodgeLaplace one_builder;
    one_builder.SetMesh(mesh_p_);
    one_builder.SetLoadFunction(f1_);
    one_builder.Compute();
    lf::assemble::COOMatrix<double> coo_mat_one_pos =
        one_builder.GetGalerkinMatrix();
    lf::assemble::COOMatrix<double> coo_mat_one(coo_mat_one_pos.rows(),
                                                coo_mat_one_pos.cols());
    coo_mat_one.setZero();
    for (Eigen::Triplet<double> triplet : coo_mat_one_pos.triplets()) {
      int col = triplet.col();
      int row = triplet.row();
      double val = triplet.value();
      coo_mat_one.AddToEntry(row, col, val);
    }

    // create mass Matrix
    projects::hldo_sphere::assemble::WhitneyOneMassMatrixProvider
        mass_matrix_provider_one;
    const lf::assemble::DofHandler &dof_handler_one =
        lf::assemble::UniformFEDofHandler(mesh_p_,
                                          {{lf::base::RefEl::kSegment(), 1}});
    const lf::assemble::size_type n_dofs_one(dof_handler_one.NumDofs());
    lf::assemble::COOMatrix<double> coo_mass_mat_one(n_dofs_one, n_dofs_one);
    coo_mass_mat_one.setZero();
    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
        0, dof_handler_one, dof_handler_one, mass_matrix_provider_one,
        coo_mass_mat_one);

    // merge whitney one for with k^2 * mass Matrix
    for (Eigen::Triplet<double> triplet : coo_mass_mat_one.triplets()) {
      int col = triplet.col();
      int row = triplet.row();
      double val = k_ * k_ * triplet.value();
      coo_mat_one.AddToEntry(row, col, val);
    }
    coo_matrix_[1] = coo_mat_one;

    // get righthandside
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi_one =
        one_builder.GetLoadVector();
    phi_[1] = phi_one;

    //*****************
    // Whitney 2 form
    //*****************

    // get the Hodge laplacian for the one form and take the negative
    projects::hldo_sphere::discretization::WhitneyTwoHodgeLaplace two_builder;
    two_builder.SetMesh(mesh_p_);
    two_builder.SetLoadFunction(f2_);
    two_builder.Compute();
    lf::assemble::COOMatrix<double> coo_mat_two_pos =
        two_builder.GetGalerkinMatrix();
    lf::assemble::COOMatrix<double> coo_mat_two(coo_mat_two_pos.rows(),
                                                coo_mat_two_pos.cols());
    coo_mat_two.setZero();
    for (Eigen::Triplet<double> triplet : coo_mat_two_pos.triplets()) {
      int col = triplet.col();
      int row = triplet.row();
      double val = triplet.value();
      coo_mat_two.AddToEntry(row, col, val);
    }

    // create mass Matrix
    projects::hldo_sphere::assemble::WhitneyTwoMassMatrixProvider
        mass_matrix_provider_two;
    const lf::assemble::DofHandler &dof_handler_two =
        lf::assemble::UniformFEDofHandler(mesh_p_,
                                          {{lf::base::RefEl::kTria(), 1}});
    const lf::assemble::size_type n_dofs_two(dof_handler_two.NumDofs());
    lf::assemble::COOMatrix<double> coo_mass_mat_two(n_dofs_two, n_dofs_two);
    coo_mass_mat_two.setZero();
    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
        0, dof_handler_two, dof_handler_two, mass_matrix_provider_two,
        coo_mass_mat_two);

    // merge whitney two for with k^2 * mass Matrix
    for (Eigen::Triplet<double> triplet : coo_mass_mat_two.triplets()) {
      // n_dofs_one contains the number of edges and hence the
      // dimension of A_{11}
      int col = triplet.col() + n_dofs_one;
      int row = triplet.row() + n_dofs_one;
      double val = k_ * k_ * triplet.value();
      coo_mat_two.AddToEntry(row, col, val);
    }

    coo_matrix_[2] = coo_mat_two;

    // get righthandside
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi_two =
        two_builder.GetLoadVector();
    phi_[2] = phi_two;
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
   * - \Delta_l u + k^2 u = f_l \qquad l = 0, 1, 2
   * @f]
   *
   */
  void Solve() {
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    for (int l = 0; l < 3; l++) {
      Eigen::SparseMatrix<double> sparse_mat = coo_matrix_[l].makeSparse();
      solver.compute(sparse_mat);
      if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Could not decompose the matrix");
      }

      mu_[l] = solver.solve(phi_[l]);
    }
  }

  /**
   * @brief Sets the mesh and creates dof_handler
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
   * @param f1 load functions in L^2_t vector valued
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
  void SetK(double k) { k_ = k; }

  /**
   * @brief returns the Loadvector
   *
   * This is the righthandside of LSE with index `index`
   *
   * @note The loadvector must be computed with `Compute` before calling
   * this function
   *
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> GetLoadVector(int index) {
    return phi_[index];
  }

  /**
   * @brief returns the Galerkin Matrix
   *
   * This is the Matrix of the LSE with index `index`
   *
   * @note The Galerkin matrix must be computed with `Compute` before
   * calling this funciton
   *
   */
  lf::assemble::COOMatrix<double> GetGalerkinMatrix(int index) {
    return coo_matrix_[index];
  }

  /**
   * @brief retunrs the basis expansion coefficiants for
   * the solution of the zero form
   *
   * The basis expansion coefficiants are with respect to the
   * barycentric basis functions on the mesh.
   *
   * @returns basis expansion coefficiants of zero form
   *
   * @note the methods `Compute` and `Solve` must be
   * called before the result is available
   *
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> GetMuZero() { return mu_[0]; }

  /**
   * @brief retunrs the basis expansion coefficiants for
   * the solution of the one form
   *
   * The first vector of the return value are the
   * coefficiants with respect to whitney one form
   * basis function. They represent an approximation of
   * @f[ u @f]
   *
   * The second vector of the return value are the
   * coefficiants with respect to to the barycentric
   * basis functions. They represent the approximation
   * for
   * @f[ p := div_{\Gamma}(u) @f]
   * on the mesh.
   *
   * @returns basis expansion coefficiants of one form
   *
   * @note the methods `Compute` and `Solve` must be
   * called before the result is available
   *
   */
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
             Eigen::Matrix<double, Eigen::Dynamic, 1>>
  GetMuOne() {
    lf::base::size_type numPoints = mesh_p_->NumEntities(2);
    lf::base::size_type numEdges = mesh_p_->NumEntities(1);
    Eigen::Matrix<double, Eigen::Dynamic, 1> u = mu_[1].head(numEdges);
    Eigen::Matrix<double, Eigen::Dynamic, 1> p = mu_[1].tail(numPoints);
    return std::make_tuple(u, p);
  }

  /**
   * @brief retunrs the basis expansion coefficiants for
   * the solution of the two form
   *
   * The first vector of the return value are the
   * coefficiants with respect to rotated
   * (90 degree colckwise) whitney one form
   * basis function. They represent an approximation of
   * @f[ j = grad_{\Gamma}(u) @f]
   *
   * The second vector of the return value are the
   * coefficiants with respect to to the piecewise constant
   * basis functions. They represent the approximation
   * for
   * @f[ u @f]
   * on the mesh.
   *
   * @returns basis expansion coefficiants of two form
   *
   * @note the methods `Compute` and `Solve` must be
   * called before the result is available
   *
   */
  std::tuple<Eigen::Matrix<double, Eigen::Dynamic, 1>,
             Eigen::Matrix<double, Eigen::Dynamic, 1>>
  GetMuTwo() {
    lf::base::size_type numCells = mesh_p_->NumEntities(0);
    lf::base::size_type numEdges = mesh_p_->NumEntities(1);
    Eigen::Matrix<double, Eigen::Dynamic, 1> j = mu_[1].head(numEdges);
    Eigen::Matrix<double, Eigen::Dynamic, 1> u = mu_[1].tail(numCells);
    return std::make_tuple(j, u);
  }

 private:
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> f0_;
  std::function<Eigen::Matrix<double, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      f1_;
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> f2_;
  std::vector<lf::assemble::COOMatrix<double>> coo_matrix_;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> phi_;
  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> mu_;
  double k_ = 1;
};

}  // namespace discretization

}  // namespace projects::hldo_sphere

#endif  // THESIS_DISCRETIZATION_HODGE_LAPLACIANS_SOURCE_PROBLEMS_H
