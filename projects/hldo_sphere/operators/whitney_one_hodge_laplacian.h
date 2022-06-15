#ifndef THESIS_DISCRETIZATION_WHITNEY_ONE_HODGE_LAPLACE_H
#define THESIS_DISCRETIZATION_WHITNEY_ONE_HODGE_LAPLACE_H

/**
 * @file whitney_one_hodge_laplace.h
 */
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
#include <whitney_one_vector_provider.h>

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere {

namespace discretization {

/**
 * @brief Computes the Galerkin LSE for the Hodge Laplacian of the whitney one
 * form
 *
 * @f[
 *    \Delta_1 u = -\mathbf{curl}_{\Gamma} \circ \text{curl}_{\Gamma} +
 *     \mathbf{grad}_{\Gamma} \circ \text{div}_{\Gamma} \\
 *     -\Delta_1 u + k^2 u =  f
 * @f]
 *
 * Basis functions are the Whitney 1-forms and the barycentric
 * coordinate functions
 *
 * @note Only triangular meshes are supported
 *
 */
class WhitneyOneHodgeLaplace {
 public:
  /**
   * @brief Constructor
   * initializes basic mesh (Octaeder with radius 1.0)
   * initializes zerovalued function f
   *
   */
  WhitneyOneHodgeLaplace()
      : coo_matrix_(lf::assemble::COOMatrix<double>(1, 1)) {
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

    // create basic function
    auto f = [](Eigen::Matrix<double, 3, 1> x) -> Eigen::Matrix<double, 3, 1> {
      return Eigen::Matrix<double, 3, 1>::Zero();
    };
    f_ = f;
  }

  /**
   * @brief Computes the Galerkin LSE
   *
   * The Galerkin matrix will be accessable with GetGalerkinMatrix()
   * The load vector will be accessable with GetLoadVector()
   *
   */
  void Compute() {
    // create element matrix provider
    projects::hldo_sphere::assemble::WhitneyOneCurlCurlMatrixProvider
        matrix_curl_provider;
    projects::hldo_sphere::assemble::WhitneyOneGradMatrixProvider
        matrix_grad_provider;
    projects::hldo_sphere::assemble::MassMatrixProvider matrix_mass_provider;

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
    lf::assemble::COOMatrix<double> coo_A_21(n_dofs_curl, n_dofs_bary);
    coo_A_21.setZero();
    lf::assemble::COOMatrix<double> coo_A_22(n_dofs_bary, n_dofs_bary);
    coo_A_22.setZero();

    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
        0, dof_handler_curl, dof_handler_curl, matrix_curl_provider, coo_A_11);

    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
        0, dof_handler_bary, dof_handler_curl, matrix_grad_provider, coo_A_21);

    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
        0, dof_handler_bary, dof_handler_bary, matrix_mass_provider, coo_A_22);

    // create full matrix
    lf::assemble::COOMatrix<double> full_matrix(n_dofs_curl + n_dofs_bary,
                                                n_dofs_curl + n_dofs_bary);
    full_matrix.setZero();

    // iterate over all triplets of the previously computed matrice and add up
    // entries
    const std::vector<Eigen::Triplet<double>> triplets_A_11 =
        coo_A_11.triplets();
    const std::vector<Eigen::Triplet<double>> triplets_A_21 =
        coo_A_21.triplets();
    const std::vector<Eigen::Triplet<double>> triplets_A_22 =
        coo_A_22.triplets();

    // Add A_11
    for (Eigen::Triplet<double> triplet : triplets_A_11) {
      int col = triplet.col();
      int row = triplet.row();
      double val = triplet.value();
      full_matrix.AddToEntry(row, col, val);
    }

    // Add A_12 and A_21
    for (Eigen::Triplet<double> triplet : triplets_A_21) {
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

    // create element vector provider
    projects::hldo_sphere::assemble::WhitneyOneVectorProvider<double>
        vector_provider(f_);

    // create load vector
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi(n_dofs_curl);
    phi.setZero();

    // assemble the global vector over entities with codim=0:
    AssembleVectorLocally(0, dof_handler_curl, vector_provider, phi);

    // create full vector
    Eigen::Matrix<double, Eigen::Dynamic, 1> full_vec(n_dofs_curl +
                                                      n_dofs_bary);
    full_vec.setZero();
    full_vec.head(n_dofs_curl) = phi;

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
   * @brief Sets the load function
   * @param f load function
   *
   * @f[
   *    \Delta_1 u = -\mathbf{curl}_{\Gamma} \circ \text{curl}_{\Gamma} +
   *     \mathbf{grad}_{\Gamma} \circ \text{div}_{\Gamma} \\
   *     -\Delta_1 u + k^2 u =  f
   * @f]
   */
  void SetLoadFunction(std::function<Eigen::Matrix<double, 3, 1>(
                           const Eigen::Matrix<double, 3, 1> &)> &f) {
    f_ = f;
  }

  /**
   * @brief returns the Loadvector
   *
   * This is the righthandside of the LSE
   *
   * @note The loadvector must be computed with Compute() before calling this
   * function
   *
   */
  Eigen::Matrix<double, Eigen::Dynamic, 1> GetLoadVector() { return phi_; }

  /**
   * @brief returns the Galerkin Matrix
   *
   * This is the Matrix of the LSE
   *
   * @note The Galerkin matrix must be computed with Compute() before calling
   * this funciton
   *
   */
  lf::assemble::COOMatrix<double> GetGalerkinMatrix() { return coo_matrix_; }

 private:
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  std::function<Eigen::Matrix<double, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      f_;
  lf::assemble::COOMatrix<double> coo_matrix_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_;
};

}  // namespace discretization

}  // namespace projects::hldo_sphere

#endif  // THESIS_DISCRETIZATION_WHITNEY_ONE_HODGE_LAPLACE_H
