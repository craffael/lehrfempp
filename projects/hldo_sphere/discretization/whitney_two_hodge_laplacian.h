#ifndef THESIS_DISCRETIZATION_WHITNEY_TWO_HODGE_LAPLACE_H
#define THESIS_DISCRETIZATION_WHITNEY_TWO_HODEE_LAPLACE_H

/**
 * @file whitney_two_hodge_laplace.h
 * @brief Class to discretise the second hodge laplacian on
 * a given mesh with a given loadfunction.
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
#include <rot_whitney_one_div_matrix_provider.h>
#include <sphere_triag_mesh_builder.h>
#include <whitney_one_mass_matrix_provider.h>
#include <whitney_one_vector_provider.h>
#include <whitney_two_vector_provider.h>

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere {

namespace discretization {

/**
 * @brief Computes the Galerkin LSE for the Hodge Laplacian of the whitney two
 * form
 *
 * @f[
 *   \Delta_2 = \text{div}_{\Gamma} \circ \mathbf{grad}_{\Gamma} \\
 *   \Delta_2 u + k^2 u =   f
 * @f]
 *
 * Basis functions used are the rotated Whitney 1-forms and the cellwise
 * constant functions
 *
 * @note Only triangular meshes are supported
 *
 */
class WhitneyTwoHodgeLaplace {
 public:
  /**
   * @brief Constructor
   * initializes basic mesh (Octaeder with radius 1.0)
   * initializes zerovalued function f
   *
   */
  WhitneyTwoHodgeLaplace()
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
    auto f = [](Eigen::Matrix<double, 3, 1> x) -> double { return 0; };
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
    projects::hldo_sphere::assemble::WhitneyOneMassMatrixProvider
        matrix_dot_provider;

    projects::hldo_sphere::assemble::RotWhitneyOneDivMatrixProvider
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
    const std::vector<Eigen::Triplet<double>> triplets_A_11 =
        coo_A_11.triplets();

    const std::vector<Eigen::Triplet<double>> triplets_A_12 =
        coo_A_12.triplets();

    // Add A_11
    for (Eigen::Triplet<double> triplet : triplets_A_11) {
      int col = triplet.col();
      int row = triplet.row();
      double val = triplet.value();
      full_matrix.AddToEntry(row, col, val);
    }

    // Add A_12 and A_21
    for (Eigen::Triplet<double> triplet : triplets_A_12) {
      int row = triplet.row();
      int col = triplet.col() + n_dofs_div;
      double val = triplet.value();

      // Add A_12
      full_matrix.AddToEntry(row, col, val);

      // Add A_21 positive because we compute the laplacian here and not the
      // negative laplacian
      full_matrix.AddToEntry(col, row, val);
    }

    coo_matrix_ = full_matrix;

    // create element vector provider
    projects::hldo_sphere::assemble::WhitneyTwoVectorProvider vector_provider(
        f_);

    // create load vector
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi(n_dofs_const);
    phi.setZero();

    // assemble the global vector over entities with codim=0:
    lf::assemble::AssembleVectorLocally(0, dof_handler_const, vector_provider,
                                        phi);

    // create full vector
    Eigen::Matrix<double, Eigen::Dynamic, 1> full_vec(n_dofs_div +
                                                      n_dofs_const);
    full_vec.setZero();
    full_vec.tail(n_dofs_const) = phi;

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
   *   \Delta_2 = \text{div}_{\Gamma} \circ \mathbf{grad}_{\Gamma} \\
   *   \Delta_2 u + k^2 u =   f
   * @f]
   */
  void SetLoadFunction(
      std::function<double(const Eigen::Matrix<double, 3, 1> &)> &f) {
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
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> f_;
  lf::assemble::COOMatrix<double> coo_matrix_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_;
};

}  // namespace discretization

}  // namespace projects::hldo_sphere

#endif  // THESIS_DISCRETIZATION_WHITNEY_TWO_HODGE_LAPLACE_H
