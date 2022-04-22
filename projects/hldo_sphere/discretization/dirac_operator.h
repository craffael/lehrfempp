#ifndef THESIS_DISCRETIZATION_DIRAC_OPERATOR_H
#define THESIS_DISCRETIZATION_DIRAC_OPERATOR_H

/**
 * @file dirac_operator.h
 * @brief Class to discretise the dirac_operator.
 */

#include <curl_element_vector_provider.h>
#include <lf/assemble/assembler.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/mesh_interface.h>
#include <lf/quad/quad.h>
#include <load_element_vector_provider.h>
#include <rot_w_one_form_div_element_matrix_provider.h>
#include <sphere_triag_mesh_builder.h>
#include <whitney_one_form_grad_element_matrix_provider.h>
#include <whitney_two_element_vector_provider.h>

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere {

namespace discretization {

/**
 * @brief Computes the Galerkin LSE for the Dirac Operator
 *
 * @note Only triangular meshes are supported
 *
 */
class DiracOperator {
 public:
  /**
   * @brief Constructor
   * creates basic mesh (Octaeder with radius 1.0)
   * creates uniform dofhandler
   * creates zerovalued function f
   *
   */
  DiracOperator()
      : coo_matrix_(lf::assemble::COOMatrix<std::complex<double>>(1, 1)) {
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
   * @brief Computes the Galerkin LSE
   *
   * The Galerkin matrix will be accessable with `get_galerkin_matrix`
   * The load vector will be accessable with `get_load_vector`
   *
   */
  void Compute() {
    // create element matrix provider
    projects::hldo_sphere::assemble::WhitneyOneFormGradElementMatrixProvider
        matrix_grad_provider;
    projects::hldo_sphere::assemble::RotWOneFormDivElementMatrixProvider
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

    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
        0, dof_handler_vert, dof_handler_edge, matrix_grad_provider, coo_A_21);

    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
        0, dof_handler_cell, dof_handler_edge, matrix_div_provider, coo_A_23_m);

    // create full matrix
    lf::assemble::COOMatrix<std::complex<double>> full_matrix(
        n_dofs_vert + n_dofs_edge + n_dofs_cell,
        n_dofs_vert + n_dofs_edge + n_dofs_cell);

    full_matrix.setZero();

    // iterate over all triplets of the previously computed matrice and add up
    // entries
    const std::vector<Eigen::Triplet<double>> triplets_A_21 =
        coo_A_21.triplets();
    const std::vector<Eigen::Triplet<double>> triplets_A_23_m =
        coo_A_23_m.triplets();

    // Add A_12 and A_21
    for (Eigen::Triplet<double> triplet : triplets_A_21) {
      int col = triplet.col();
      int row = triplet.row();
      std::complex<double> val = std::complex<double>(triplet.value(), 0.);
      // A_21
      full_matrix.AddToEntry(row + n_dofs_vert, col, val);
      // A_12
      full_matrix.AddToEntry(col, row + n_dofs_vert, val);
    }

    // Add A_23 and A_32
    for (Eigen::Triplet<double> triplet : triplets_A_23_m) {
      int col = triplet.col();
      int row = triplet.row();
      std::complex<double> val = std::complex<double>(-triplet.value(), 0.);

      // Add A_23
      full_matrix.AddToEntry(row + n_dofs_vert, col + n_dofs_vert + n_dofs_edge,
                             val);

      // Add A_32
      full_matrix.AddToEntry(col + n_dofs_vert + n_dofs_edge, row + n_dofs_vert,
                             val);
    }

    coo_matrix_ = full_matrix;

    // define quad rule with sufficiantly high degree since the
    // baricentric coordinate functions and whitney one form basis functions
    // have degree 1
    lf::quad::QuadRule quadrule{lf::quad::make_TriaQR_EdgeMidpointRule()};

    // create element vector provider
    projects::hldo_sphere::assemble::LoadElementVectorProvider
        vector_provider_0(f0_);

    projects::hldo_sphere::assemble::CurlElementVectorProvider
        vector_provider_1(f1_, quadrule);

    projects::hldo_sphere::assemble::WhitneyTwoElementVectorProvider
        vector_provider_2(f2_, quadrule);

    // create load vector
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi0(n_dofs_vert);
    phi0.setZero();
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi1(n_dofs_edge);
    phi1.setZero();
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi2(n_dofs_cell);
    phi2.setZero();

    // assemble the global vector over entities with codim=0:
    AssembleVectorLocally(0, dof_handler_vert, vector_provider_0, phi0);
    AssembleVectorLocally(0, dof_handler_edge, vector_provider_1, phi1);
    AssembleVectorLocally(0, dof_handler_cell, vector_provider_2, phi2);

    // create full vector
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> full_vec(
        n_dofs_vert + n_dofs_edge + n_dofs_cell);
    full_vec.setZero();
    full_vec.head(n_dofs_vert) = phi0.cast<std::complex<double>>();
    full_vec.segment(n_dofs_vert, n_dofs_edge) =
        phi1.cast<std::complex<double>>();
    full_vec.tail(n_dofs_cell) = phi2.cast<std::complex<double>>();

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
   * @brief returns the Loadvector
   *
   * This is the righthandside of the LSE
   *
   * @note The loadvector must be computed with `Compute` before calling
   * this function
   *
   */
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> GetLoadVector() {
    return phi_;
  }

  /**
   * @brief returns the Galerkin Matrix
   *
   * This is the Matrix of the LSE
   *
   * @note The Galerkin matrix must be computed with `Compute` before
   * calling this funciton
   *
   */
  lf::assemble::COOMatrix<std::complex<double>> GetGalerkinMatrix() {
    return coo_matrix_;
  }

 private:
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> f0_;
  std::function<Eigen::Matrix<double, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      f1_;
  std::function<double(const Eigen::Matrix<double, 3, 1> &)> f2_;
  lf::assemble::COOMatrix<std::complex<double>> coo_matrix_;
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> phi_;
};

}  // namespace discretization

}  // namespace projects::hldo_sphere

#endif  // THESIS_DISCRETIZATION_WHITNEY_ONE_HODGE_LAPLACE_H
