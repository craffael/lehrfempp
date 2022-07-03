#ifndef HLDO_SPHERE_OPERATORS_WHITNEY_ZERO_HODGE_LAPLACE_H
#define HLDO_SPHERE_OPERATORS_WHITNEY_ZERO_HODEE_LAPLACE_H

/**
 * @file whitney_zero_hodge_laplace.h
 * @brief An element matrix provider for one vector valued piecewise linear
 * basis functions and one piecewise constant basis function
 */

#include <laplace_matrix_provider.h>
#include <lf/assemble/assembler.h>
#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/hybrid2d/mesh.h>
#include <lf/mesh/hybrid2d/mesh_factory.h>
#include <lf/mesh/mesh_interface.h>
#include <load_vector_provider.h>
#include <sphere_triag_mesh_builder.h>

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere {

namespace operators {

/**
 * @brief Computes the Galerkin LSE for the Hodge Laplacian of the whitney zero
 * form
 *
 * @f[
 *   -\Delta_0 = -\text{div}_{\Gamma} \circ \mathbf{grad}_{\Gamma} \\
 *   -\Delta_0 u  =   f
 * @f]
 *
 * As basis functions, the barycentric basis functions are used
 *
 * @note Only triangular meshes are supported
 *
 */
class WhitneyZeroHodgeLaplace {
 public:
  /**
   * @brief Basic Constructor
   *
   * initializes basic mesh (Octaeder with radius 1.0)
   * initializes zerovalued function f
   *
   */
  WhitneyZeroHodgeLaplace()
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
   * The load vector will be accessable with GetLoadVector
   *
   */
  void Compute() {
    // create element matrix provider
    projects::hldo_sphere::assemble::LaplaceMatrixProvider matrix_provider;

    const lf::assemble::DofHandler& dof_handler =
        lf::assemble::UniformFEDofHandler(mesh_p_,
                                          {{lf::base::RefEl::kPoint(), 1}});

    // create COO Matrix
    const lf::assemble::size_type n_dofs(dof_handler.NumDofs());

    lf::assemble::COOMatrix<double> coo_mat(n_dofs, n_dofs);
    coo_mat.setZero();

    lf::assemble::AssembleMatrixLocally<lf::assemble::COOMatrix<double>>(
        0, dof_handler, dof_handler, matrix_provider, coo_mat);

    // create element vector provider
    projects::hldo_sphere::assemble::LoadVectorProvider<double> vector_provider{
        f_};

    // create load vector
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi(n_dofs);
    phi.setZero();

    // assemble the global vector over entities with codim=0:
    AssembleVectorLocally(0, dof_handler, vector_provider, phi);

    // set the class attributes
    phi_ = phi;
    coo_matrix_ = coo_mat;
  }

  /**
   * @brief Sets the mesh for later computations
   *
   * @param mesh_p pointer to the mesh
   *
   * requries all cells in the mesh are triangles
   * requries mesh global dimension to be 3
   */
  void SetMesh(std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
    // check if cells are triagles
    for (const lf::mesh::Entity* tria : mesh_p->Entities(0)) {
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
   * @brief Sets the load function for later computations
   * @param f load function
   *
   * @f[
   *   -\Delta_0 = \text{div}_{\Gamma} \circ \mathbf{grad}_{\Gamma} \\
   *   -\Delta_0 u + k^2 u =   f
   * @f]
   */
  void SetLoadFunction(
      std::function<double(const Eigen::Matrix<double, 3, 1>&)>& f) {
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
  std::function<double(const Eigen::Matrix<double, 3, 1>&)> f_;
  lf::assemble::COOMatrix<double> coo_matrix_;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_;
};

}  // namespace operators

}  // namespace projects::hldo_sphere

#endif  // HLDO_SPHERE_OPERATORS_WHITNEY_ZERO_HODGE_LAPLACE_H
