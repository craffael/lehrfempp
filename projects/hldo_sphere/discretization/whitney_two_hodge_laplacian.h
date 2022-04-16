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
#include <sphere_triag_mesh_builder.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace discretization {

/**
 * @brief Computes the Galerkin LSE for the Hodge Laplacian of the whitney two
 * form
 *
 * @f[
 * -\Delta_2 u = div_{\Gamma} \circ grad_{\Gamma}
 * @f]
 *
 * Basis functions are the rotated Whitney 1-forms and the cellwise
 * constant functions
 *
 * @note Only triangular meshes are supported
 *
 */
class WhitneyTwoHodgeLaplace {
 public:
  /**
   * @brief Constructor
   * creates basic mesh (Octaeder with radius 1.0)
   * creates uniform dofhandler
   * creates zerovalued function f
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
    auto f = [](Eigen::VectorXd x) -> double { return 0; };
    f_ = f;
  }

  /**
   * @brief Computes the Galerkin LSE
   *
   * The Galerkin matrix will be accessable with `get_galerkin_matrix`
   * The load vector will be accessable with `get_load_vector`
   *
   */
  void Compute();

  /**
   * @brief Sets the mesh and creates dof_handler
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
   * @brief Sets the load function
   * @param f load function
   *
   * @f[
   * -\Delta_0 u = -div_{\Gamma} \circ \bm{grad}_{\Gamma} \li
   *  u := \Delta_0^{-1} f
   * @f]
   */
  void SetLoadFunction(std::function<double(const Eigen::Vector3d&)> f) {
    f_ = f;
  }

  /**
   * @brief returns the Loadvector
   *
   * This is the righthandside of the LSE
   *
   * @note The loadvector must be computed with `Compute` before calling this
   * function
   *
   */
  Eigen::VectorXd GetLoadVector() { return phi_; }

  /**
   * @brief returns the Galerkin Matrix
   *
   * This is the Matrix of the LSE
   *
   * @note The Galerkin matrix must be computed with `Compute` before calling
   * this funciton
   *
   */
  lf::assemble::COOMatrix<double> GetGalerkinMatrix() { return coo_matrix_; }

 private:
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  std::function<double(const Eigen::Vector3d&)> f_;
  lf::assemble::COOMatrix<double> coo_matrix_;
  Eigen::VectorXd phi_;
};

}  // namespace discretization

}  // namespace projects::hldo_sphere

#endif  // THESIS_DISCRETIZATION_WHITNEY_TWO_HODGE_LAPLACE_H
