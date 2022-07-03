#ifndef HLDO_SPHERE_OPERATORS_DIRAC_OPERATOR_H
#define HLDO_SPHERE_OPERATORS_DIRAC_OPERATOR_H

/**
 * @file dirac_operator.h
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
#include <rot_whitney_one_div_matrix_provider.h>
#include <sphere_triag_mesh_builder.h>
#include <whitney_one_grad_matrix_provider.h>
#include <whitney_one_vector_provider.h>
#include <whitney_two_vector_provider.h>

#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere {

namespace operators {

using complex = std::complex<double>;

/**
 * @brief Computes the Galerkin LSE for the Dirac Operator and the load vector
 *
 * @f[
 *   \begin{pmatrix}
 *       & \int\limits_{\partial \mathbf{S}} \mathbf{u} \cdot
 * \mathbf{grad}_{\Gamma} v \, dS & \\ \int_{\partial \mathbf{S}}
 * \mathbf{grad}_{\Gamma} u \cdot \mathbf{v} \, dS & & \int\limits_{\partial
 * \mathbf{S}} \mu \ \ \text{curl}_{\Gamma} \mathbf{v} \, dS  \\ &
 * \int\limits_{\partial \mathbf{S}} \text{curl}_{\Gamma} \mathbf{u} \ \ \nu
 * \, dS & \end{pmatrix}
 *   ,
 *   \begin{pmatrix}
 *      \int\limits_{\partial \mathbf{S}} f v \, dS \\
 *      \int\limits_{\partial \mathbf{S}} \mathbf{f} \cdot \mathbf{v} \, dS \\
 *      \int\limits_{\partial \mathbf{S}} \varphi \nu \, dS
 *   \end{pmatrix}
 *   @f]
 *
 *
 * @note Only triangular meshes are supported
 *
 */
class DiracOperator {
 public:
  /**
   * @brief Constructor
   * creates basic mesh (Octaeder with radius 1.0)
   * and zero valued functions f
   *
   */
  DiracOperator() : coo_matrix_(lf::assemble::COOMatrix<complex>(1, 1)) {
    // create mesh factory for basic mesh only because empty values do not
    // complile
    std::unique_ptr<lf::mesh::MeshFactory> factory =
        std::make_unique<lf::mesh::hybrid2d::MeshFactory>(3);
    std::shared_ptr<projects::hldo_sphere::mesh::SphereTriagMeshBuilder>
        sphere = std::make_shared<
            projects::hldo_sphere::mesh::SphereTriagMeshBuilder>(
            std::move(factory));
    sphere->setRefinementLevel(0);
    sphere->setRadius(1);
    mesh_p_ = sphere->Build();

    // create basic function everyting 0 valued by default
    auto f_0 = [](Eigen::Matrix<double, 3, 1> x) -> complex { return 0; };
    auto f_1 =
        [](Eigen::Matrix<double, 3, 1> x) -> Eigen::Matrix<complex, 3, 1> {
      return Eigen::Matrix<complex, 3, 1>::Zero();
    };
    auto f_2 = [](Eigen::Matrix<double, 3, 1> x) -> complex { return 0; };
    f0_ = f_0;
    f1_ = f_1;
    f2_ = f_2;
  }

  /**
   * @brief Computes the Galerkin LSE
   *
   * The Galerkin matrix can be accessable with GetGalerkinMatrix()
   * The load vector can be accessable with GetLoadVector()
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
  void SetMesh(std::shared_ptr<const lf::mesh::Mesh> mesh_p);

  /**
   * @brief Sets the load functions
   * @param f0 load function in @f$ L^2 @f$
   * @param f1 load functions in @f$ L^2_t @f$ vector valued
   * @param f2 load functions in @f$ L^2 @f$
   */
  void SetLoadFunctions(
      std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f0,
      std::function<
          Eigen::Matrix<complex, 3, 1>(const Eigen::Matrix<double, 3, 1> &)>
          f1,
      std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f2) {
    f0_ = f0;
    f1_ = f1;
    f2_ = f2;
  }

  /**
   * @brief returns the Loadvector
   *
   * This is the righthandside of the LSE
   *
   * @note The loadvector must be computed with Compute() before calling
   * this function
   *
   */
  Eigen::Matrix<complex, Eigen::Dynamic, 1> GetLoadVector() { return phi_; }

  /**
   * @brief returns the Galerkin Matrix
   *
   * This is the Matrix of the LSE
   *
   * @note The Galerkin matrix must be computed with Compute() before
   * calling this funciton
   *
   */
  lf::assemble::COOMatrix<complex> GetGalerkinMatrix() { return coo_matrix_; }

 private:
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;
  std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f0_;
  std::function<Eigen::Matrix<complex, 3, 1>(
      const Eigen::Matrix<double, 3, 1> &)>
      f1_;
  std::function<complex(const Eigen::Matrix<double, 3, 1> &)> f2_;
  lf::assemble::COOMatrix<complex> coo_matrix_;
  Eigen::Matrix<complex, Eigen::Dynamic, 1> phi_;
};

}  // namespace operators

}  // namespace projects::hldo_sphere

#endif  // HLDO_SPHERE_OPERATORS_WHITNEY_ONE_HODGE_LAPLACE_H
