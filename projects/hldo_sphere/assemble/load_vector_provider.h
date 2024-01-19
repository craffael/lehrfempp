#ifndef HLDO_SPHERE_ASSEMBLE_LOAD_VECTOR_PROVIDER_H
#define HLDO_SPHERE_ASSEMBLE_LOAD_VECTOR_PROVIDER_H

/**
 * @file load_vector_provider.h
 */

#include <lf/mesh/entity.h>
#include <lf/mesh/utils/mesh_data_set.h>
#include <lf/quad/quad.h>
#include <lf/uscalfe/lagr_fe.h>

#include <Eigen/Dense>
#include <functional>
#include <utility>

namespace projects::hldo_sphere {

namespace assemble {

/**
 *
 * @brief Element Vector Provider for scalar valued load functions using
 * picewise linear barycentric basis functions.
 *
 * The element matrix provider works in a 3 dimensional world with 2
 * dimensional triangular cells.
 *
 * @tparam SCALAR either double or std::complex
 *
 * The locally evaluated linear form is
 * @f[
 * v \mapsto \int\limits_{K}\ f(\mathbf{x}) \cdot v(\mathbf{x}) \, dx, \quad f,
 * v \in H^1(\Omega)
 * @f]
 *
 * Details regarding the mathematical derivations can be found in the thesis
 * `Hodge-Laplacians and Dirac Operators on the Surface of the 3-Sphere`
 * section 4.2.3.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_VECTOR_PROVIDER of the function
 * lf::assemble::AssembleVectorLocally().
 *
 * @note Only triangular meshes are supported
 *
 */
template <typename SCALAR>
class LoadVectorProvider {
 public:
  /**
   *
   * @brief Constructor setting the load function
   *
   * @param f a SCALAR valued functor defined on the Mesh on which Eval()
   * will be called
   *
   */
  LoadVectorProvider(const std::function<SCALAR(const Eigen::Vector3d &)> &f)
      : f_(f) {}

  /**
   * @brief Compute the element vector for some cell on the mesh
   * @param entity the mesh cell to compute the element vector for
   * @returns The element vector of the given cell
   *
   * @note Only triangular cells are supported
   */
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> Eval(
      const lf::mesh::Entity &entity) const {
    const auto *const geom = entity.Geometry();

    // Only triangles are supported
    LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                  "Unsupported cell type " << entity.RefEl());

    // Compute the global vertex coordinates
    Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> vertices =
        geom->Global(entity.RefEl().NodeCoords());

    // define quad rule with sufficiantly high degree since the
    // baricentric coordinate functions have degree 1
    lf::quad::QuadRule quadrule{lf::quad::make_TriaQR_EdgeMidpointRule()};

    const lf::uscalfe::FeLagrangeO1Tria<SCALAR> hat_func;
    const auto lambda_hat = [&](Eigen::Vector2d x_hat)
        -> Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> {
      return hat_func.EvalReferenceShapeFunctions(x_hat);
    };

    // returns evaluation of @f$ \lambda @f$ for at a given point on the
    // reference triangle f is evaluated at the global coordinates
    const auto f_tilde_hat =
        [&](Eigen::Vector2d x_hat) -> Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> {
      Eigen::MatrixXd x = geom->Global(x_hat);
      return lambda_hat(x_hat) * f_(x);
    };

    // Compute the elements of the vector with given quadrature rule
    Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> element_vector =
        Eigen::Matrix<SCALAR, 3, 1>::Zero();
    const Eigen::MatrixXd points = quadrule.Points();
    const Eigen::VectorXd weights =
        (geom->IntegrationElement(points).array() * quadrule.Weights().array())
            .matrix();
    for (lf::base::size_type n = 0; n < quadrule.NumPoints(); ++n) {
      const Eigen::Matrix<SCALAR, 3, 1> f_eval = f_tilde_hat(points.col(n));
      element_vector += weights[n] * f_eval;
    }
    return element_vector;
  }

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
  const std::function<SCALAR(const Eigen::Vector3d &)> f_;
};

}  // end namespace assemble

}  // namespace projects::hldo_sphere

#endif  // HLDO_SPHERE_ASSEMBLE_LOAD_VECTOR_PROVIDER_H
