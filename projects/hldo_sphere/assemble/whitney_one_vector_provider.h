#ifndef THESIS_ASSEMBLE_WHITNEY_ONE_VECTOR_PROVIDER_H
#define THESIS_ASSEMBLE_WHITNEY_ONE_VECTOR_PROVIDER_H

/**
 * @file whitney_one_vector_provider.h
 *
 */

#include <lf/base/base.h>
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
 * @brief Element vector provider for whitney one forms
 *
 * The linear form is given by
 * @f[
 *  (\mathbf{v}) \mapsto \int\limits_K \mathbf{f} \cdot \mathbf{v}
 * \,\mathrm{d}x, \quad \mathbf{f}, \mathbf{v} \in
 * \mathbf{H}(\text{curl}_{\Gamma}, \partial\mathbf{S})
 * @f]
 *
 * The element vector provider works in a 3 dimensional world with 2
 * dimensional triangular cells.
 *
 * Basis functions are the Whitney 1-forms, surface edge elements for v
 *
 * The whitney 1-forms, surface edge elements are associated with
 * edges and defined as
 *
 * @f[
 *  b_i = s_i (\lambda_i \mathbf{grad}_{\Gamma}(\lambda_{i+1}) - \lambda_{i+1}
 * \mathbf{grad}_{\Gamma}(\lambda_{i}))
 * @f]
 *
 * with @f$ \lambda_i @f$ barycentric basis functions and @f$ s_i @f$ is a sign
 * of the function based on the relative orientation of the edge in the mesh.
 * @note This class complies with the type requirements for the template
 * argument ENTITY_VECTOR_PROVIDER of the function
 * lf::assemble::AssembleVectorLocally().
 *
 * @note Only triangular meshes are supported
 */
template <typename SCALAR>
class WhitneyOneVectorProvider {
 public:
  /**
   * @brief Constructor
   *
   * @tparam SCALAR representing the scalar field, on which the codomain is
   * based on
   *
   * @param f a tangential vector field functor required to be defined on the
   * mesh
   */
  WhitneyOneVectorProvider(
      const std::function<Eigen::Matrix<SCALAR, 3, 1>(const Eigen::Vector3d &)>
          &f)
      : f_(std::move(f)) {}

  /**
   * @brief Compute the element vector for a given triangle off the mesh
   * @param entity The mesh triangle to compute the element vector for
   * @returns The element vector of the given triangle
   *
   * @note Only triangular cells are supported
   */
  Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> Eval(
      const lf::mesh::Entity &entity) const {
    LF_ASSERT_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                  "Only Triangles are supported no " << entity.RefEl());

    const auto *const geom = entity.Geometry();

    // Compute the global vertex coordinates
    Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());

    // define quad rule with sufficiantly high degree
    lf::quad::QuadRule quadrule{lf::quad::make_TriaQR_EdgeMidpointRule()};

    // Construct the whitney one form basis functions
    const lf::uscalfe::FeLagrangeO1Tria<SCALAR> hat_func;
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> ref_grads =
        hat_func.GradientsReferenceShapeFunctions(Eigen::VectorXd::Zero(2))
            .transpose();
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> J_inv_trans =
        geom->JacobianInverseGramian(Eigen::VectorXd::Zero(2));
    const Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> grads =
        J_inv_trans * ref_grads;

    // get edge orientations
    auto edgeOrientations = entity.RelativeOrientations();
    Eigen::Vector3d s;
    s(0) = lf::mesh::to_sign(edgeOrientations[0]);
    s(1) = lf::mesh::to_sign(edgeOrientations[1]);
    s(2) = lf::mesh::to_sign(edgeOrientations[2]);

    const auto b_hat = [&](Eigen::Vector2d x_hat)
        -> Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> {
      Eigen::Matrix<SCALAR, 3, 3> b_hat;
      Eigen::Matrix<SCALAR, 3, 1> lambda_hat =
          hat_func.EvalReferenceShapeFunctions(x_hat);
      for (int i = 0; i < 3; i++) {
        int ip1 = (i + 1) % 3;

        b_hat.col(i) = s(i) * (lambda_hat(i) * grads.col(ip1) -
                               lambda_hat(ip1) * grads.col(i));
      }

      return b_hat;
    };

    // returns evaluation of @f$\textbf{f} * \textbf{b}@f$ for at a given point
    // on the reference triangle f is evaluated at the global coordinates of the
    // triangle radially projected on the sphere
    const auto f_tilde_hat =
        [&](Eigen::Vector2d x_hat) -> Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> {
      Eigen::Vector3d x = geom->Global(x_hat);
      Eigen::Matrix<SCALAR, Eigen::Dynamic, Eigen::Dynamic> b_hat_out =
          b_hat(x_hat);
      Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> result =
          b_hat(x_hat).transpose() * f_(x);
      return result;
    };

    // Compute the elements of the vector with given quadrature rule
    Eigen::Matrix<SCALAR, 3, 1> element_vector;
    element_vector.setZero();
    const Eigen::MatrixXd points = quadrule.Points();
    const Eigen::VectorXd weights =
        (geom->IntegrationElement(points).array() * quadrule.Weights().array())
            .matrix();
    for (lf::base::size_type n = 0; n < quadrule.NumPoints(); ++n) {
      const Eigen::Matrix<SCALAR, Eigen::Dynamic, 1> f_eval =
          f_tilde_hat(points.col(n));
      element_vector += weights[n] * f_eval;
    }
    return element_vector;
  }

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
  const std::function<Eigen::Matrix<SCALAR, 3, 1>(const Eigen::Vector3d &)> f_;
};

}  // end namespace assemble

}  // namespace projects::hldo_sphere

#endif  // THESIS_ASSEMBLE_WHITNEY_ONE_VECTOR_PROVIDER_H
