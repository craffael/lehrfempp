#ifndef HLDO_SPHERE_ASSEMBLE_WHITNEY_TWO_VECTOR_PROVIDER_H
#define HLDO_SPHERE_ASSEMBLE_WHITNEY_TWO_VECTOR_PROVIDER_H

/**
 * @file whitney_two_vector_provider,h
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
 * @brief Element vector provider for piecewise constant whitney two forms
 *
 * @tparam SCALAR codomain type of the laodfunction
 *
 * The linear form is given by
 * @f[
 *  \int\limits_K\!  f\, q \,\mathrm{d}x, \quad f, q \in L^2(\partial
 * \mathbb{S})
 * @f]
 *
 * As basis functions for q we use the cellwise constant functions.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_VECTOR_PROVIDER of the function
 * lf::assemble::AssembleVectorLocally().
 *
 * @note Only triangular meshes are supported
 */
template <typename SCALAR>
class WhitneyTwoVectorProvider {
 public:
  /**
   * @brief Constructor
   * @param f A scalar valued function defined on the surface of the sphere
   * linear form
   */
  WhitneyTwoVectorProvider(
      const std::function<SCALAR(const Eigen::Vector3d &)> &f)
      : f_(std::move(f)) {}

  /**
   * @brief Compute the element vector for some trinagle of the mesh
   * @param entity The mesh triangle to compute the element vector for
   * @returns The element vector of the given triangle
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
    Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());

    // define quad rule with sufficiantly high degree since the
    // baricentric coordinate functions have degree 1
    lf::quad::QuadRule quadrule{lf::quad::make_TriaQR_EdgeMidpointRule()};

    // Compute the elements of the vector with given quadrature rule
    const Eigen::MatrixXd points = geom->Global(quadrule.Points());
    const Eigen::VectorXd weights =
        (geom->IntegrationElement(quadrule.Points()).array() *
         quadrule.Weights().array())
            .matrix();

    SCALAR sum = 0;
    for (lf::base::size_type n = 0; n < quadrule.NumPoints(); ++n) {
      Eigen::VectorXd x = points.col(n);
      sum += weights[n] * f_(x);
    }

    Eigen::Matrix<SCALAR, 1, 1> element_vector;
    element_vector(0) = sum;

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

#endif  // HLDO_SPHERE_ASSEMBLE_WHITNEY_TWO_VECTOR_PROVIDER_H
