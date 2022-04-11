#ifndef THESIS_ASSEMBLE_WHITNEY_TWO_ELEMENT_VECTOR_PROVIDER_H
#define THESIS_ASSEMBLE_WHITNEY_TWO_ELEMENT_VECTOR_PROVIDER_H

/**
 * @file whitney_two_element_vector_provider,h
 * @brief Assembly routines for the right hand side function of
 * the thrid Dirac opeartor
 */

#include <lf/mesh/entity.h>
#include <lf/mesh/utils/mesh_data_set.h>
#include <lf/quad/quad_rule.h>

#include <Eigen/Dense>
#include <functional>
#include <utility>

namespace projects::hldo_sphere {

namespace assemble {

/**
 * @brief Element vector provider for piecewise constant whitney two forms
 *
 * @note Only triangular meshes are supported
 *
 * The linear form is given by
 * @f[
 *  \int_K\!  f q \,\mathrm{d}x, \quad f, q \in L^2(\partial \mathbb{S})
 * @f]
 */
class WhitneyTwoElementVectorProvider {
 public:
  /**
   * @brief Constructor
   * @param f A scalar valued function
   * @param quadrule The quadrature rule used to compute the integral in the
   * linear form
   */
  WhitneyTwoElementVectorProvider(
      std::function<double(const Eigen::Vector3d &)> f,
      lf::quad::QuadRule quadrule);

  /**
   * @brief Compute the element vector for some entity of the mesh
   * @param entity The mesh entity to compute the element vector for
   * @returns The element vector of the given entity
   *
   * @note This function currently only works for triangles
   */
  Eigen::VectorXd Eval(const lf::mesh::Entity &entity) const;

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
  const std::function<double(const Eigen::Vector3d &)> f_;
  const lf::quad::QuadRule quadrule_;
};

}  // end namespace assemble

}  // namespace projects::hldo_sphere

#endif  // THESIS_ASSEMBLE_WHITNEY_TWO_ELEMENT_VECTOR_PROVIDER_H
