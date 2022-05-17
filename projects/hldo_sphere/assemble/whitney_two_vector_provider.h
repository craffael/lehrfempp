#ifndef THESIS_ASSEMBLE_WHITNEY_TWO_VECTOR_PROVIDER_H
#define THESIS_ASSEMBLE_WHITNEY_TWO_VECTOR_PROVIDER_H

/**
 * @file whitney_two_vector_provider,h
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
class WhitneyTwoVectorProvider {
 public:
  /**
   * @brief Constructor
   * @param f A scalar valued function defined on the surface of the sphere
   * linear form
   */
  WhitneyTwoVectorProvider(std::function<double(const Eigen::Vector3d &)> f);

  /**
   * @brief Compute the element vector for some trinagle of the mesh
   * @param entity The mesh triangle to compute the element vector for
   * @returns The element vector of the given triangle
   *
   * @note Only triangular cells are supported
   */
  Eigen::VectorXd Eval(const lf::mesh::Entity &entity) const;

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
  const std::function<double(const Eigen::Vector3d &)> f_;
};

}  // end namespace assemble

}  // namespace projects::hldo_sphere

#endif  // THESIS_ASSEMBLE_WHITNEY_TWO_VECTOR_PROVIDER_H
