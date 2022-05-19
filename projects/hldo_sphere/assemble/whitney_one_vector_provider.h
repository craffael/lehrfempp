#ifndef THESIS_ASSEMBLE_WHITNEY_ONE_VECTOR_PROVIDER_H
#define THESIS_ASSEMBLE_WHITNEY_ONE_VECTOR_PROVIDER_H

/**
 * @file whitney_one_vector_provider.h
 *
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
 * @brief Element vector provider for whitney one forms
 *
 * The linear form is given by
 * @f[
 *  \int\limits_K f \cdot v \,\mathrm{d}x, \quad f, v \in
 * \mathbf{H}(\text{curl}_{\Gamma}, \partial\mathbf{S})
 * @f]
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
class WhitneyOneVectorProvider {
 public:
  /**
   * @brief Constructor
   * @param f a tangential vector field on the sphere
   */
  WhitneyOneVectorProvider(
      std::function<Eigen::Vector3d(const Eigen::Vector3d &)> f);

  /**
   * @brief Compute the element vector for a given triangle off the mesh
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
  const std::function<Eigen::Vector3d(const Eigen::Vector3d &)> f_;
};

}  // end namespace assemble

}  // namespace projects::hldo_sphere

#endif  // THESIS_ASSEMBLE_WHITNEY_ONE_VECTOR_PROVIDER_H
