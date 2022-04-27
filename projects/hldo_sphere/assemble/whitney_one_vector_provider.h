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
 *  \int_K\ f \cdot v \,\mathrm{d}x, \quad f, v \in \bm{H}(div_{\Gamma},
 * \partial\mathbb{S})
 * @f]
 * The basis functions used for the discretization are the
 * whitney one forms
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
