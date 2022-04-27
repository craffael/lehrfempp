#ifndef THESIS_ASSEMBLE_WHITNEY_ONE_MASS_MATRIX_PROVIDER_H
#define THESIS_ASSEMBLE_WHITNEY_ONE_MASS_MATRIX_PROVIDER_H

/**
 * @file whitney_one_matrix_provider.h
 */

#include <lf/mesh/entity.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace assemble {

/**
 * @brief Element matrix provider for the whitney one mass matrix
 * that is the dot product of two whitney one basis functions
 *
 * @f[
 * \int_{\Omega} \bm{u} \cdot \bm{v} dx
 * @f]
 *
 * As basis functions, the Whitney 1-forms, surface edge elements are
 * used
 *
 * @note Only triangular meshes are supported
 *
 */
class WhitneyOneMassMatrixProvider {
 public:
  /**
   * @brief Constructor
   */
  WhitneyOneMassMatrixProvider(){};

  /**
   * @brief Compute the element matrix for a given triangle of a mesh
   * @param entity The mesh triangle to compute the element matrix for
   * @returns The 3 by 3 element matrix of the triangle
   *
   * @note Only triangular cells are supported
   */
  Eigen::MatrixXd Eval(const lf::mesh::Entity &entity) const;

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
};

}  // end namespace assemble

}  // namespace projects::hldo_sphere

#endif  // THESIS_ASSEMBLE_WHITNEY_ONE_MASS_MATRIX_PROVIDER_H
