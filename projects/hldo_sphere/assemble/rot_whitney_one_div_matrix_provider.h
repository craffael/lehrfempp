#ifndef THESIS_ASSEMBLE_ROT_WHITNEY_ONE_DIV_MATRIX_PROVIDER_H
#define THESIS_ASSEMBLE_ROT_WHITNEY_ONE_DIV_MATRIX_PROVIDER_H

/**
 * @file rot_whitney_one_div_matrix_provider.h
 */

#include <lf/mesh/entity.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace assemble {

/**
 * @brief Element matrix provider for rotated whitney one form
 *
 * @f[
 * \int div_{\Gamma}(v) u dx
 * @f]
 *
 * As basis functions, the rotatedrotated  Whitney 1-forms, surface edge
 * elements are used for v
 * and the cellwise constant basis functions are used for u
 *
 * @note Only triangular meshes are supported
 *
 */
class RotWhitneyOneDivMatrixProvider {
 public:
  /**
   * @brief Constructor
   */
  RotWhitneyOneDivMatrixProvider(){};

  /**
   * @brief Compute the element matrix for some triangle of a mesh
   * @param entity The mesh triangle to compute the element matrix for
   * @returns The 3 by 1 element matrix of the triangle
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

#endif  // THESIS_ASSEMBLE_ROT_WHITNEY_ONE_DIV_MATRIX_PROVIDER_H
