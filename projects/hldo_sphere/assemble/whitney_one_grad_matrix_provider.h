#ifndef THESIS_ASSEMBLE_WHITNEY_ONE_GRAD_MATRIX_PROVIDER_H
#define THESIS_ASSEMBLE_WHITNEY_ONE_GRAD_MATRIX_PROVIDER_H

/**
 * @file whitney_one_grad_matrix_provider.h
 */

#include <lf/mesh/entity.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace assemble {

/**
 * @brief Element matrix provider for the bilinear form
 *
 * @f[
 * \int grad_{\Gamma}(u) v dx
 * @f]
 *
 * Basis functions are the Whitney 1-forms, surface edge elements for v
 * and the barycentric basis function for u
 *
 * @note Only triangluar meshes are supported
 *
 */
class WhitneyOneGradMatrixProvider {
 public:
  /**
   * @brief Constructor
   */
  WhitneyOneGradMatrixProvider(){};

  /**
   * @brief Compute the element matrix for some cell of a mesh
   * @param entity The mesh cell to compute the element matrix for
   * @returns The 3 by 3 element matrix of the cell
   *
   * @note Only triangluar cells are supported
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

#endif  // THESIS_ASSEMBLE_WHITNEY_ONE_GRAD_MATRIX_PROVIDER_H
