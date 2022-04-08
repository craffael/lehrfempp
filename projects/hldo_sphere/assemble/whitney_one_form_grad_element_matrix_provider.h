#ifndef THESIS_ASSEMBLE_WHITNEY_ONE_FORM_GRAD_ELEMENT_MATRIX_PROVIDER_H
#define THESIS_ASSEMBLE_WHITNEY_ONE_FORM_GRAD_ELEMENT_MATRIX_PROVIDER_H

/**
 * @file whitney_one_form_curl_element_matrix_provider.h
 * @brief An element matrix provider for a vector valued piecewise linear basis
 * function and a picewise linear function
 * @f[
 * \int grad_{\Gamma}(u) v dx
 * @f]
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
 */
class WhitneyOneFormGradElementMatrixProvider {
 public:
  /**
   * @brief Constructor
   */
  WhitneyOneFormGradElementMatrixProvider(){};

  /**
   * @brief Compute the element matrix for some cell of a mesh
   * @param entity The mesh cell to compute the element matrix for
   * @returns The 3 by 3 element matrix of the cell
   *
   * @note This function only works for triangles
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

#endif  // THESIS_ASSEMBLE_WHITNEY_ONE_FORM_GRAD_ELEMENT_MATRIX_PROVIDER_H
