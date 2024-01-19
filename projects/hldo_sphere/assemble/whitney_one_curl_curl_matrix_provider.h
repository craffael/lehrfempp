#ifndef HLDO_SPHERE_ASSEMBLE_WHITNEY_ONE_CURL_CURL_MATRIX_PROVIDER_H
#define HLDO_SPHERE_ASSEMBLE_WHITNEY_ONE_CURL_CURL_MATRIX_PROVIDER_H

/**
 * @file whitney_one_form_curl_element_matrix_provider.h
 */

#include <lf/mesh/entity.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace assemble {

/**
 * @brief An elment matrix provider computed with the Whitney one forms, vector
 * valued basis function
 *
 * The element matrix provider works in a 3 dimensional world with 2
 * dimensional triangular cells.
 *
 * @f[
 * (\mathbf{u}, \mathbf{v})  \mapsto \int\limits_K
 * \text{curl}_{\Gamma}(\mathbf{u}) \ \text{curl}_{\Gamma}(\mathbf{v}) dx
 * @f]
 *
 * The whitney 1-forms, surface edge elements are associated with
 * edges and defined as
 *
 * @f[
 *  \mathbf{b}_i = s_i (\lambda_i \mathbf{grad}_{\Gamma}(\lambda_{i+1}) -
 * \lambda_{i+1} \mathbf{grad}_{\Gamma}(\lambda_{i}))
 * @f]
 *
 * Details regarding the mathematical derivations can be found in the thesis
 * `Hodge-Laplacians and Dirac Operators on the Surface of the 3-Sphere`
 * section 4.2.4.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * @note Only triangular meshes are supported
 *
 */
class WhitneyOneCurlCurlMatrixProvider {
 public:
  /**
   * @brief Constructor
   */
  WhitneyOneCurlCurlMatrixProvider(){};

  /**
   * @brief Compute the element matrix for some triangle of a mesh
   * @param entity The mesh triangle to compute the element matrix for
   * @returns The 3 by 3 element matrix of the triangle
   *
   * @note Currently, only triangular cells are supported
   */
  static Eigen::MatrixXd Eval(const lf::mesh::Entity &entity);

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
};

}  // end namespace assemble

}  // namespace projects::hldo_sphere

#endif  // HLDO_SPHERE_ASSEMBLE_WHITNEY_ONE_CURL_CURL_MATRIX_PROVIDER_H
