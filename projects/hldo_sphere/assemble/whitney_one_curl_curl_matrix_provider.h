#ifndef THESIS_ASSEMBLE_WHITNEY_ONE_CURL_CURL_MATRIX_PROVIDER_H
#define THESIS_ASSEMBLE_WHITNEY_ONE_CURL_CURL_MATRIX_PROVIDER_H

/**
 * @file whitney_one_form_curl_element_matrix_provider.h
 */

#include <lf/mesh/entity.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace assemble {

/**
 * @brief Element matrix provider for the bilinear form
 *
 * @f[
 * \int curl_{\Gamma}(u) curl_{\Gamma}(v) dx
 * @f]
 *
 * As basis functions, the Whitney 1-forms, surface edge elements are used
 *
 * @note Currently, only triangular meshes are supported
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
  Eigen::MatrixXd Eval(const lf::mesh::Entity &entity) const;

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
};

}  // end namespace assemble

}  // namespace projects::hldo_sphere

#endif  // THESIS_ASSEMBLE_WHITNEY_ONE_CURL_CURL_MATRIX_PROVIDER_H
