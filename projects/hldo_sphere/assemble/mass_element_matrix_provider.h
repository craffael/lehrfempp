#ifndef THESIS_ASSEMBLE_MASS_ELEMENT_MATRIX_PROVIDER_H
#define THESIS_ASSEMBLE_MASS_ELEMENT_MATRIX_PROVIDER_H

/**
 * @file mass_element_matrix_provider.h
 * @brief An element matrix provider for vector valued piecewise linear basis
 * functions and the bilinear form
 * @f[
 * \int_{\Omega} u v dx
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
 * \int_{\Omega} u v dx
 * @f]
 *
 * Basis functions are the baricentric basis functions
 *
 * @note Only triangular meshes are supported
 *
 */
class MassElementMatrixProvider {
 public:
  /**
   * @brief Constructor
   */
  MassElementMatrixProvider(){};

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

#endif  // THESIS_ASSEMBLE_LAPLACE_ELEMENT_MATRIX_PROVIDER_H
