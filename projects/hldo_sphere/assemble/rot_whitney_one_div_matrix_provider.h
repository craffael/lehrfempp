#ifndef HLDO_SPHERE_ASSEMBLE_ROT_WHITNEY_ONE_DIV_MATRIX_PROVIDER_H
#define HLDO_SPHERE_ASSEMBLE_ROT_WHITNEY_ONE_DIV_MATRIX_PROVIDER_H

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
 * (\mathbf{v}, u) \mapsto \int\limits_K \text{div}_{\Gamma}(\mathbf{v})\, u
 * \ dx
 * @f]
 *
 * The element matrix provider works in a 3 dimensional world with 2
 * dimensional triangular cells.
 *
 * As basis functions, the rotated  Whitney 1-forms, surface edge
 * elements are used for @f$\mathbf{v}@f$
 * and the cellwise constant basis functions are used for u
 *
 * The whitney 1-forms, surface edge elements are associated with
 * edges and defined as
 *
 * @f[
 *  \mathbf{b}_i = s_i (\lambda_i \mathbf{grad}_{\Gamma}(\lambda_{i+1}) -
 * \lambda_{i+1} \mathbf{grad}_{\Gamma}(\lambda_{i}))
 * @f]
 *
 * with @f$ \lambda_i @f$ barycentric basis functions and @f$ s_i @f$ is a sign
 * of the function based on the relative orientation of the edge in the mesh.
 *
 * The element matrices are equivalent to the element Matrices for
 * @f[
 * (\mathbf{v}', u) \mapsto \int\limits_K \text{curl}_{\Gamma}(\mathbf{v}')\, u
 * \ dx
 * @f]
 * with the normal whintey one basis functions for @f$\mathbf{v}'@f$
 *
 * Details regarding the mathematical derivations can be found in the thesis
 * `Hodge-Laplacians and Dirac Operators on the Surface of the 3-Sphere`
 * section 4.2.8.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
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
  static Eigen::MatrixXd Eval(const lf::mesh::Entity &entity);

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
};

}  // end namespace assemble

}  // namespace projects::hldo_sphere

#endif  // HLDO_SPHERE_ASSEMBLE_ROT_WHITNEY_ONE_DIV_MATRIX_PROVIDER_H
