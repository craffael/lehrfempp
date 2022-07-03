#ifndef HLDO_SPHERE_ASSEMBLE_WHITNEY_ONE_GRAD_MATRIX_PROVIDER_H
#define HLDO_SPHERE_ASSEMBLE_WHITNEY_ONE_GRAD_MATRIX_PROVIDER_H

/**
 * @file whitney_one_grad_matrix_provider.h
 */

#include <lf/mesh/entity.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace assemble {

/**
 * @brief Element matrix provider Whitney one forms, surface vector fields
 *
 * The element matrix provider works in a 3 dimensional world with 2
 * dimensional triangular cells. And evaluates the bilinear form
 * @f[
 * (u, \mathbf{v}) \mapsto \int\limits_K \mathbf{grad}_{\Gamma}(u) \, \mathbf{v}
 * \ dx
 * @f]
 *
 * Basis functions are the Whitney 1-forms, surface edge elements for
 * @f$\mathbf{v}@f$ and the barycentric basis functions for @f$u@f$
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
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
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

#endif  // HLDO_SPHERE_ASSEMBLE_WHITNEY_ONE_GRAD_MATRIX_PROVIDER_H
