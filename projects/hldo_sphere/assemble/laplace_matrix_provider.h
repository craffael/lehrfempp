#ifndef THESIS_ASSEMBLE_LAPLACE_MATRIX_PROVIDER_H
#define THESIS_ASSEMBLE_LAPLACE_MATRIX_PROVIDER_H

/**
 * @file laplace_matrix_provider.h
 *
 */

#include <lf/mesh/entity.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace assemble {

/**
 * @brief An element matrix provider for piecewise linear (barycentric) basis
 * functions in a 3 dimensional world with 2 dimensional triangular cells
 *
 * The locally evaluated bilinear form is
 * @f[
 * (u,v) \mapsto \int_{K} \mathbf{grad}_{\Gamma}(u) \cdot
 * \mathbf{grad}_{\Gamma}(v) dx
 * @f]
 *
 * Where \f$ \mathbf{grad}_{\Gamma} \f$ denotes the tangential
 * gradient on the triangle.
 *
 *
 * @note Only triangular meshes are supported
 *
 */
class LaplaceMatrixProvider {
 public:
  /**
   * @brief Constructor
   */
  LaplaceMatrixProvider(){};

  /**
   * @brief Compute the element matrix for a given triangle of a mesh
   * @param entity The mesh triangle on which the element matrix will be
   * computed
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

#endif  // THESIS_ASSEMBLE_LAPLACE_MATRIX_PROVIDER_H
