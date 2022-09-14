#ifndef HLDO_SPHERE_ASSEMBLE_MASS_MATRIX_PROVIDER_H
#define HLDO_SPHERE_ASSEMBLE_MASS_MATRIX_PROVIDER_H

/**
 * @file mass_matrix_provider.h
 */

#include <lf/mesh/entity.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace assemble {

/**
 * @brief Element Matrix Provider for the mass matrix using picewise linear
 * barycentric basis functions.
 *
 * The element matrix provider works in a 3 dimensional world with 2
 * dimensional triangular cells.
 *
 * The bilinear form is locally evaluated for the basis functions.
 * The locally evaluated bilinear form is
 * @f[
 * (u,v) \mapsto \int\limits_{K} u \, v\, dx
 * @f]
 *
 * Details regarding the mathematical derivations can be found in the thesis
 * `Hodge-Laplacians and Dirac Operators on the Surface of the 3-Sphere`
 * section 4.2.2.
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * @note Only triangular meshes are supported
 *
 */
class MassMatrixProvider {
 public:
  /**
   * @brief Constructor
   */
  MassMatrixProvider(){};

  /**
   * @brief Compute the element matrix for a given cell of a mesh
   * @param entity The mesh cell on which the element matrix will be
   * computed
   * @returns The 3 by 3 element matrix of the cell
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

#endif  // HLDO_SPHERE_ASSEMBLE_LAPLACE_MATRIX_PROVIDER_H
