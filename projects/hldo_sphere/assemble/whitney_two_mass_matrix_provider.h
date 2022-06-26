#ifndef THESIS_ASSEMBLE_WHITNEY_TWO_MASS_MATRIX_PROVIDER_H
#define THESIS_ASSEMBLE_WHITNEY_TWO_MASS_MATRIX_PROVIDER_H

/**
 * @file whitney_two_mass_matrix_provider.h
 */

#include <lf/mesh/entity.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace assemble {

/**
 * @brief Element matrix provider for the bilinear form with cellwise constant
 * basis functions
 *
 * @f[
 * (u, v)  \mapsto \int\limits_{K} u \, v \ dx
 * @f]
 *
 * The element matrix provider works in a 3 dimensional world with 2
 * dimensional triangular cells.
 *
 * Basis functions are the piecewise constant functions
 * Hence the retuned matrix is of size 1x1 and contains the area of the cell
 *
 * @note This class complies with the type requirements for the template
 * argument ENTITY_MATRIX_PROVIDER of the function
 * lf::assemble::AssembleMatrixLocally().
 *
 * @note Only triangular meshes are supported
 *
 */
class WhitneyTwoMassMatrixProvider {
 public:
  /**
   * @brief Constructor
   */
  WhitneyTwoMassMatrixProvider(){};

  /**
   * @brief Compute the element matrix for some cell of a mesh
   * @param entity The mesh cell to compute the element matrix for
   * @returns The 1 by 1 element matrix of the cell
   *
   * @note Only triangluar meshes are supported
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

#endif  // THESIS_ASSEMBLE_WHITNEY_TWO_MASS_MATRIX_PROVIDER_H
