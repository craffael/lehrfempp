#ifndef LF_LD_TWONORM_H
#define LF_LD_TWONORM_H
/**
 * @file
 * @brief Conmputation of the L2-norm of a piecewise linear finite element
 * function in different ways
 * @author Ralf Hiptmair
 * @date   April 2019
 * @copyright MIT License
 */

#include <cmath>
#include "lecturedemo.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"
#include "lf/uscalfe/uscalfe.h"

namespace lecturedemo {

/** @brief Computation of L2-norm of a piecewise linear finite element function
 *          on a triangular mesh via the Galerkin mass matrix
 *   @parm dofh \ref DofHandler object providing the mesh and local-to-global
 *         indexing
 *   @param uvec basis expansion coefficient vector for FE function
 *   @return L2-norm of FE function
 */
double l2normByMassMatrix(const lf::assemble::DofHandler &dofh,
                          const Eigen::VectorXd &uvec);

/** @brief Computation of L2-norm of a piecewise linear finite element function
 *          on a triangular mesh via direct edge based quadrature.
 *   @parm dofh \ref DofHandler object providing the mesh and local-to-global
 *         indexing
 *   @param uvec basis expansion coefficient vector for FE function
 *   @return L2-norm of FE function
 */
double l2normByQuadrature(const lf::assemble::DofHandler &dofh,
                          const Eigen::VectorXd &uvec);

/** @brief Computation of L2-norm of a piecewise linear finite element function
 *          on a triangular mesh via direct edge based quadrature.
 *  @parm dofh \ref DofHandler object providing the mesh and local-to-global
 *        indexing
 *  @param uvec basis expansion coefficient vector for FE function
 *  @return L2-norm of FE function
 */
double l2normByMeshFunction(
    const std::shared_ptr<lf::uscalfe::UniformScalarFESpace<double>> &fe_space,
    const Eigen::VectorXd &uvec);

/**
 * @brief Driver routine for demos for LehrFEM++ matrix/vector
 * assembly functions
 */
void lecturedemotwonorm();

}  // namespace lecturedemo

#endif
