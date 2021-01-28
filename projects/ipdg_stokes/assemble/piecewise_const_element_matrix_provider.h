#ifndef THESIS_ASSEMBLE_PIECEWISE_CONST_ELEMENT_MATRIX_PROVIDER_H
#define THESIS_ASSEMBLE_PIECEWISE_CONST_ELEMENT_MATRIX_PROVIDER_H

/**
 * @file piecewise_const_element_matrix_provider.h
 * @brief An element matrix provider for vector valued divergence free piecewise
 * constant basis functions
 */

#include <lf/mesh/utils/mesh_data_set.h>

#include <Eigen/Dense>

namespace projects::ipdg_stokes {

namespace assemble {

/**
 * @brief Element matrix provider for the stokes system
 *
 * As basis functions, the curl of the standard linear Lagrangian FEM hat
 * functions is used in order to guarantee divergence-freeness of the solution.
 *
 * @note Currently, only triangular meshes are supported
 *
 * The (global) bilinear form is given by
 * @f[
 *  (u, v) \mapsto \sum_{e \in \mathcal{E}(\mathcal{M})} \int_e\!
 * \frac{\sigma}{|e|}\left[\!\left[ u \cdot \tau_e \right]\!\right] :
 * \left[\!\left[ v \cdot \tau_e \right]\!\right] \,\mathrm{d}x
 * @f]
 * To enable local assembly of the system matrix, auxilary dofs are introduced
 * on the edges of the mesh which are then forced to be equal to the jump of
 * @f$u@f$ over that edge. As we use Raviart-Thomas basis functions, the normal
 * component of the function is continuous over the mesh edges. For this reason,
 * the jump terms only consist of the jump in the tangential direction over an
 * edge.
 */
class PiecewiseConstElementMatrixProvider {
 public:
  /**
   * @brief Constructor
   * @param sigma The stabilization coefficient of the IPDG method
   * @param boundary A MeshDataSet marking the boundary edges
   * @param modified If true, the modified penalty term is used missing the
   * scaling with 1/|e_n|
   */
  PiecewiseConstElementMatrixProvider(
      double sigma, const lf::mesh::utils::MeshDataSet<bool> &boundary,
      bool modified = false);

  /**
   * @brief Compute the element matrix for some entity of a mesh
   * @param entity The mesh entity to compute the element matrix for
   * @returns The 6 by 6 element matrix of the given entity
   *
   * @note This function only works for triangles
   */
  Eigen::MatrixXd Eval(const lf::mesh::Entity &entity) const;

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
  const double sigma_;
  const lf::mesh::utils::MeshDataSet<bool> &boundary_;
  const bool modified_;
};

}  // end namespace assemble

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_ASSEMBLE_PIECEWISE_CONST_ELEMENT_MATRIX_PROVIDER_H
