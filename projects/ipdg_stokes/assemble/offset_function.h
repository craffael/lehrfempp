#ifndef THESIS_ASSEMBLE_OFFSET_FUNCTION_H
#define THESIS_ASSEMBLE_OFFSET_FUNCTION_H

/**
 * @file offset_function.h
 * @brief Provides functionality to compute the offset function for in-out flow
 * boundary conditions
 */

#include <functional>

#include <lf/assemble/dofhandler.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/mesh_data_set.h>

namespace projects::ipdg_stokes {

namespace assemble {

/**
 * @brief Local assembler for a matrix mapping boundary basis function
 * coefficients to jumps
 *
 * This assembler builds a matrix mapping the basis function coefficients f
 * nodes on the mesh boundary to the normal component of the jump over the
 * boundary edges. This is used for assembling the offset function, as the
 * dirichlet data is given in terms of boundary jumps, but the actual solution
 * uses the nodal values as basis function coefficients.
 */
class PiecewiseBoundaryNormalJumpAssembler {
 public:
  /**
   * @brief Constructor
   * @param mesh A shared pointer to the mesh on which the PDE is solved
   * @param boundary A MeshDataSet marking the boundary entities
   */
  PiecewiseBoundaryNormalJumpAssembler(
      const std::shared_ptr<const lf::mesh::Mesh> &mesh,
      const lf::mesh::utils::MeshDataSet<bool> &boundary);

  /**
   * @brief Compute the element matrix for some triangle of the mesh
   * @param entity The mesh entity to compute the element matrix for
   * @returns The element matrix for the given entity
   */
  Eigen::MatrixXd Eval(const lf::mesh::Entity &entity) const;

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
  const std::shared_ptr<const lf::mesh::Mesh> mesh_;
  const lf::mesh::utils::MeshDataSet<bool> &boundary_;
};

/**
 * @brief Compute the offset function for given boundary conditions
 * @param mesh A shared pointer the the mesh instance used
 * @param boundary A MeshDataSet marking the boundary elements of the mesh
 * @param dofh The DOF handler for the FE space used
 * @param dirichlet_data A function mapping boundary entities to dirichlet
 * values for the flow velocity
 * @param A The full system matrix of the discretized problem
 * @returns The offset function for the given boundary conditions
 */
Eigen::VectorXd createOffsetFunction(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::mesh::utils::MeshDataSet<bool> &boundary,
    const lf::assemble::DofHandler &dofh,
    const std::function<Eigen::Vector2d(const lf::mesh::Entity &)>
        &dirichlet_data,
    const Eigen::SparseMatrix<double> &A);

}  // end namespace assemble

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_ASSEMBLE_OFFSET_FUNCTION_H
