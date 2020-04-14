#ifndef THESIS_POST_PROCESSING_SOLUTION_TO_MESH_DATA_SET_H
#define THESIS_POST_PROCESSING_SOLUTION_TO_MESH_DATA_SET_H

/**
 * @file solution_to_mesh_data_set.h
 * @brief Mappings from solution vectors of the LSE to mesh data sets
 * representing physical properties
 */

#include <lf/assemble/dofhandler.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>

#include <Eigen/Dense>

namespace projects::ipdg_stokes {

namespace post_processing {

/**
 * @brief Extract the basis function coefficients from the solution vector
 * @param mesh A shared pointer to the corresponding mesh
 * @param dofh The DOF handler for the used FE space
 * @param solution The solution vector obtained from solving the PDE
 * @returns A mesh data set mapping nodes of the mesh to the basis function
 * coefficient associated with them
 */
lf::mesh::utils::CodimMeshDataSet<double> extractBasisFunctionCoefficients(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh, const Eigen::VectorXd &solution);

/**
 * @brief Extract the flow velocity on cells from the solution vector
 * @param mesh A shared pointer to the corresponding mesh
 * @param dofh The DOF handler for the used FE space
 * @param solution The solution vector obtained from solving the PDE
 * @returns A mesh data set mapping cells of the mesh to the flow velocity on
 * them
 */
lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> extractVelocity(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh, const Eigen::VectorXd &solution);

}  // end namespace post_processing

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_POST_PROCESSING_SOLUTION_TO_MESH_DATA_SET_H
