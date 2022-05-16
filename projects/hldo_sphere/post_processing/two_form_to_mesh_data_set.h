#ifndef THESIS_POST_PROCESSING_TWO_FORM_TO_MESH_DATA_SET_H
#define THESIS_POST_PROCESSING_TWO_FORM_TO_MESH_DATA_SET_H

/**
 * @file two_form_to_mesh_data_set.h
 * @brief The piecewise constant resultfunction to
 * cell data
 */

#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/codim_mesh_data_set.h>

#include <Eigen/Dense>

namespace projects::hldo_sphere {

namespace post_processing {

/**
 * @brief Extract the basis function coefficients from the solution vector
 * @param mesh A shared pointer to the corresponding mesh
 * @param solution The solution vector obtained from solving the PDE
 * @returns A mesh data set mapping cells of the mesh to the basis function
 * coefficient associated with them
 */
lf::mesh::utils::CodimMeshDataSet<double> extractSolution(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const Eigen::VectorXd &solution);

}  // end namespace post_processing

}  // namespace projects::hldo_sphere

#endif  // THESIS_POST_PROCESSING_TWO_FORM_TO_MESH_DATA_SET_H
