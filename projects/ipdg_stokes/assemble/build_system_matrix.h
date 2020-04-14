#ifndef THESIS_ASSEMBLE_BUILD_SYSTEM_MATRIX_H
#define THESIS_ASSEMBLE_BUILD_SYSTEM_MATRIX_H

/**
 * @file build_system_matrix.h
 * @brief Functions to build the system matrix for the stokes system with given
 * Dirichlet boundary conditions
 */

#include <Eigen/Dense>
#include <functional>
#include <tuple>
#include <utility>

#include <lf/assemble/coomatrix.h>
#include <lf/assemble/dofhandler.h>
#include <lf/mesh/entity.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/mesh_data_set.h>
#include <lf/quad/quad_rule.h>

namespace projects::ipdg_stokes {

namespace assemble {

/**
 * @brief Build the system matrix for the stokes system with no flow boundary
 * conditions
 *
 * @param mesh A shared pointer to some mesh instance
 * @param dofh The dof handler for the given mesh
 * @param f A functor mapping 2D-coordinates to volumetric forces
 * @param dirichlet_data A functor mapping entities on the mesh boundary to
 * Dirichlet data for the velocity field
 * @param sigma The stabilization factor
 * @param quadrule The quadrature rule used to integrate the volumetric forces
 * over an element
 * @param modified_penalty If true, uses a modified penalty term where the jumps
 * are not scaled by 1/|e_n|
 * @return A tuple with the first value containing the assembles COOMatrix of
 * the stokes system and the second value being equal to the right hand side of
 * the LSE
 */
std::tuple<lf::assemble::COOMatrix<double>, Eigen::VectorXd>
buildSystemMatrixNoFlow(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh,
    const std::function<Eigen::Vector2d(const Eigen::Vector2d &)> &f,
    const std::function<Eigen::Vector2d(const lf::mesh::Entity &)>
        &dirichlet_data,
    double sigma, const lf::quad::QuadRule &quadrule,
    bool modified_penalty = false);

/**
 * @brief Build the system matrix for the stokes system with in- and out flow
 * boundary conditions
 *
 * @param mesh A shared pointer to some mesh instance
 * @param dofh The dof handler for the given mesh
 * @param f A functor mapping 2D-coordinates to volumetric forces
 * @param dirichlet_data A functor mapping entities on the mesh boundary to
 * Dirichlet data for the velocity field
 * @param sigma The stabilization factor
 * @param quadrule The quadrature rule used to integrate the volumetric forces
 * over an element
 * @param modified_penalty If true, uses a modified penalty term where the jumps
 * are not scaled by 1/|e_n|
 * @return A tuple with the first value containing the assembles COOMatrix of
 * the stokes system, the second value being equal to the right hand side of the
 * LSE and the third value contains the offset function
 */
std::tuple<lf::assemble::COOMatrix<double>, Eigen::VectorXd, Eigen::VectorXd>
buildSystemMatrixInOutFlow(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh,
    const std::function<Eigen::Vector2d(const Eigen::Vector2d &)> &f,
    const std::function<Eigen::Vector2d(const lf::mesh::Entity &)>
        &dirichlet_data,
    double sigma, const lf::quad::QuadRule &quadrule,
    bool modified_penalty = false);

}  // end namespace assemble

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_ASSEMBLE_BUILD_SYSTEM_MATRIX_H
