#ifndef THESIS_MESH_UTILS_H
#define THESIS_MESH_UTILS_H

/**
 * @file utils.h
 * @brief Various utility functions for meshes not fitting in any specific
 * category
 */

#include <lf/mesh/entity.h>

namespace projects::ipdg_stokes {

namespace mesh {

/**
 * @brief Compute the outward pointing normals of a triangle
 * @param entity The triangle for which the outward pointing normals should be
 * computed
 * @returns A matrix containing the outward pointing normals in its columns
 *
 * The n-th column of the matrix returned corresponds to the outward pointing normal
 * on the n-th edge of the given triangle.
 */
Eigen::Matrix<double, 2, 3> computeOutwardNormals(
    const lf::mesh::Entity &entity);

/**
 * @brief Compute the tangentials onto a triangle
 * @param entity The triangle for which the tangentials should be computed
 * @returns A matrix containing the tangentials in its columns
 *
 * The n-th column of the matrix returned corresponds to the tangential vector
 * on the n-th edge of the given triangle. The tangentials on an edge are
 * guaranteed to point in opposite directions when taken from the left
 * or the right adjacent triangle.
 */
Eigen::Matrix<double, 2, 3> computeTangentials(const lf::mesh::Entity &entity);

}  // end namespace mesh

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_MESH_UTILS_H
