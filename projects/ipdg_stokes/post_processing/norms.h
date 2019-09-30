#ifndef THESIS_POST_PROCESSING_NORMS_H
#define THESIS_POST_PROCESSING_NORMS_H

/**
 * @file norms.h
 * @brief Functions to compute norms of discontinuous vector valued functions
 * over a mesh
 */

#include <lf/mesh/entity.h>
#include <lf/mesh/mesh.h>

#include <functional>

namespace projects::ipdg_stokes {

namespace post_processing {

/**
 * @breif Compute the @f$L^2@f$-norm of a vector valued function over a mesh
 * @param mesh A shared pointer to the mesh over which the function is defined
 * @param f The function to take the norm from
 * @param quadrule_order The order of the quadrature rule to use to integrate
 * over the elements. Defaults to 50
 * @returns The @f$L^2@f$-norm of the function over the mesh
 *
 * This function computes the value of the following formula
 * @f[
 *  ||f||_0^2 = \sum_{K\in\mathcal{M}} \int_K\! ||f(x)||^2 \,\mathrm{d}x
 * @f]
 */
double L2norm(const std::shared_ptr<const lf::mesh::Mesh> &mesh,
              const std::function<Eigen::Vector2d(const lf::mesh::Entity &,
                                                  const Eigen::Vector2d &)> &f,
              unsigned quadrule_order = 50);

/**
 * @brief Compute the DG-norm of a vector valued function over a mesh
 * @param mesh A shared pointer to the mesh over which the function is defined
 * @param f The function to take the norm from
 * @param f_grad The gradient of the function to take the norm from
 * @param quadrule_order The order of the quadrature rule to use to integrate
 * over the elements. Defaults to 50
 * @returns The @f$H^1@f$-norm of the function over the mesh
 *
 * This functoin computes the following formula
 * @f[
 *  ||f||_{1,h}^2 = \sum_{K\in\mathcal{M}} \int_K\! ||\nabla f(x)||^2
 * \,\mathrm{d}x + \sum_{e\in\mathcal{E}(\mathcal{M})} \int_e\!
 * \frac{1}{|e|}||\left[\!\left[ f(x) \otimes n \right]\!\right]||^2
 * \,\mathrm{d}x
 * @f]
 */
double DGnorm(const std::shared_ptr<const lf::mesh::Mesh> &mesh,
              const std::function<Eigen::Vector2d(const lf::mesh::Entity &,
                                                  const Eigen::Vector2d &)> &f,
              const std::function<Eigen::Matrix2d(
                  const lf::mesh::Entity &, const Eigen::Vector2d &)> &f_grad,
              unsigned quadrule_order = 50);

}  // end namespace post_processing

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_POST_PROCESSING_NORMS_H
