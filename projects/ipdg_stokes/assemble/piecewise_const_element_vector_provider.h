#ifndef THESIS_ASSEMBLE_PIECEWISE_CONST_ELEMENT_VECTOR_PROVIDER_H
#define THESIS_ASSEMBLE_PIECEWISE_CONST_ELEMENT_VECTOR_PROVIDER_H

/**
 * @file piecewise_const_element_vector_provider,h
 * @brief Assembly routines for the right hand side of the discretized Stokes
 * system
 */

#include <lf/mesh/entity.h>
#include <lf/mesh/utils/mesh_data_set.h>
#include <lf/quad/quad_rule.h>

#include <Eigen/Dense>
#include <functional>
#include <utility>

namespace projects::ipdg_stokes {

namespace assemble {

/**
 * @brief Element vector provider for the stokes system
 *
 * @note Currently, only triangular meshes are supported
 *
 * The linear form is given by
 * @f[
 *  v \mapsto \int_K\! f \cdot v \,\mathrm{d}x +
 * \sum_{e\in\mathcal{E}^\partial(\mathcal{M})} \int_e\! \frac{\sigma}{|e|}
 * \left(g_D\cdot\tau_e\right) \left[\!\left[v\cdot\tau_e\right]\!\right]
 * \,\mathrm{d}x
 * @f]
 */
class PiecewiseConstElementVectorProvider {
 public:
  /**
   * @brief Constructor
   * @param sigma The stabilization parameter for IP-DG
   * @param f A function mapping coordinates to the volumetric forces
   * @param quadrule The quadrature rule used to compute the integral in the
   * linear form
   * @param boundary A mesh data set marking the boundary edges
   * @param dirichlet_data A mesh data set containing the dirichlet data for the
   * velocity on the boundary edges
   */
  PiecewiseConstElementVectorProvider(
      double sigma, std::function<Eigen::Vector2d(const Eigen::Vector2d &)> f,
      lf::quad::QuadRule quadrule,
      const lf::mesh::utils::MeshDataSet<bool> &boundary,
      const lf::mesh::utils::MeshDataSet<Eigen::Vector2d> &dirichlet_data);

  /**
   * @brief Compute the element vector for some entity of the mesh
   * @param entity The mesh entity to compute the element vector for
   * @returns The element vector of the given entity
   *
   * @note This function currently only works for triangles
   */
  Eigen::VectorXd Eval(const lf::mesh::Entity &entity) const;

  /**
   * @brief All entities are regarded as active
   */
  bool isActive(const lf::mesh::Entity &entity) const { return true; }

 private:
  const double sigma_;
  const std::function<Eigen::Vector2d(const Eigen::Vector2d &)> f_;
  const lf::quad::QuadRule quadrule_;
  const lf::mesh::utils::MeshDataSet<bool> &boundary_;
  const lf::mesh::utils::MeshDataSet<Eigen::Vector2d> &dirichlet_data_;
};

}  // end namespace assemble

}  // end namespace projects::ipdg_stokes

#endif  // THESIS_ASSEMBLE_PIECEWISE_CONST_ELEMENT_VECTOR_PROVIDER_H
