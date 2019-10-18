#include "norms.h"

#include <iostream>
#include <vector>

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <lf/quad/quad_rule.h>

#include <utils.h>

namespace projects::ipdg_stokes::post_processing {

double L2norm(const std::shared_ptr<const lf::mesh::Mesh> &mesh,
              const std::function<Eigen::Vector2d(const lf::mesh::Entity &,
                                                  const Eigen::Vector2d &)> &f,
              unsigned quadrule_order) {
  const lf::quad::QuadRule quadrule =
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), quadrule_order);
  double norm2 = 0;
  for (const auto entity : mesh->Entities(0)) {
    const auto geom = entity->Geometry();
    const Eigen::MatrixXd points = geom->Global(quadrule.Points());
    const Eigen::VectorXd weights =
        (geom->IntegrationElement(quadrule.Points()).array() *
         quadrule.Weights().array())
            .matrix();
    for (lf::base::size_type i = 0; i < quadrule.NumPoints(); ++i) {
      norm2 += weights[i] * f(*entity, points.col(i)).squaredNorm();
    }
  }
  return std::sqrt(norm2);
}

double DGnorm(const std::shared_ptr<const lf::mesh::Mesh> &mesh,
              const std::function<Eigen::Vector2d(const lf::mesh::Entity &,
                                                  const Eigen::Vector2d &)> &f,
              const std::function<Eigen::Matrix2d(
                  const lf::mesh::Entity &, const Eigen::Vector2d &)> &f_grad,
              unsigned quadrule_order) {
  double norm2 = 0;
  // Compute the DG norm on the elements
  const lf::quad::QuadRule quadrule_triangle =
      lf::quad::make_QuadRule(lf::base::RefEl::kTria(), quadrule_order);
  for (const auto entity : mesh->Entities(0)) {
    const auto geom = entity->Geometry();
    const Eigen::MatrixXd points = geom->Global(quadrule_triangle.Points());
    const Eigen::VectorXd weights =
        (geom->IntegrationElement(quadrule_triangle.Points()).array() *
         quadrule_triangle.Weights().array())
            .matrix();
    for (lf::base::size_type i = 0; i < quadrule_triangle.NumPoints(); ++i) {
      norm2 += weights[i] * f_grad(*entity, points.col(i)).squaredNorm();
    }
  }
  // Compute the DG norm on the edges
  const lf::quad::QuadRule quadrule_segment =
      lf::quad::make_QuadRule(lf::base::RefEl::kSegment(), quadrule_order);
  std::vector<Eigen::Matrix2d> jumps(mesh->NumEntities(1),
                                     Eigen::Matrix2d::Zero());
  // We have to manually set the jump at the edges of the mesh to zero
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1);
  for (const auto entity : mesh->Entities(0)) {
    const auto edges = entity->SubEntities(1);
    const auto normals =
        projects::ipdg_stokes::mesh::computeOutwardNormals(*entity);
    for (int i = 0; i < 3; ++i) {
      const auto &edge = *edges[i];
      if (!boundary(edge)) {
        const auto geom = edge.Geometry();
        const Eigen::MatrixXd points = geom->Global(quadrule_segment.Points());
        const Eigen::VectorXd weights =
            (geom->IntegrationElement(quadrule_segment.Points()).array() *
             quadrule_segment.Weights().array())
                .matrix() /
            lf::geometry::Volume(*geom);
        for (lf::base::size_type p = 0; p < quadrule_segment.NumPoints(); ++p) {
          jumps[mesh->Index(edge)] += weights[p] * f(*entity, points.col(p)) *
                                      normals.col(i).transpose();
        }
      }
    }
  }
  for (const auto &jump : jumps) {
    norm2 += jump.squaredNorm();
  }
  return std::sqrt(norm2);
}

}  // namespace projects::ipdg_stokes::post_processing
