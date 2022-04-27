#include "whitney_two_vector_provider.h"

#include <lf/base/base.h>
#include <lf/quad/quad.h>
#include <lf/quad/quad_rule.h>
#include <lf/uscalfe/lagr_fe.h>

#include <iostream>

namespace projects::hldo_sphere::assemble {

WhitneyTwoVectorProvider::WhitneyTwoVectorProvider(
    std::function<double(const Eigen::Vector3d &)> f)
    : f_(std::move(f)) {}

Eigen::VectorXd WhitneyTwoVectorProvider::Eval(
    const lf::mesh::Entity &entity) const {
  const auto *const geom = entity.Geometry();

  // Only triangles are supported
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << entity.RefEl());

  // Compute the global vertex coordinates
  Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());

  double eps = 1e-20;
  LF_ASSERT_MSG(vertices.col(0).norm() - vertices.col(1).norm() < eps &&
                    vertices.col(0).norm() - vertices.col(2).norm() < eps,
                "The norms of the vertices have to be equal");

  double r = vertices.col(0).norm();

  // define quad rule with sufficiantly high degree since the
  // baricentric coordinate functions have degree 1
  lf::quad::QuadRule quadrule{lf::quad::make_TriaQR_EdgeMidpointRule()};

  // Compute the elements of the vector with given quadrature rule
  const Eigen::MatrixXd points = geom->Global(quadrule.Points());
  const Eigen::VectorXd weights =
      (geom->IntegrationElement(quadrule.Points()).array() *
       quadrule.Weights().array())
          .matrix();

  double sum = 0;
  for (lf::base::size_type n = 0; n < quadrule.NumPoints(); ++n) {
    Eigen::VectorXd x = geom->Global(points.col(n));
    sum += weights[n] * f_(x / x.norm() * r);
  }

  Eigen::VectorXd element_vector(1);
  element_vector(0) = sum;

  return element_vector;
}

}  // namespace projects::hldo_sphere::assemble
