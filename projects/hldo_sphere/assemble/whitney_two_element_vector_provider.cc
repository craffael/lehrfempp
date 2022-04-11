#include "whitney_two_element_vector_provider.h"

#include <lf/base/base.h>
#include <lf/uscalfe/lagr_fe.h>

#include <iostream>

namespace projects::hldo_sphere::assemble {

WhitneyTwoElementVectorProvider::WhitneyTwoElementVectorProvider(
    std::function<double(const Eigen::Vector3d &)> f,
    lf::quad::QuadRule quadrule)
    : f_(std::move(f)), quadrule_(std::move(quadrule)) {}

Eigen::VectorXd WhitneyTwoElementVectorProvider::Eval(
    const lf::mesh::Entity &entity) const {
  const auto *const geom = entity.Geometry();

  // Compute the global vertex coordinates
  Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());

  double eps = 1e-20;
  LF_ASSERT_MSG(vertices.col(0).norm() - vertices.col(1).norm() < eps &&
                    vertices.col(0).norm() - vertices.col(2).norm() < eps,
                "The norms of the vertices have to be equal");

  double r = vertices.col(0).norm();

  // Compute the elements of the vector with given quadrature rule
  const Eigen::MatrixXd points = geom->Global(quadrule_.Points());
  const Eigen::VectorXd weights =
      (geom->IntegrationElement(quadrule_.Points()).array() *
       quadrule_.Weights().array())
          .matrix();

  double sum = 0;
  for (lf::base::size_type n = 0; n < quadrule_.NumPoints(); ++n) {
    Eigen::VectorXd x = geom->Global(points.col(n));
    sum += weights[n] * f_(x / x.norm() * r);
  }

  Eigen::VectorXd element_vector(1);
  element_vector(0) = sum;

  return element_vector;
}

}  // namespace projects::hldo_sphere::assemble
