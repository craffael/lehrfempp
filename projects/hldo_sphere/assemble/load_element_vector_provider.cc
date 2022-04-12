#include "load_element_vector_provider.h"

#include <lf/base/base.h>
#include <lf/quad/quad.h>
#include <lf/quad/quad_rule.h>
#include <lf/uscalfe/lagr_fe.h>

#include <iostream>

namespace projects::hldo_sphere::assemble {

LoadElementVectorProvider::LoadElementVectorProvider(
    std::function<double(const Eigen::Vector3d &)> f)
    : f_(f) {}

Eigen::VectorXd LoadElementVectorProvider::Eval(
    const lf::mesh::Entity &entity) const {
  const auto *const geom = entity.Geometry();

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

  const lf::uscalfe::FeLagrangeO1Tria<double> hat_func;
  const auto lambda_hat = [&](Eigen::Vector2d x_hat) -> Eigen::MatrixXd {
    return hat_func.EvalReferenceShapeFunctions(x_hat);
  };

  const auto f_tilde_hat = [&](Eigen::Vector2d x_hat) -> Eigen::MatrixXd {
    Eigen::Vector3d x = geom->Global(x_hat);
    const Eigen::Vector3d x_scaled = x / x.norm() * r;
    return lambda_hat(x_hat) * f_(x_scaled);
  };

  // Compute the elements of the vector with given quadrature rule
  Eigen::Vector3d element_vector = Eigen::Vector3d::Zero();
  const Eigen::MatrixXd points = quadrule.Points();
  const Eigen::VectorXd weights =
      (geom->IntegrationElement(points).array() * quadrule.Weights().array())
          .matrix();
  for (lf::base::size_type n = 0; n < quadrule.NumPoints(); ++n) {
    const Eigen::Vector3d f_eval = f_tilde_hat(points.col(n));
    element_vector += weights[n] * f_eval;
  }
  return element_vector;
}

}  // namespace projects::hldo_sphere::assemble
