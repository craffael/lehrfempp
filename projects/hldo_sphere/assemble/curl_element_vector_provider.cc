#include "curl_element_vector_provider.h"

#include <lf/base/base.h>
#include <lf/uscalfe/lagr_fe.h>

#include <iostream>

namespace projects::hldo_sphere::assemble {

CurlElementVectorProvider::CurlElementVectorProvider(
    std::function<Eigen::Vector3d(const Eigen::Vector3d &)> f,
    lf::quad::QuadRule quadrule)
    : f_(std::move(f)), quadrule_(std::move(quadrule)) {}

Eigen::VectorXd CurlElementVectorProvider::Eval(
    const lf::mesh::Entity &entity) const {
  const auto *const geom = entity.Geometry();

  // Compute the global vertex coordinates
  Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());

  double eps = 1e-20;
  LF_ASSERT_MSG(vertices.col(0).norm() - vertices.col(1).norm() < eps &&
                    vertices.col(0).norm() - vertices.col(2).norm() < eps,
                "The norms of the vertices have to be equal");

  double r = vertices.col(0).norm();

  // Construct the whitney one form basis functions
  const lf::uscalfe::FeLagrangeO1Tria<double> hat_func;
  const Eigen::MatrixXd ref_grads =
      hat_func.GradientsReferenceShapeFunctions(Eigen::VectorXd::Zero(2))
          .transpose();
  const Eigen::MatrixXd J_inv_trans =
      geom->JacobianInverseGramian(Eigen::VectorXd::Zero(2));
  const Eigen::MatrixXd grads = J_inv_trans * ref_grads;

  // get edge orientations
  auto edgeOrientations = entity.RelativeOrientations();
  Eigen::Vector3d s;
  s(0) = lf::mesh::to_sign(edgeOrientations[0]);
  s(1) = lf::mesh::to_sign(edgeOrientations[1]);
  s(2) = lf::mesh::to_sign(edgeOrientations[2]);

  const auto b_hat = [&](Eigen::Vector2d x_hat) -> Eigen::MatrixXd {
    Eigen::MatrixXd b_hat(3, 3);
    Eigen::VectorXd lambda_hat = hat_func.EvalReferenceShapeFunctions(x_hat);
    for (int i = 0; i < 3; i++) {
      int ip1 = (i + 1) % 3;
      b_hat.col(i) = s(i) * (lambda_hat(i) * grads.col(ip1) -
                             lambda_hat(ip1) * grads.col(i));
    }
    return b_hat;
  };

  const auto f_tilde_hat = [&](Eigen::Vector2d x_hat) -> Eigen::Vector3d {
    Eigen::Vector3d x = geom->Global(x_hat);
    const Eigen::Vector3d x_scaled = x / x.norm() * r;
    Eigen::Vector3d result = b_hat(x_hat).transpose() * f_(x_scaled);
    return result;
  };

  // Compute the elements of the vector with given quadrature rule
  Eigen::Vector3d element_vector;
  const Eigen::MatrixXd points = quadrule_.Points();
  const Eigen::VectorXd weights =
      (geom->IntegrationElement(points).array() * quadrule_.Weights().array())
          .matrix();
  for (lf::base::size_type n = 0; n < quadrule_.NumPoints(); ++n) {
    const Eigen::Vector3d f_eval = f_tilde_hat(points.col(n));
    element_vector += weights[n] * f_eval;
  }
  return element_vector;
}

}  // namespace projects::hldo_sphere::assemble
