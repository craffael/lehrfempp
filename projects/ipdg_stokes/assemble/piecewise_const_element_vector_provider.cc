#include "piecewise_const_element_vector_provider.h"

#include <lf/base/base.h>
#include <lf/uscalfe/lagr_fe.h>

#include <iostream>

#include <utils.h>

namespace projects::ipdg_stokes::assemble {

PiecewiseConstElementVectorProvider::PiecewiseConstElementVectorProvider(
    double sigma, std::function<Eigen::Vector2d(const Eigen::Vector2d &)> f,
    lf::quad::QuadRule quadrule,
    const lf::mesh::utils::MeshDataSet<bool> &boundary,
    const lf::mesh::utils::MeshDataSet<Eigen::Vector2d> &dirichlet_data)
    : sigma_(sigma),
      f_(std::move(f)),
      quadrule_(std::move(quadrule)),
      boundary_(boundary),
      dirichlet_data_(dirichlet_data) {}

Eigen::VectorXd PiecewiseConstElementVectorProvider::Eval(
    const lf::mesh::Entity &entity) const {
  const auto geom = entity.Geometry();
  // Construct the basis functions from the curl of the standard hat functions
  const lf::uscalfe::FeLagrangeO1Tria<double> hat_func;
  const Eigen::MatrixXd ref_grads =
      hat_func.GradientsReferenceShapeFunctions(Eigen::VectorXd::Zero(2))
          .transpose();
  const Eigen::MatrixXd J_inv_trans =
      geom->JacobianInverseGramian(Eigen::VectorXd::Zero(2));
  const Eigen::MatrixXd grads = J_inv_trans * ref_grads;
  Eigen::MatrixXd basis_funct(2, 3);
  basis_funct << grads.row(1), -grads.row(0);
  // Compute the elements of the vector via the given quadrature rule
  Eigen::VectorXd element_vector = Eigen::VectorXd::Zero(6);
  const Eigen::MatrixXd points = geom->Global(quadrule_.Points());
  const Eigen::VectorXd weights =
      (geom->IntegrationElement(quadrule_.Points()).array() *
       quadrule_.Weights().array())
          .matrix();
  for (lf::base::size_type n = 0; n < quadrule_.NumPoints(); ++n) {
    const auto f_eval = f_(points.col(n));
    element_vector.head(3) += basis_funct.transpose() * weights[n] * f_eval;
  }
  // Weakly impose the dirichlet boundary conditions
  Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());
  // Use the vertex coordinates to compute the local normals on the edges
  const Eigen::Matrix<double, 2, 3> normals =
      projects::ipdg_stokes::mesh::computeOutwardNormals(entity);
  // Iterate over all edges of the current entity and set the boundary
  // conditions where necessary
  const auto edges = entity.SubEntities(1);
  for (int edge_idx = 0; edge_idx < 3; ++edge_idx) {
    const auto &edge = *edges[edge_idx];
    if (boundary_(edge)) {
      Eigen::Vector2d t;
      t << -normals(1, edge_idx), normals(0, edge_idx);
      element_vector[edge_idx + 3] = dirichlet_data_(edge).transpose() * t;
    }
  }
  return element_vector;
}

}  // end namespace projects::ipdg_stokes::assemble
