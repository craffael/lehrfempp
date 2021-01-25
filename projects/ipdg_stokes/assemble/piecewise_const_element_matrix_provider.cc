#include "piecewise_const_element_matrix_provider.h"

#include <lf/uscalfe/lagr_fe.h>
#include <utils.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::ipdg_stokes::assemble {

PiecewiseConstElementMatrixProvider::PiecewiseConstElementMatrixProvider(
    double sigma, const lf::mesh::utils::MeshDataSet<bool> &boundary,
    bool modified)
    : sigma_(sigma), boundary_(boundary), modified_(modified) {}

Eigen::MatrixXd PiecewiseConstElementMatrixProvider::Eval(
    const lf::mesh::Entity &entity) const {
  // Get the geometry of the entity
  const auto geom = entity.Geometry();
  // Compute the global vertex coordinates
  Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());
  // Use the vertex coordinates to compute the local normals on the edges
  const Eigen::Matrix<double, 2, 3> normals =
      projects::ipdg_stokes::mesh::computeOutwardNormals(entity);
  // Construct the basis functions from the curl of the standard hat functions
  const lf::uscalfe::FeLagrangeO1Tria<double> hat_func;
  const Eigen::MatrixXd ref_grads =
      hat_func.GradientsReferenceShapeFunctions(Eigen::VectorXd::Zero(2))
          .transpose();
  const Eigen::MatrixXd J_inv_trans =
      geom->JacobianInverseGramian(Eigen::VectorXd::Zero(2));
  const Eigen::Matrix<double, 2, 3> grads = J_inv_trans * ref_grads;
  Eigen::Matrix<double, 2, 3> basis_funct;
  basis_funct << grads.row(1), -grads.row(0);
  // Fill the element matrix
  Eigen::MatrixXd elem_mat = Eigen::MatrixXd::Zero(6, 6);
  // Upper left block
  elem_mat.block(0, 0, 3, 3).setZero();
  // Lower left and upper right blocks
  const auto edges = entity.SubEntities(1);
  for (int basis_func_idx = 0; basis_func_idx < 3; ++basis_func_idx) {
    for (int edge_idx = 0; edge_idx < 3; ++edge_idx) {
      Eigen::Vector2d t;
      t << -normals(1, edge_idx), normals(0, edge_idx);
      elem_mat(edge_idx + 3, basis_func_idx) =
          basis_funct.col(basis_func_idx).transpose() * t;
    }
  }
  elem_mat.block(0, 3, 3, 3) = elem_mat.block(3, 0, 3, 3).transpose();
  // Lower right block
  for (int edge_idx = 0; edge_idx < 3; ++edge_idx) {
    const auto &edge = *edges[edge_idx];
    elem_mat(edge_idx + 3, edge_idx + 3) = -1 / sigma_;
    if (modified_) {
      elem_mat(edge_idx + 3, edge_idx + 3) /=
          lf::geometry::Volume(*(edge.Geometry()));
    }
    // If the edge is not on the boundary, the entry in the full system matrix
    // will be the sum of two entries of element matrices but we want it still
    // to be equal to 1/sigma. For this reason, we divide the local contribution
    // by 2
    if (!boundary_(edge)) {
      elem_mat(edge_idx + 3, edge_idx + 3) /= 2;
    }
  }
  return elem_mat;
}

}  // end namespace projects::ipdg_stokes::assemble
