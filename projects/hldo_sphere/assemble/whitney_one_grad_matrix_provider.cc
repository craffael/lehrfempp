#include "whitney_one_grad_matrix_provider.h"

#include <lf/uscalfe/lagr_fe.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::assemble {

Eigen::MatrixXd WhitneyOneGradMatrixProvider::Eval(
    const lf::mesh::Entity &entity) const {
  // Only triangles are supported
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << entity.RefEl());

  // Get the geometry of the entity
  const auto *geom = entity.Geometry();

  // Compute the global vertex coordinates
  Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());

  // Construct the basis functions from curl of the standard hat functions
  const lf::uscalfe::FeLagrangeO1Tria<double> hat_func;
  // The gradients are constant on the triangle
  const Eigen::MatrixXd ref_grads =
      hat_func.GradientsReferenceShapeFunctions(Eigen::VectorXd::Zero(2))
          .transpose();
  // The JacobianInverseGramian is constant on the triangle
  const Eigen::MatrixXd J_inv_trans =
      geom->JacobianInverseGramian(Eigen::VectorXd::Zero(2));

  // get the gradients
  const Eigen::MatrixXd grad = J_inv_trans * ref_grads;

  // create whitney 1-Form basis functions only the constants
  Eigen::MatrixXd bs(grad.rows(), 3);
  bs.col(0) = grad.col(1) - grad.col(0);
  bs.col(1) = grad.col(2) - grad.col(1);
  bs.col(2) = grad.col(0) - grad.col(2);

  // correct for orientation
  auto edgeOrientations = entity.RelativeOrientations();
  bs.col(0) *= lf::mesh::to_sign(edgeOrientations[0]);
  bs.col(1) *= lf::mesh::to_sign(edgeOrientations[1]);
  bs.col(2) *= lf::mesh::to_sign(edgeOrientations[2]);

  // fill the element matrix
  Eigen::MatrixXd elem_mat(3, 3);
  elem_mat = bs.transpose() * grad;

  // correct for integral
  elem_mat *= lf::geometry::Volume(*geom) / 3.0;

  return elem_mat;
}

}  // namespace projects::hldo_sphere::assemble
