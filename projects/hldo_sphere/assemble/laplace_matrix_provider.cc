#include "laplace_matrix_provider.h"

#include <lf/uscalfe/lagr_fe.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::assemble {

Eigen::MatrixXd LaplaceMatrixProvider::Eval(const lf::mesh::Entity &entity) {
  // Only triangles are supported
  LF_VERIFY_MSG(entity.RefEl() == lf::base::RefEl::kTria(),
                "Unsupported cell type " << entity.RefEl());

  // Get the geometry of the entity
  const auto *geom = entity.Geometry();

  // Compute the global vertex coordinates
  const Eigen::MatrixXd vertices = geom->Global(entity.RefEl().NodeCoords());

  // Construct the basis functions from curl of the standard hat functions
  const lf::uscalfe::FeLagrangeO1Tria<double> hat_func;

  // The gradients are constant on the triangle
  const Eigen::MatrixXd ref_grads =
      hat_func.GradientsReferenceShapeFunctions(Eigen::VectorXd::Zero(2))
          .transpose();

  // The JacobianInverseGramian is constant on the triangle
  const Eigen::MatrixXd J_inv_trans =
      geom->JacobianInverseGramian(Eigen::VectorXd::Zero(2));

  // Compute gradients
  Eigen::MatrixXd grads = J_inv_trans * ref_grads;

  // Compute the element matrix for the product of baricentric functions
  Eigen::MatrixXd elem_mat(3, 3);

  elem_mat = grads.transpose() * grads;

  // correct for integral
  elem_mat *= lf::geometry::Volume(*geom);

  return elem_mat;
}

}  // namespace projects::hldo_sphere::assemble
