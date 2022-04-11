#include "rot_w_one_form_dot_element_matrix_provider.h"

#include <lf/uscalfe/lagr_fe.h>

#include <cmath>
#include <iomanip>
#include <iostream>

namespace projects::hldo_sphere::assemble {

Eigen::MatrixXd RotWOneFormDotElementMatrixProvider::Eval(
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

  // Compute gradients
  Eigen::MatrixXd grads = J_inv_trans * ref_grads;

  // correct for orientation
  auto edgeOrientations = entity.RelativeOrientations();
  Eigen::VectorXd s(3);
  s << lf::mesh::to_sign(edgeOrientations[0]),
      lf::mesh::to_sign(edgeOrientations[1]),
      lf::mesh::to_sign(edgeOrientations[2]);

  // Compute the element matrix for the product of baricentric functions
  std::vector<Eigen::Matrix3d> elem_mats(4);
  // clang-format off
  elem_mats[0] <<   2, 1, 1,
                    1, 2, 1,
                    1, 1, 2;
  elem_mats[1] <<   1, 1, 2,
                    2, 1, 1,
                    1, 2, 1;
  elem_mats[1] *= -1;
  elem_mats[2] <<   1, 2, 1,
                    1, 1, 2,
                    2, 1, 1;
  elem_mats[2] *= -1;
  elem_mats[3] <<   2, 1, 1,
                    1, 2, 1,
                    1, 1, 2;
  // clang-format on

  for (int i = 0; i < 4; i++) {
    elem_mats[i] *= lf::geometry::Volume(*geom) / 12.;
  }

  // index modulo three
  auto index = [](int i) -> int { return i % 3; };

  // compute the contributions of gradients for the first term
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      elem_mats[0](i, k) *=
          (grads.col(index(i + 1)).transpose() * grads.col(index(k + 1)));
      std::cout << "\ngrad(" << index(i + 1) << ") times "
                << "\ngrad(" << index(k + 1) << ") equals "
                << (grads.col(index(i + 1)).transpose() *
                    grads.col(index(k + 1)))
                << "\n\n\n";
      elem_mats[0](i, k) *= s(i) * s(k);
    }
  }

  // compute the contributions of gradients for the second term
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      elem_mats[1](i, k) *=
          (grads.col(index(i + 1)).transpose() * grads.col(index(k)));
      elem_mats[1](i, k) *= s(i) * s(k);
    }
  }

  // compute the contributions of gradients for the third term
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      elem_mats[2](i, k) *=
          (grads.col(index(i)).transpose() * grads.col(index(k + 1)));
      elem_mats[2](i, k) *= s(i) * s(k);
    }
  }

  // compute the contributions of gradients for the fourth term
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      elem_mats[3](i, k) *=
          (grads.col(index(i)).transpose() * grads.col(index(k)));
      elem_mats[3](i, k) *= s(i) * s(k);
    }
  }

  // sum up all the four contributions
  Eigen::Matrix3d elem_mat = Eigen::MatrixXd::Zero(3, 3);
  for (int i = 0; i < 4; i++) {
    elem_mat += elem_mats[i];
  }

  return elem_mat;
}

}  // namespace projects::hldo_sphere::assemble
