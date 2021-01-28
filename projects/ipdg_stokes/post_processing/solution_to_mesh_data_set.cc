
#include "solution_to_mesh_data_set.h"

#include <lf/uscalfe/lagr_fe.h>

namespace projects::ipdg_stokes::post_processing {

lf::mesh::utils::CodimMeshDataSet<double> extractBasisFunctionCoefficients(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh, const Eigen::VectorXd &solution) {
  lf::mesh::utils::CodimMeshDataSet<double> coefficients(mesh, 2);
  for (const auto *const node : mesh->Entities(2)) {
    coefficients(*node) = solution[dofh.GlobalDofIndices(*node)[0]];
  }
  return coefficients;
}

lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> extractVelocity(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const lf::assemble::DofHandler &dofh, const Eigen::VectorXd &solution) {
  lf::mesh::utils::CodimMeshDataSet<Eigen::Vector2d> velocity(mesh, 0);
  for (const auto *const cell : mesh->Entities(0)) {
    const auto *geom = cell->Geometry();
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
    const auto idx = dofh.GlobalDofIndices(*cell);
    velocity(*cell) = solution[idx[0]] * basis_funct.col(0) +
                      solution[idx[1]] * basis_funct.col(1) +
                      solution[idx[2]] * basis_funct.col(2);
  }
  return velocity;
}

}  // end namespace projects::ipdg_stokes::post_processing
