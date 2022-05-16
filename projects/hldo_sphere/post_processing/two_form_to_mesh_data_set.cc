#include "two_form_to_mesh_data_set.h"

namespace projects::hldo_sphere::post_processing {

lf::mesh::utils::CodimMeshDataSet<double> extractSolution(
    const std::shared_ptr<const lf::mesh::Mesh> &mesh,
    const Eigen::VectorXd &solution) {
  lf::mesh::utils::CodimMeshDataSet<double> coefficients(mesh, 0);
  for (const lf::mesh::Entity *const cell : mesh->Entities(0)) {
    coefficients(*cell) = solution[mesh->Index(*cell)];
  }
  return coefficients;
}

}  // namespace projects::hldo_sphere::post_processing
