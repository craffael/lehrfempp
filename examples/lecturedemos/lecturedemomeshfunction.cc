/**
 * @file
 * @brief Simple demonstrations of mesh refinement
 * @author Ralf Hiptmair
 * @date   March 2019
 * @copyright MIT License
 */

#include "lecturedemomeshfunction.h"

namespace lecturedemo {

void lecturedemomeshfunction() {
  std::cout << "LehrFEM++ mesh function demos" << std::endl;
  // (I) Prepare meshes
  // Obtain hybrid test mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(0);
  // Create sequence of meshes by regular refinement
  const unsigned int refsteps = 5;
  std::shared_ptr<lf::refinement::MeshHierarchy> multi_mesh_p =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh_p, refsteps);
  lf::refinement::MeshHierarchy& multi_mesh{*multi_mesh_p};

  // A function for testing
  auto f = [](Eigen::Vector2d x) -> double {
    return std::sin(x[0] * x[0] - x[1]);
  };
  // A vectorfield for testing
  auto vf = [](Eigen::Vector2d x) {
    return Eigen::Vector2d(-x[1], x[0]);
  };
  // A coefficient function, matrix-valued
  auto cf = [](Eigen::Vector2d x) -> Eigen::Matrix2d {
    return (Eigen::Matrix2d() << x[0] * x[0], x[0] * x[1], x[0] * x[1],
            x[1] * x[1])
        .finished();
  };
  // Number of levels
  lf::base::size_type L = multi_mesh.NumLevels();

  // Loop over all meshes
  for (int level = 0; level < L; ++level) {
    std::shared_ptr<const lf::mesh::Mesh> lev_mesh_p =
        multi_mesh.getMesh(level);
    // Build FE space
    std::shared_ptr<const lf::uscalfe::FeSpaceLagrangeO2<double>> fe_space_p =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO2<double>>(lev_mesh_p);
    // Compute FE nodal interpolant
    lf::mesh::utils::MeshFunctionGlobal mf_f(f);
    auto coeff_vec{lf::fe::NodalProjection(*fe_space_p, mf_f)};
    std::cout << "Level " << level << ", integral = "
              << integrateCoeffgradUhVf(fe_space_p, coeff_vec, cf, vf)
              << std::endl;
  }
}

}  // namespace lecturedemo
