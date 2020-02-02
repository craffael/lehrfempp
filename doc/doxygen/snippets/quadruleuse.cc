/**
 * @file
 * @brief Doxygen snippets to show usage of LehrFEM++'s quadrature support
 * @author Ralf Hiptmair
 * @date Jan 28, 2020
 * @copyright MIT License
 */

#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/quad/quad.h>
#include <iostream>

//! [quadrule]
template <typename FUNCTOR>
auto localQuadFunction(const lf::mesh::Mesh &mesh,
                       const lf::quad::QuadRule &quadrule, FUNCTOR &&f) {
  // Variable for summing the result
  using value_t = std::invoke_result_t<FUNCTOR, Eigen::VectorXd>;
  value_t sum_var{};
  // Loop over entities of co-dimension 0 = cells
  for (const lf::mesh::Entity *entity : mesh.Entities(0)) {
    // Check matching of reference element
    LF_VERIFY_MSG(entity->RefEl() == quadrule.RefEl(),
                  "Mismatch of reference element for " << entity->RefEl());
    // Obtain geometry information for entity
    const lf::geometry::Geometry &geo{*entity->Geometry()};
    // Number of quadrature points
    const int P = quadrule.NumPoints();
    // Quadrature points
    const Eigen::MatrixXd zeta_ref{quadrule.Points()};
    // Map quadrature points to physical/world coordinates
    const Eigen::MatrixXd zeta{geo.Global(zeta_ref)};
    // Quadrature weights
    const Eigen::VectorXd w_ref{quadrule.Weights()};
    // Gramian determinants
    const Eigen::VectorXd gram_dets{geo.IntegrationElement(zeta_ref)};
    // Iterate over the quadrature points
    for (int l = 0; l < P; ++l) {
      sum_var += w_ref[l] * f(zeta.col(l)) * gram_dets[l];
    }
  }
  return sum_var;
}
//! [quadrule]

void doquad() {
  // Obtain a triangular mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p{
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3)};
  // Cell-based composite quadrature covering all cells
  double integral_val =
      localQuadFunction(*mesh_p, lf::quad::make_TriaQR_EdgeMidpointRule(),
                        [](const Eigen::VectorXd &x) -> double {
                          return (x[0] * x[0] + x[1] * x[1]);
                        });
  std::cout << "Integral value = " << integral_val << std::endl;
}
