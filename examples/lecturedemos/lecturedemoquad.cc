/**
 * @file
 * @brief Simple LehrFEM++ demo and sample codes for the use of numerical
 * quadrature; codes for lecture document
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include "lecturedemoquad.h"

namespace lecturedemo {

void lecturedemoquad() {
  std::cout << "LehrFEM++ demo: composite quadrature on mesh" << std::endl;
  // Obtain an affine hybrid mesh from the collection of LehrFEM++'s
  // built-in meshes
  std::shared_ptr<lf::mesh::Mesh> mesh_p{
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(5)};
  /* SAM_LISTING_BEGIN_1 */
  // Function to be integrated
  auto f = [](const Eigen::VectorXd& x) -> double {
    return (x[0] * x[0] + x[1] * x[1]);
  };
  // Cell-based composite quadrature covering all cells
  double integral_val = localQuadFunction(
      *mesh_p,
      {{lf::base::RefEl::kTria(), lf::quad::make_TriaQR_EdgeMidpointRule()},
       {lf::base::RefEl::kQuad(), lf::quad::make_QuadQR_P4O4()}},
      f, 0);

  /* SAM_LISTING_END_1 */
  std::cout << "Integral value = " << integral_val << std::endl;
}

}  // namespace lecturedemo
