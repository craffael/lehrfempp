/** @file geometry_output.cc
 * Demo for outputting information about geometry
 */

#include <iostream>
#include "lf/geometry/geometry.h"
#include "lf/geometry/point.h"
#include "lf/geometry/print_info.h"

int main() {

  using namespace lf::geometry;
  lf::geometry::Geometry::output_ctrl_ = 100;
  std::cout << "Output of information of geometry" << std::endl;


  // Create point
  lf::geometry::Point test_point(Eigen::Vector3d(1,0,0));

  // Print point
  PrintInfo(test_point, std::cout);

  std::cout << test_point << std::endl;

  return 0L;
}
