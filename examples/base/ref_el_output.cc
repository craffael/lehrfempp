/** @file ref_el_output.cc
 * Demo for outputting information about reference elements
 */

#include <iostream>
#include "lf/base/base.h"

int main(int argc, char** argv) {
  // Node
  auto re_node = lf::base::RefEl::kPoint();
  std::cout << "Information about node: \n" << re_node << std::endl;

  // Segment
  auto re_seg = lf::base::RefEl::kSegment();
  std::cout << "\nInformation about segment: \n" << re_seg << std::endl;

  // Triangle
  auto re_tria = lf::base::RefEl::kTria();
  std::cout << "\nInformation about triangle: \n" << re_tria << std::endl;

  // Qaud
  auto re_quad = lf::base::RefEl::kQuad();
  std::cout << "\nInformation about quad: \n" << re_quad << std::endl;

  return 0L;
}
