/** @file ref_el_output.cc
 * Demo for outputting information about reference elements
 */

#include <iostream>

#include "lf/base/base.h"

int main() {
  // Node
  auto re_node = lf::base::RefEl::kPoint();
  std::cout << "Information about node: \n" << re_node << '\n';

  // Segment
  auto re_seg = lf::base::RefEl::kSegment();
  std::cout << "\nInformation about segment: \n" << re_seg << '\n';

  // Triangle
  auto re_tria = lf::base::RefEl::kTria();
  std::cout << "\nInformation about triangle: \n" << re_tria << '\n';

  // Qaud
  auto re_quad = lf::base::RefEl::kQuad();
  std::cout << "\nInformation about quad: \n" << re_quad << '\n';

  return 0L;
}
