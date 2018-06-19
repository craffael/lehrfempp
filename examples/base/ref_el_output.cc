/** @file ref_el_output.cc
 * Demo for outputting information about reference elements 
 */

#include <iostream>
#include "lf/base/base.h"

int main(int argc,char **argv) {
  std::cout << "Output of information on reference elements" << std::endl;

  auto re_node = lf::base::RefEl::kPoint();
  std::cout << "Information about node: " << re_node << std::endl;

  auto re_seg = lf::base::RefEl::kSegment();
  std::cout << "Information about segment: " << re_seg << std::endl;
							
  auto re_tria = lf::base::RefEl::kTria();
  std::cout << "Information about triangle: " << re_tria << std::endl;
							
  auto re_quad = lf::base::RefEl::kQuad();
  std::cout << "Information about quad: " << re_quad << std::endl;
  return 0L; 
}
