/** @file ref_el_output.cc
 * Demo for outputting information about reference elements
 */

#include <iostream>
#include "lf/base/base.h"
#include "lf/base/comm.h"
namespace ci = lf::base::ci;

int main(int argc, char** argv) {
  std::cout << "Output of information on reference elements" << std::endl;
  ci::Init(argc, argv);
  ci::Add()("help,h", "Print this help message.");
  ci::ParseCommandLine();

  if (ci::Help()) {
    std::cout << "Try for instance: --RefEl_ctrl 2\n";
    return 0L;
  }

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
