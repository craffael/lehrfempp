/**
 * @file
 * @brief Outputs location of global shape functions as managed by a Dofhandler
 * @author Ralf Hiptmair
 * @date   October 2018
 * @copyright MIT License
 */

#include <boost/program_options.hpp>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2dp/hybrid2dp.h>

int main(int argc, char** argv) {
  // define allowed command line arguments:
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()
  ("help", "-ndof_node <N> -ndof_edge <N> -ndof_tria <N> -ndof_quad <N>")
  ("ndof_node", po::value<int>()->default_value(1), "No of dofs on nodes")
  ("ndof_edge", po::value<int>()->default_value(2), "No of dofs on edges")
  ("ndof_tria", po::value<int>()->default_value(1), "No of dofs on triangles")
  ("ndof_quad", po::value<int>()->default_value(4), "Mp of dofs on quadrilaterals")
  ;
  // clang-format on
}
