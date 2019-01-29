/**
 * @file
 * @brief Driver function for simple LehrFEM++ demo
 * @author Ralf Hiptmair
 * @date   January 2019
 * @copyright MIT License
 */

#include <boost/program_options.hpp>
#include "lecturedemoassemble.h"
#include "lecturedemodof.h"
#include "lecturedemomesh.h"
#include "lecturedemoquad.h"

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()
  ("help,h", "This message")
  ("demo_number,d", po::value<int>()->default_value(0), "Selector for demo code");
  // clang-format on
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  if (vm.count("help") > 0) {
    std::cout << desc << std::endl;
    std::cout << "No arg: run all demos" << std::endl;
    std::cout << "N = 1: demo of LehrFEM++ mesh capabilities" << std::endl;
    std::cout << "N = 2: demo of LehrFEM++ DofHandler capabilities"
              << std::endl;
    std::cout << "N = 3: demo of LehrFEM++ assembly of LSE" << std::endl;
    std::cout << "N = 4: demo of numerical quadrature in LehrFEM++"
              << std::endl;
    std::cout << "N = 5: demo of solving a Dirichlet BVP" << std::endl;
  } else {
    int selector = vm["demo_number"].as<int>();
    if ((selector == 1) || (selector == 0)) {
      lecturedemo::lecturedemomesh();
    }
    if ((selector == 2) || (selector == 0)) {
      lecturedemo::lecturedemodof();
    }
    if ((selector == 3) || (selector == 0)) {
      lecturedemo::lecturedemoassemble();
    }
    if ((selector == 4) || (selector == 0)) {
      lecturedemo::lecturedemoquad();
    }
    if ((selector == 5) || (selector == 0)) {
      lecturedemo::lecturedemoDirichlet();
    }
  }

  return 0L;
}
