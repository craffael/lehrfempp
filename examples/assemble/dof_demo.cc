/**
 * @file
 * @brief Outputs location of global shape functions as managed by a Dofhandler
 * @author Ralf Hiptmair
 * @date   October 2018
 * @copyright MIT License
 */

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>

#include <boost/program_options.hpp>

#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

int main(int argc, char **argv) {
  // The following code is modeled after the example from
  // https://theboostcpplibraries.com/boost.program_options
  // and defines allowed command line arguments:
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  // clang-format off
  desc.add_options()
  ("help,h", "--ndof_node <N> --ndof_edge <N> --ndof_tria <N> --ndof_quad <N>")
  ("ndof_node,n", po::value<int>()->default_value(1), "No of dofs on nodes")
  ("ndof_edge,e", po::value<int>()->default_value(2), "No of dofs on edges")
  ("ndof_tria,t", po::value<int>()->default_value(1), "No of dofs on triangles")
  ("ndof_quad,q", po::value<int>()->default_value(4), "Mp of dofs on quadrilaterals");
  // clang-format on
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  if (vm.count("help") > 0) {
    std::cout << desc << '\n';
  } else {
    // Retrieve number of degrees of freedom for each entity type from command
    // line arguments
    lf::base::size_type ndof_node = 1;
    if (vm.count("ndof_node") > 0) {
      ndof_node = vm["ndof_node"].as<int>();
    }
    lf::base::size_type ndof_edge = 2;
    if (vm.count("ndof_edge") > 0) {
      ndof_edge = vm["ndof_edge"].as<int>();
    }
    lf::base::size_type ndof_tria = 1;
    if (vm.count("ndof_tria") > 0) {
      ndof_tria = vm["ndof_tria"].as<int>();
    }
    lf::base::size_type ndof_quad = 4;
    if (vm.count("ndof_quad") > 0) {
      ndof_quad = vm["ndof_quad"].as<int>();
    }
    std::cout << "LehrFEM++ demo: assignment of global shape functions"
              << '\n';
    std::cout << "#dof/vertex = " << ndof_node << '\n';
    std::cout << "#dof/edge = " << ndof_edge << '\n';
    std::cout << "#dof/triangle = " << ndof_tria << '\n';
    std::cout << "#dof/quadrilateral = " << ndof_quad << '\n';

    // Build a mesh comprising two cells
    const std::shared_ptr<lf::mesh::Mesh> mesh_p =
        lf::mesh::test_utils::GenerateHybrid2DTestMesh(2);
    // Output information about the mesh
    lf::mesh::utils::PrintInfo(std::cout, *mesh_p);

    // Create a dof handler object describing a uniform distribution
    // of shape functions
    const lf::assemble::UniformFEDofHandler dof_handler(
        mesh_p, {{lf::base::RefEl::kPoint(), ndof_node},
                 {lf::base::RefEl::kSegment(), ndof_edge},
                 {lf::base::RefEl::kTria(), ndof_tria},
                 {lf::base::RefEl::kQuad(), ndof_quad}});
    // Copious output of information about dof handler
    PrintInfo(std::cout, dof_handler, 30);
  }
  return 0L;
}  // end main
