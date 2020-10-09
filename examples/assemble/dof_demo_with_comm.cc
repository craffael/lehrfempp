/**
 * @file
 * @brief Outputs location of global shape functions as managed by a Dofhandler
 * @author Ralf Hiptmair
 * @date   October 2018
 * @copyright MIT License
 */

#include <lf/base/base.h>

#include <boost/program_options.hpp>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include "lf/assemble/assemble.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

int main(int argc, char** argv) {
  // Add options
  namespace po = boost::program_options;

  po::options_description desc;

  // clang-format off
  desc.add_options()
    ("help,h", "Display this help message")
    ("ndof_node,n", po::value<int>()->default_value(1), "No of dofs on nodes")
    ("ndof_edge,e", po::value<int>()->default_value(2), "No of dofs on edges")
    ("ndof_tria,t", po::value<int>()->default_value(1), "No of dofs on triangles")
    ("ndof_quad,q", po::value<int>()->default_value(4), "No of dofs on quadrilaterals");
  // clang-format on

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }
  vm.notify();

  // Retrieve number of degrees of freedom for each entity type from command
  // line arguments
  lf::base::size_type ndof_node = vm["ndof_node"].as<int>();
  lf::base::size_type ndof_edge = vm["ndof_edge"].as<int>();
  lf::base::size_type ndof_tria = vm["ndof_tria"].as<int>();
  lf::base::size_type ndof_quad = vm["ndof_quad"].as<int>();
  std::cout << "LehrFEM++ demo: assignment of global shape functions"
            << std::endl;
  std::cout << "#dof/vertex = " << ndof_node << std::endl;
  std::cout << "#dof/edge = " << ndof_edge << std::endl;
  std::cout << "#dof/triangle = " << ndof_tria << std::endl;
  std::cout << "#dof/quadrilateral = " << ndof_quad << std::endl;

  // Build a mesh comprising two cells
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(2);
  // Output information about the mesh
  lf::mesh::utils::PrintInfo(std::cout, *mesh_p);

  // Create a dof handler object describing a uniform distribution
  // of shape functions
  lf::assemble::UniformFEDofHandler dof_handler(
      mesh_p, {{lf::base::RefEl::kPoint(), ndof_node},
               {lf::base::RefEl::kSegment(), ndof_edge},
               {lf::base::RefEl::kTria(), ndof_tria},
               {lf::base::RefEl::kQuad(), ndof_quad}});
  lf::assemble::DofHandler::output_ctrl_ = 30;
  std::cout << dof_handler << std::endl;

  return 0;
}  // end main
