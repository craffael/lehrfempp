/**
 * @file
 * @brief Outputs location of global shape functions as managed by a Dofhandler
 * @author Ralf Hiptmair
 * @date   October 2018
 * @copyright MIT License
 */

#include <lf/base/base.h>
namespace ci = lf::base::ci; // avoid typing lf::base all the time

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include "lf/assemble/assemble.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

int main(int argc, char** argv) {
  ci::Init(argc, argv); // Initialise input reader with command line arguments
  // Add options
  ci::Add("help,h", "--ndof_node <N> --ndof_edge <N> --ndof_tria <N> --ndof_quad <N>");
  ci::Add<int>("ndof_node,n", "No of dofs on nodes", 1); // Default value 2
  ci::Add<int>("ndof_edge,e", "No of dofs on edges", 2);
  ci::Add<int>("ndof_tria,t", "No of dofs on triangles", 1);
  ci::Add<int>("ndof_quad,q", "Mp of dofs on quadrilaterals", 4);
  // Alternatively:
  //namespace po = boost::program_options;
  //ci::Add()
  //("help,h", "--ndof_node <N> --ndof_edge <N> --ndof_tria <N> --ndof_quad <N>")
  //("ndof_node,n", po::value<int>()->default_value(1), "No of dofs on nodes")
  //("ndof_edge,e", po::value<int>()->default_value(2), "No of dofs on edges")
  //("ndof_tria,t", po::value<int>()->default_value(1), "No of dofs on triangles")
  //("ndof_quad,q", po::value<int>()->default_value(4), "Mp of dofs on quadrilaterals");
  // clang-format on
  ci::ParseCommandLine();
  if (ci::Help()) {
    return 0;
  } else {
    // Retrieve number of degrees of freedom for each entity type from command
    // line arguments
    // Instead of:
    //lf::base::size_type ndof_node = 1;
    //if (vm.count("ndof_node") > 0) {
    //  ndof_node = vm["ndof_node"].as<int>();
    //}
    // Simply use the following. Note that the default value has already been set.
    lf::base::size_type ndof_node = ci::Get<int>("ndof_node"); // does this int to size_type cast work??
    // If we wouldn't have set a default value we could call it like this:
    //ndof_node = ci::Get<lf::base::size_type>("ndof_node", 1);
    lf::base::size_type ndof_edge = ci::Get<int>("ndof_edge");
    lf::base::size_type ndof_tria = ci::Get<int>("ndof_tria");
    lf::base::size_type ndof_quad = ci::Get<int>("ndof_quad");
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
    lf::mesh::utils::printinfo_ctrl = 100;
    lf::mesh::Entity::output_ctrl_ = 0;
    lf::mesh::utils::PrintInfo(*mesh_p, std::cout);

    // Create a dof handler object describing a uniform distribution
    // of shape functions
    lf::assemble::UniformFEDofHandler dof_handler(
        mesh_p, {{lf::base::RefEl::kPoint(), ndof_node},
                 {lf::base::RefEl::kSegment(), ndof_edge},
                 {lf::base::RefEl::kTria(), ndof_tria},
                 {lf::base::RefEl::kQuad(), ndof_quad}});
    lf::assemble::DofHandler::output_ctrl_ = 30;
    std::cout << dof_handler << std::endl;
  }
}  // end main
