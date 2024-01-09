/** @file refinement_demo.cc
 * @brief example code for building and refining a mesh
 *
 * This code generates a test mesh and performs a user specified number of
 * steps of regular uniform refinement.
 *
 * In the end the information about all meshes created in the process
 * is stored in the form of MATLAB functions.
 */

#include <lf/io/io.h>
#include <lf/refinement/mesh_hierarchy.h>
#include <lf/refinement/refutils.h>

#include <boost/program_options.hpp>
#include <iostream>
#include <sstream>

#include "lf/base/base.h"
#include "lf/mesh/hybrid2d/hybrid2d.h"
#include "lf/mesh/test_utils/check_mesh_completeness.h"
#include "lf/mesh/test_utils/test_meshes.h"
#include "lf/mesh/utils/utils.h"

// Global variable for selection of refinement type
// 0 -> regular global refinement
// 1 -> global baryccentric refinement
// 2 -> global marking of edges
// 3 -> local marking of edges with midpoints in [0,1]^2

int main(int argc, char **argv) {
  using size_type = lf::base::size_type;
  using lf::io::TikzOutputCtrl;
  namespace po = boost::program_options;

  // setup boost program options:
  po::options_description desc(
      "LehrFEM++ demo of mesh construction and refinement");

  // clang-format off
  desc.add_options()
    ("help,h", "Display this help message")
    ("refselector", po::value<int>()->default_value(0), "Selector for refinement method (0 [default] to 3)")
    ("refsteps", po::value<int>()->default_value(3), "Number of refinement steps")
    ("basename", po::value<std::string>()->default_value("out"), "Matlab basename for the output files")
  ;
  // clang-format on

  po::variables_map vm;
  try {
    po::store(po::parse_config_file<char>("setup.vars", desc), vm);
  } catch (const po::reading_file &) {
    std::cout << "No file `setup.vars` specifying control variables" << '\n';
  }
  po::store(po::parse_command_line(argc, argv, desc), vm);

  // check for the help option (-h or --help)
  if (vm.count("help") > 0) {
    std::cout << desc << '\n';
    return 1;
  }

  // get the value for refselector
  auto refselector = vm["refselector"].as<int>();
  auto refsteps = vm["refsteps"].as<int>();

  // Generate hybrid test mesh and obtain a pointer to it
  const std::shared_ptr<lf::mesh::Mesh> mesh_ptr =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh();
  const lf::mesh::Mesh &mesh = *mesh_ptr;

  // Output information about the mesh
  std::cout << "########################################" << '\n';
  std::cout << "Output of PrintInfo() for test mesh" << '\n';
  lf::mesh::utils::PrintInfo(std::cout, mesh);

  // Build mesh hierarchy
  std::cout << "########################################" << '\n';
  std::cout << "Initialization of data structure for mesh hierarchy" << '\n';
  std::unique_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::refinement::MeshHierarchy multi_mesh(mesh_ptr,
                                           std::move(mesh_factory_ptr));

  // lambda functions for (local) marking
  // First, a global marker
  const std::function<bool(const lf::mesh::Mesh &, const lf::mesh::Entity &)>
      allmarker = [](const lf::mesh::Mesh & /*mesh*/, const lf::mesh::Entity &
                     /*edge*/) -> bool { return true; };
  // Mark edges whose center lies inside a square
  const std::function<bool(const lf::mesh::Mesh &,
                           const lf::mesh::Entity &edge)>
      locmarker = [](const lf::mesh::Mesh & /*mesh*/,
                     const lf::mesh::Entity &edge) -> bool {
    Eigen::MatrixXd ref_c(1, 1);
    ref_c(0, 0) = 0.5;
    Eigen::VectorXd c(edge.Geometry()->Global(ref_c));
    return ((c[0] >= 0.0) && (c[0] < 1.0) && (c[1] >= 0.0) && (c[1] < 1.0));
  };

  // Main refinement loop
  std::cout << "########################################" << '\n';

  std::cout << "No of refinement steps: " << refsteps;
  std::cout << "Entering main refinement loop" << '\n';
  for (std::size_t refstep = 0; refstep < refsteps; refstep++) {
    std::cout << "#### Refinement step " << refstep + 1 << '\n';
    // Depending on the value of the control variable do different types
    // of refinement
    switch (refselector) {
      case 0: {
        // Global regular refinement
        std::cout << "#### global regular refinement" << '\n';
        multi_mesh.RefineRegular();
        break;
      }
      case 1: {
        // Global barycentric refinement
        std::cout << "#### global barycentric refinement" << '\n';
        multi_mesh.RefineRegular(lf::refinement::RefPat::rp_barycentric);
        break;
      }
      case 2: {
        // Global regular refinement
        std::cout << "#### global bisection refinement" << '\n';
        multi_mesh.MarkEdges(allmarker);
        multi_mesh.RefineMarked();
        break;
      }
      case 3: {
        // Local refinement
        std::cout << "#### local bisection refinement" << '\n';
        multi_mesh.MarkEdges(locmarker);
        multi_mesh.RefineMarked();
        break;
      }
      default: {
        LF_VERIFY_MSG(false, "Unknown refinement selector");
        break;
      }
    }  // end switch refselector

    // Inquire about number of existing levels
    const size_type n_levels = multi_mesh.NumLevels();
    // Obtain pointer to mesh on finest level
    const std::shared_ptr<const lf::mesh::Mesh> mesh =
        multi_mesh.getMesh(n_levels - 1);
    // Print number of entities of various co-dimensions
    std::cout << "#### Mesh on level " << n_levels - 1 << ": "
              << mesh->NumEntities(2) << " nodes, " << mesh->NumEntities(1)
              << " edges, " << mesh->NumEntities(0) << " cells," << '\n';
    std::stringstream level_asc;
    level_asc << refstep;

    writeTikZ(*mesh, std::string("refinement_mesh") + level_asc.str() + ".txt",
              TikzOutputCtrl::RenderCells | TikzOutputCtrl::CellNumbering |
                  TikzOutputCtrl::VerticeNumbering |
                  TikzOutputCtrl::NodeNumbering |
                  TikzOutputCtrl::EdgeNumbering);

    const lf::io::VtkWriter vtk_writer(
        mesh, "refinement_mesh" + level_asc.str() + ".vtk");
  }

  // Generate  MATLAB functions that provide a description of all
  // levels of the mesh hierarchy
  auto basename = vm["basename"].as<std::string>();
  lf::refinement::WriteMatlab(multi_mesh, basename);

  return 0;
}
